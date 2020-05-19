# extract the solvent only
function make_solvent(solutefn::String, outputfn::String)

    soluteio = open(solutefn, "r");
    solutelines = readlines(soluteio);
    close(soluteio);
    solutelines = solutelines[startswith.(solutelines, "ATOM")];

    IONS = ["K"; "CL"; "MGH"; "NA"; "Na+"; "Cl-"];
    SOL = ["SOL", "WAT", "HOH"];

    outputlines = Array{String, 1}(undef, 0);

    # get the waters from the solute frame
    for i = 1:length(solutelines)
        substr = split(solutelines[i]);
        substr[4] in vcat(SOL, IONS) ? append!(outputlines, solutelines[i:i]) : nothing;
    end

    fout = open(outputfn, "w");
    for i = 1:length(outputlines)
        write(fout, outputlines[i]);
        write(fout, "\n");
    end
    close(fout);

    return nothing;
end



# sticking the ions
function make_solvent(solutefn::String, sourcefn::String, outputfn::String)

    soluteio = open(solutefn, "r");
    solutelines = readlines(soluteio);
    close(soluteio);
    solutelines = solutelines[startswith.(solutelines, "ATOM")];

    sourceio = open(sourcefn, "r");
    sourcelines = readlines(sourceio);
    close(sourceio);
    sourcelines = sourcelines[startswith.(sourcelines, "ATOM")];

    IONS = ["K"; "CL"; "MGH"; "NA"];
    SOL = ["SOL", "WAT", "HOH"];

    outputlines = Array{String, 1}(undef, 0);

    # get the waters from the solute frame
    for i = 1:length(solutelines)
        substr = split(solutelines[i]);
        substr[4] in SOL ? append!(outputlines, solutelines[i:i]) : nothing;
    end

    # add the ions from the source
    for i = 1:length(sourcelines)
        substr = split(sourcelines[i]);
        substr[4] in IONS ? append!(outputlines, sourcelines[i:i]) : nothing;
    end

    fout = open(outputfn, "w");
    for i = 1:length(outputlines)
        write(fout, outputlines[i]);
        write(fout, "\n");
    end
    close(fout);

    return nothing;
end




# sticking everything (water and ions) within the solute envelope
function make_solvent(solutefn::String, sourcefn::String, ξ::Float64, outputfn::String)

    @info("--- Making Solvent: Processing solute and bulk ...");
    soluteio = open(solutefn, "r");
    solutelines = readlines(soluteio);
    close(soluteio);
    solutelines = solutelines[startswith.(solutelines, "ATOM")];
    __, solutepos = SimplyPDB(solutefn; waters=true, ions=true);

    sourceio = open(sourcefn, "r");
    sourcelines = readlines(sourceio);
    close(sourceio);
    sourcelines = sourcelines[startswith.(sourcelines, "ATOM")];
    __, sourcepos = SimplyPDB(sourcefn);

    IONS = ["K"; "CL"; "MGH"; "NA"];
    SOL = ["SOL"; "WAT"; "HOH"];

    # count solute atoms
    ndist = 0;
    for i = 1:length(solutelines)
        substr = split(solutelines[i]);
        substr[4] in [SOL; IONS] ? nothing : ndist += 1;
    end

    # Construct the RNA gaussians
    Σ = Matrix{Float64}(I, 3, 3);
    dist = Array{MvNormal, 1}(undef, ndist);
    for i = 1:ndist
        dist[i] = MvNormal(solutepos[i, :], Σ);
    end

    # define cutoff
    cutoff = 10.0 * pdf(MvNormal(zeros(3), Σ), [ξ; 0.0; 0.0]);
    @info("--- Making Solvent: cutoff = $cutoff ...");

    outputlines = Array{String, 1}(undef, 0);
    # get the waters from the solute frame
    for i = (ndist+1):length(solutelines)
        r = solutepos[i, :];
        sum(map(x -> pdf(x, r), dist)) <= cutoff ? append!(outputlines, solutelines[i:i]) : nothing;
    end

    # i = ndist + 1;
    # while i <= length(solutelines)
    #     substr = split(solutelines[i]);
    #     r = solutepos[i, :];
    #     if sum(map(x -> pdf(x, r), dist)) <= cutoff
    #         if substr[4] in SOL
    #             append!(outputlines, solutelines[i:i+2]);
    #             i += 3;
    #         else
    #             append!(outputlines, [solutelines[i]]);
    #             i += 1;
    #         end
    #     else
    #         i += 1;
    #     end
    # end


    # fill in the blank of the molecules.
    for i = 1:length(sourcelines)
        r = sourcepos[i, :];
        sum(map(x -> pdf(x, r), dist)) >= cutoff ? append!(outputlines, sourcelines[i:i]) : nothing;
    end
    # i = 1;
    # while i <= length(sourcelines)
    #     substr = split(sourcelines[i]);
    #     r = sourcepos[i, :];
    #     if sum(map(x -> pdf(x, r), dist)) >= cutoff
    #         if substr[4] in SOL
    #             append!(outputlines, sourcelines[i:i+2]);
    #             i += 3;
    #         else
    #             append!(outputlines, sourcelines[i:i]);
    #             i += 1;
    #         end
    #     else
    #         i += 1;
    #     end
    # end


    fout = open(outputfn, "w");
    for i = 1:length(outputlines)
        write(fout, outputlines[i]);
        write(fout, "\n");
    end
    close(fout);

    return nothing;
end
