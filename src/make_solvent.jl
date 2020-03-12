

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
        substr[4] in SOL ? append!(outputlines, [solutelines[i]]) : nothing;
    end

    # add the ions from the source
    for i = 1:length(sourcelines)
        substr = split(sourcelines[i]);
        substr[4] in IONS ? append!(outputlines, [sourcelines[i]]) : nothing;
    end

    fout = open(outputfn, "w");
    for i = 1:length(outputlines)
        write(fout, outputlines[i]);
        write(fout, "\n");
    end
    close(fout);

    return nothing;
end
