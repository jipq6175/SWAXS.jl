
# batch swaxs

using DelimitedFiles, Distributed
cd("G:\\My Drive\\18. Github Repo\\SWAXS.jl");
NPROCS = Sys.CPU_THREADS;
nprocs() > NPROCS ? rmprocs(workers()[end - (nprocs() - NPROCS) + 1: end]) : addprocs(NPROCS - nprocs());

@everywhere include(".\\src\\SWAXS.jl");
@everywhere using .SWAXS



function solvent_batch(dir::AbstractString, soluteprefix::AbstractString, sourcefn::AbstractString; solventprefix::AbstractString="solvent")

    filelist = readdir(dir);
    filelist = filelist[endswith.(filelist, ".pdb")];

    solutelist = filelist[startswith.(filelist, soluteprefix)];
    nsolute = length(solutelist);
    source = joinpath(dir, sourcefn);

    for i = 1:nsolute
        @info("Making solvent for $(solutelist[i]) ... ");
        outputfn = joinpath(dir, solventprefix*"$i.pdb");
        make_solvent(joinpath(dir, solutelist[i]), source, outputfn);
    end

    return nothing;
end






function swaxs_batch(dir::AbstractString, q::AbstractVector; J::Signed=1500, soluteprefix::AbstractString="frame", solventprefix::AbstractString="bulk")

    filelist = readdir(dir);
    filelist = filelist[endswith.(filelist, ".pdb")];

    solutelist = filelist[startswith.(filelist, soluteprefix)];
    solventlist = filelist[startswith.(filelist, solventprefix)];

    nsolute = length(solutelist);
    nsolvent = length(solventlist);

    if nsolvent == 1
        @info("Using one single solvent frame for all the solutes ...");
        solventlist = fill(solventlist[1], nsolute);
    elseif nsolvent == 0
        @warn("No solvent found ... ");
        return nothing;
    elseif nsolvent != nsolute
        @warn("Solvent and solute are not matching ... ");
        return nothing;
    else
        @info("Solute and solvent match ...");
    end

    @info("Starting to compute the SWAXS profiles of $nsolute frames .....");
    for i = 1:nsolute


        datfile = joinpath(dir, solutelist[i][1:end-4] * ".dat");
        if isfile(datfile)
            @warn("The SWAXS profile for $(solutelist[i]) exists, skipping this frame ...");
        else
            t = @elapsed mat = PDBSWAXS(joinpath(dir, solutelist[i]), joinpath(dir, solventlist[i]), q; J=J);
            writedlm(datfile, [q mat]);
            @info("The SWAXS profile of the $(solutelist[i]) was computed using solvent $(solventlist[i]) with $t seconds ...");
        end

    end

    return nothing;
end


dir = "C:\\Users\\Yen-Lin\\Desktop\\YenAATT_120mM_NA";


# Batch solvent
@time solvent_batch(dir, "frame", "bulk.pdb");















#q = collect(0.0:0.005:1.5);
#swaxs_batch(dir, q; J=1500, soluteprefix="frame", solventprefix="solvent");
