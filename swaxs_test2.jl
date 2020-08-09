using DelimitedFiles, Distributed, Dates

cd("G:\\My Drive\\18. Github Repo\\SWAXS.jl");
NPROCS = Sys.CPU_THREADS;
nprocs() > NPROCS ? rmprocs(workers()[end - (nprocs() - NPROCS) + 1: end]) : addprocs(NPROCS - nprocs());

@everywhere include(".\\src\\SWAXS.jl");
@everywhere using .SWAXS


function solvent_batch(dir::AbstractString, solutefn::AbstractString, sourceprefix::AbstractString; solventprefix::AbstractString="solvent", d::Float64=10)

    filelist = readdir(dir);
    filelist = filelist[endswith.(filelist, ".pdb")];

    sourcelist = filelist[startswith.(filelist, sourceprefix)];
    nsource = length(sourcelist);
    solute = joinpath(dir, solutefn);

    for i = 1:nsource
        s = sourcelist[i];
        @info("Making solvent using $s ... ");
        outputfn = joinpath(dir, replace(s, sourceprefix => solventprefix));
        isfile(outputfn) ? @warn("$outputfn exists, skipped ... ") : make_solvent(solute, joinpath(dir, s), d, outputfn);
    end

    return nothing;
end


dir = "C:\\Users\\Yen-Lin\\Desktop\\solvent";

# make solvent
solvent_batch(dir, "solute.pdb", "frame"; solventprefix="solvent", d=2.0);

# compute the buffer-subttracted
filelist = readdir(dir);
solventlist = joinpath.(dir, filelist[startswith.(filelist, "solvent")]);
q = collect(0.0:0.025:1.5);

@info("Computing SWAXS ... ");
@info("Starting Time: $(now())");
@time intensity = PDBSWAXS(joinpath(dir, "solute.pdb"), solventlist, q; J=1500);
writedlm(joinpath(dir, "profile_d2.dat"), [q intensity]);
