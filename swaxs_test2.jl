using DelimitedFiles, Distributed, Dates

cd("G:\\My Drive\\18. Github Repo\\SWAXS.jl");
NPROCS = Sys.CPU_THREADS;
nprocs() > NPROCS ? rmprocs(workers()[end - (nprocs() - NPROCS) + 1: end]) : addprocs(NPROCS - nprocs());

@everywhere include(".\\src\\SWAXS.jl");
@everywhere using .SWAXS





@everywhere dir = "C:\\Users\\Yen-Lin\\Desktop\\solvent";

# make solvent
@time solvent_batch(dir, joinpath(dir,"cluster2.pdb"), "frame"; solventprefix="solvent", d=10.0);

# compute the buffer-subttracted
filelist = readdir(dir);
solventlist = joinpath.(dir, filelist[startswith.(filelist, "solvent")]);
q = collect(0.0:0.01:1.5);

@info("Computing SWAXS ... ");
@info("Starting Time: $(now())");
@time intensity = PDBSWAXS(joinpath(dir, "cluster2.pdb"), solventlist, q; J=1500);
writedlm(joinpath(dir, "cluster2.dat"), [q intensity]);
