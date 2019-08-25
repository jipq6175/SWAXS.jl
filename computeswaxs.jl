#swaxs.jl

# setting up user's environment

cd("G:\\My Drive\\18. Github Repo\\SWAXS.jl\\");


using Distributed, DelimitedFiles
NPROCS = Sys.CPU_THREADS;
@info("SWAXS: Setting up parallel workers ...");
nprocs() > NPROCS ? rmprocs(workers()[end - (nprocs() - NPROCS) + 1: end]) : addprocs(NPROCS - nprocs());

@everywhere include("swaxs.jl");
@everywhere using .SWAXS

q = collect(0.0:0.002:1.4);

#m = mrc_reader("..\\swaxs_testdata\\dna_denss1.mrc");
#@time DenSWAXS(m, q);

v = readvox("..\\swaxs_testdata\\cat.binvox");
@time intensity = ShapeSWAXS(v, 1.0, q);
writedlm(".\\swaxs_testdata\\cat.dat", [q intensity 0.02 * intensity]);
