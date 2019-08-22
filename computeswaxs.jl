#swaxs.jl

# setting up user's environment


using Distributed
NPROCS = Sys.CPU_THREADS;
@info("SWAXS: Setting up parallel workers ...");
nprocs() > NPROCS ? rmprocs(workers()[end - (nprocs() - NPROCS) + 1: end]) : addprocs(NPROCS - nprocs());

@everywhere include("swaxs.jl");
@everywhere using .SWAXS
