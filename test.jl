

using DelimitedFiles, Distributed, Plots
cd("G:\\My Drive\\18. Github Repo\\SWAXS.jl");
NPROCS = Sys.CPU_THREADS;
nprocs() > NPROCS ? rmprocs(workers()[end - (nprocs() - NPROCS) + 1: end]) : addprocs(NPROCS - nprocs());

@everywhere include(".\\src\\SWAXS.jl");
@everywhere using .SWAXS


q = collect(0.0:0.01:1.5);

# rna12
@time st = PDBSWAXS("solute.pdb", q; J=1200);
@time sv = PDBSWAXS("solvent.pdb", q; J=1200);
@time rna = PDBSWAXS("solute.pdb", q; J=1200, waters=false);
@time rna_bs = PDBSWAXS("solute.pdb", "solvent.pdb", q; J=1200);


scatter(q, rna, yscale=:log10)
plot!(q, rna_bs, yscale=:log10)


plot(q, st, yscale=:log10)
plot!(q, sv, yscale=:log10)


# dna25
@time dna = PDBSWAXS("frame1.pdb", q; J=1200, waters=false);
@time dna_bs = PDBSWAXS("frame1.pdb", "frame1_solvent.pdb", q; J=1200);



scatter(q, dna, yscale=:log10);
plot!(q, dna_bs, yscale=:log10, lw=2.0)
plot!(q, dna_bs2, yscale=:log10, lw=3.0)
