

using DelimitedFiles, Plots


# q, rna, rna25k400, rna25na100, rna25mg10 = load("rna25.jld", "q", "rna", "rna25k400", "rna25na100", "rna25mg10");
# q, dna, dna25k400, dna25mg05, dna25mg2 = load("dna25.jld", "q", "dna", "dna25k400", "dna25mg05", "dna25mg2");

dir = "C:\\Users\\Yen-Lin\\Box Sync\\01_WorkData\\NA25WAXS\\YenMixDNA_400mM_K";

datlist = readdir(dir; join=true);
datlist = datlist[endswith.(datlist, ".dat")];

q = collect(0.0:0.005:1.5);
n = length(datlist);
data = Array{Float64, 2}(undef, length(q), n);

for i = 1:n
    m = readdlm(datlist[i]);
    data[:, i] = m[:, 2];
end


plot(q, data[:,1:end], yscale=:log10)
