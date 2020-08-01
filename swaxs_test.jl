
# batch swaxs

using DelimitedFiles, Distributed
cd("G:\\My Drive\\18. Github Repo\\SWAXS.jl");
NPROCS = Sys.CPU_THREADS;
nprocs() > NPROCS ? rmprocs(workers()[end - (nprocs() - NPROCS) + 1: end]) : addprocs(NPROCS - nprocs());

@everywhere include(".\\src\\SWAXS.jl");
@everywhere using .SWAXS



## testing

dir = "C:\\Users\\yc225\\Desktop\\test";
@time sys = SimplyPDB(joinpath(dir, "frame1.pdb"));
@time A = SWAXS._SA(sys, 1.0);
buffers = readdir(dir);
buffers = joinpath.(dir, buffers[startswith.(buffers, "frame")]);
@time A = SWAXS._SA(buffers[1:5], 1.0);
z = zeros(2, 1500, 101);
for i = 1:101
    @info(" i = $i ... ");
    @time z[:, :, i] = SWAXS._SA(SimplyPDB(buffers[i]), 1.0);
end

v = zeros(2, 1500, 101);
for i = 1:101
    @info(" i = $i ... ");
    @time v[:, :, i] = SWAXS._SA(SimplyPDB(buffers[i]), 0.3);
end

w = zeros(2, 1500, 101);
for i = 1:101
    @info(" i = $i ... ");
    @time w[:, :, i] = SWAXS._SA(SimplyPDB(buffers[i]), 0.7);
end

x = zeros(2, 1500, 101);
for i = 1:101
    @info(" i = $i ... ");
    @time x[:, :, i] = SWAXS._SA(SimplyPDB(buffers[i]), 1.25);
end

y = zeros(2, 1500, 101);
for i = 1:101
    @info(" i = $i ... ");
    @time y[:, :, i] = SWAXS._SA(SimplyPDB(buffers[i]), 1.5);
end


amp = cumsum(z, dims=3);
for i = 1:101
    amp[:, :, i] = amp[:, :, i] ./ i;
end
fintz = sum(amp[1,:,:].^2 + amp[2, :, :].^2, dims=1)[1, :];

amp = cumsum(v, dims=3);
for i = 1:101
    amp[:, :, i] = amp[:, :, i] ./ i;
end
fintv = sum(amp[1,:,:].^2 + amp[2, :, :].^2, dims=1)[1, :];

amp = cumsum(w, dims=3);
for i = 1:101
    amp[:, :, i] = amp[:, :, i] ./ i;
end
fintw = sum(amp[1,:,:].^2 + amp[2, :, :].^2, dims=1)[1, :];

amp = cumsum(x, dims=3);
for i = 1:101
    amp[:, :, i] = amp[:, :, i] ./ i;
end
fintx = sum(amp[1,:,:].^2 + amp[2, :, :].^2, dims=1)[1, :];

amp = cumsum(z, dims=3);
for i = 1:101
    amp[:, :, i] = amp[:, :, i] ./ i;
end
finty = sum(amp[1,:,:].^2 + amp[2, :, :].^2, dims=1)[1, :];


## average intensity at q = 1.0 by the solvents
TITLE = font(20, "Times");
TICKS = font(14, "Times");
GUIDE = font(16, "Times");
LEGEN = font(14, "Times");

#plot!(collect(1:101), fintv./1e8, lw=4.0, lab="q = 0.3");
p = plot(collect(1:101), fintw./1e8, lw=4.0, lab="q = 0.7");
plot!(collect(1:101), fintz./1e8, lw=4.0, lab="q = 1.0");
plot!(collect(1:101), fintx./1e8, lw=4.0, lab="q = 1.25");
plot!(collect(1:101), finty./1e8, lw=4.0, lab="q = 1.5");
ylims!(0, 3);
ylabel!("I(q=1.0) (A.U.)");
xlabel!("# of Solvent Frames");
title!("Estimated Intensity from Bulk");
plot!(size=(600, 600), dpi=600, grid=false, legend=:topright, xtickfont=TICKS, ytickfont=TICKS, titlefont=TITLE, guidefont=GUIDE, legendfont=LEGEN, framestyle=:box)
savefig(p, "solvent_frames.png");
savefig(p, "solvent_frames.pdf");


## subtracts
sintv = zeros(101);
sintw = zeros(101);
sintx = zeros(101);
sinty = zeros(101);
sintz = zeros(101);
# for i = 2:101
#     subtracted = z[:,:,1:i] .- mean(z[:,:,1:i], dims=3);
#     sint[i] = sum(sum(subtracted[1, :, :].^2 + subtracted[2, :, :].^2, dims=1)[1, :]) * (i + 1) / i / (i - 1);
# end

for i = 2:101
    subtracted = v[:,:,1:i] .- mean(v[:,:,1:i], dims=3);
    tmp = mean(subtracted, dims=2)[:, 1, :];
    sintv[i] = sum(tmp[1, :].^2 + tmp[2, :].^2) * (i + 1) / i / (i - 1);
end
for i = 2:101
    subtracted = w[:,:,1:i] .- mean(w[:,:,1:i], dims=3);
    tmp = mean(subtracted, dims=2)[:, 1, :];
    sintw[i] = sum(tmp[1, :].^2 + tmp[2, :].^2) * (i + 1) / i / (i - 1);
end
for i = 2:101
    subtracted = x[:,:,1:i] .- mean(x[:,:,1:i], dims=3);
    tmp = mean(subtracted, dims=2)[:, 1, :];
    sintx[i] = sum(tmp[1, :].^2 + tmp[2, :].^2) * (i + 1) / i / (i - 1);
end
for i = 2:101
    subtracted = y[:,:,1:i] .- mean(y[:,:,1:i], dims=3);
    tmp = mean(subtracted, dims=2)[:, 1, :];
    sinty[i] = sum(tmp[1, :].^2 + tmp[2, :].^2) * (i + 1) / i / (i - 1);
end
for i = 2:101
    subtracted = z[:,:,1:i] .- mean(z[:,:,1:i], dims=3);
    tmp = mean(subtracted, dims=2)[:, 1, :];
    sintz[i] = sum(tmp[1, :].^2 + tmp[2, :].^2) * (i + 1) / i / (i - 1);
end

p = plot(collect(2:101), sintv[2:101]./1e8, lw=4.0, lab="q = 0.3");
plot!(collect(2:101), sintw[2:101]./1e8, lw=4.0, lab="q = 0.7");
plot!(collect(2:101), sintz[2:101]./1e8, lw=4.0, lab="q = 1.0");
plot!(collect(2:101), sintx[2:101]./1e8, lw=4.0, lab="q = 1.25");
plot!(collect(2:101), sinty[2:101]./1e8, lw=4.0, lab="q = 1.5");
ylims!(0, 5e-6);
ylabel!("I(q-1.0) (A.U.)");
xlabel!("# of Solvent Frames");
title!("Estimated Intensity from Solvent Fluctuation");
plot!(size=(600, 600), dpi=600, grid=false, legend=:topright, xtickfont=TICKS, ytickfont=TICKS, titlefont=TITLE, guidefont=GUIDE, legendfont=LEGEN, framestyle=:box)
savefig(p, "solvent_fluctuation.png");
savefig(p, "solvent_fluctuation.pdf");


@time intensity = PDBSWAXS(joinpath(dir, "solute1.pdb"), buffers[1:5], collect(0.2:0.2:1.0); J = 1000);





# function solvent_batch(dir::AbstractString, soluteprefix::AbstractString, sourcefn::AbstractString; solventprefix::AbstractString="solvent")
#
#     filelist = readdir(dir);
#     filelist = filelist[endswith.(filelist, ".pdb")];
#
#     solutelist = filelist[startswith.(filelist, soluteprefix)];
#     nsolute = length(solutelist);
#     source = joinpath(dir, sourcefn);
#
#     for i = 1:nsolute
#         s = solutelist[i];
#         @info("Making solvent for $s ... ");
#         outputfn = joinpath(dir, replace(s, soluteprefix => solventprefix));
#         isfile(outputfn) ? @warn("$outputfn exists, skipped ... ") : make_solvent(joinpath(dir, s), source, 2.0, outputfn);
#     end
#
#     return nothing;
# end
#
#
#
#
#
#
# function swaxs_batch(dir::AbstractString, q::AbstractVector; J::Signed=1500, soluteprefix::AbstractString="frame", solventprefix::AbstractString="bulk")
#
#     filelist = readdir(dir);
#     filelist = filelist[endswith.(filelist, ".pdb")];
#
#     solutelist = filelist[startswith.(filelist, soluteprefix)];
#     solventlist = filelist[startswith.(filelist, solventprefix)];
#
#     nsolute = length(solutelist);
#     nsolvent = length(solventlist);
#
#     if nsolvent == 1
#         @info("Using one single solvent frame for all the solutes ...");
#         solventlist = fill(solventlist[1], nsolute);
#     elseif nsolvent == 0
#         @warn("No solvent found ... ");
#         return nothing;
#     elseif nsolvent != nsolute
#         @warn("Solvent and solute are not matching ... ");
#         return nothing;
#     else
#         @info("Solute and solvent match ...");
#     end
#
#     @info("Starting to compute the SWAXS profiles of $nsolute frames .....");
#     for i = 1:nsolute
#
#
#         datfile = joinpath(dir, solutelist[i][1:end-4] * ".dat");
#         if isfile(datfile)
#             @warn("The SWAXS profile for $(solutelist[i]) exists, skipping this frame ...");
#         else
#             t = @elapsed mat = PDBSWAXS(joinpath(dir, solutelist[i]), joinpath(dir, solventlist[i]), q; J=J);
#             #t = @elapsed mat = PDBSWAXS(joinpath(dir, solutelist[i]), q; J=J, waters=false, ions=false);
#             writedlm(datfile, [q mat]);
#             @info("The SWAXS profile of {$(solutelist[i]), $(solventlist[i])} is done with $(round(t, digits=2)) seconds ...");
#         end
#
#     end
#
#     return nothing;
# end
#
#
# dir = "C:\\Users\\Yen-Lin\\Box Sync\\01_WorkData\\TPXWAXS\\TripnUAU_200mMNaCl";
#
# # q = collect(0.0:0.02:1.5);
# # solute = SimplyPDB(joinpath(dir, "frame1.pdb"));
# # solvent = SimplyPDB(joinpath(dir, "solvent1.pdb"));
# # @time mg = PDBSWAXS(joinpath(dir, "frame1.pdb"), joinpath(dir, "solvent1.pdb"), q; J=1000);
#
#
# # Batch solvent
# @time solvent_batch(dir, "frame", "bulk.pdb");
#
#
# q = collect(0.0:0.005:1.5);
# swaxs_batch(dir, q; J=1750, soluteprefix="frame", solventprefix="solvent");
