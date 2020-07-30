
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
