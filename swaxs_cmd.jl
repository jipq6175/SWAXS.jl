#swaxs.jl

# setting up user's environment
include(".\\src\\parser.jl");

function main()

    args = parse_commandline();

    qmin = args["qmin"];
    qmax = args["qmax"];
    qspacing = args["qspacing"];
    q = collect(qmin:qspacing:qmax);
    J = args["J"];

    density = args["density"];
    pdb = args["pdb"];
    solute = args["solute"];
    solvent = args["solvent"];
    binvox = args["binvox"];
    solvent_density = args["solvent_density"];
    voxel_density = args["voxel_density"];
    density_cutoff = args["density_cutoff"];
    output = args["output"];

    @info("--- SWAXS: Please wait ... ");
    if !isnothing(density) && isnothing(pdb) && isnothing(solute) && isnothing(solvent) && isnothing(binvox)
        @info("--- SWAXS: Computing SWAXS (J=$J) using electron density file: $density, with sden=$solvent_density cutoff=$density_cutoff. ");
        m = mrc_reader(density);
        t = @elapsed intensity = DenSWAXS(m, q; density_cutoff=density_cutoff, J=J, sden=solvent_density);

    elseif isnothing(density) && !isnothing(pdb) && isnothing(solute) && isnothing(solvent) && isnothing(binvox)
        @info("--- SWAXS: Computing SWAXS (J=$J) using single PDB file: $pdb. ");
        t = @elapsed intensity = PDBSWAXS(pdb, q; J=J);

    elseif isnothing(density) && isnothing(pdb) && !isnothing(solute) && !isnothing(solvent) && isnothing(binvox)
        @info("--- SWAXS: Computing SWAXS (J=$J) using solute: $solute and solvent: $solvent. ");
        t = @elapsed intensity = PDBSWAXS(solute, solvent, q; J=J);

    elseif isnothing(density) && isnothing(pdb) && isnothing(solute) && isnothing(solvent) && !isnothing(binvox)
        @info("--- SWAXS: Computing SWAXS (J=$J) using shape file: $binvox, with voxel_density=$voxel_density, sden=$solvent_density.");
        v = readvox(binvox);
        t = @elapsed intensity = ShapeSWAXS(v, voxel_density, q; J=J, sd=solvent_density);

    else
        @warn("--- SWAXS: Observed conflicting options or didn't parse a pair of solute and solvent files. ");
        return nothing;
    end

    writedlm(output * ".dat", [q intensity]);
    @info("--- SWAXS: SWAXS program completed successfully: elapsed time = $t seconds with $(args["npr"]) cores.")
    @info("--- SWAXS: Removing parallel workers ...");
    rmprocs(workers());
    return nothing;
end



# toplevel code
args = parse_commandline();

# check the options
NPROCS = args["npr"];
@info("--- SWAXS: Setting up parallel workers ...");
nprocs() > NPROCS ? rmprocs(workers()[end - (nprocs() - NPROCS) + 1: end]) : addprocs(NPROCS - nprocs());
@everywhere include(".\\src\\SWAXS.jl");
@everywhere using .SWAXS

main();





# cd("G:\\My Drive\\18. Github Repo\\SWAXS.jl\\");
# #
# using Distributed, DelimitedFiles
# NPROCS = Sys.CPU_THREADS;
# @info("SWAXS: Setting up parallel workers ...");
# nprocs() > NPROCS ? rmprocs(workers()[end - (nprocs() - NPROCS) + 1: end]) : addprocs(NPROCS - nprocs());
# #
# @everywhere include("swaxs_module.jl");
# @everywhere using .SWAXS
#
# q = collect(0.0:0.02:1.0);
#
# # m = mrc_reader("..\\swaxs_testdata\\dna_denss1.mrc");
# # @time DenSWAXS(m, q);
#
#
#
# solutefn = "..\\swaxs_testdata\\H2helix30mMMg_state00001_solute.pdb";
# solventfn = "..\\swaxs_testdata\\H2helix30mMMg_state00001_solvent.pdb";
# @time intensity = PDBSWAXS(solutefn, solventfn, q);
# dir = "C:\\Users\\Yen-Lin\\Desktop\\New folder";
# mat = readdlm(joinpath(dir, "dna_a.dat"));
# q = mat[:, 1];
# v = readvox(joinpath(dir, "wine_bottle.binvox"));
# @time intensity = ShapeSWAXS(v, 0.5, q);
# writedlm(joinpath(dir, "wine_bottle.dat"), [q intensity 0.02*intensity]);
