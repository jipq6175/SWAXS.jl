#swaxs.jl

# setting up user's environment
using ArgParse, Distributed, DelimitedFiles

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--density", "-d"
            help = "Electron density map file: .mrc or .map"
            arg_type = String
        "--pdb", "-p"
            help = "Single file containing atomic coordinates: .pdb"
            arg_type = String
        "--solute", "-t"
            help = "The solute.pdb file"
            arg_type = String
        "--solvent", "-v"
            help = "The solvent.pdb file"
            arg_type = String
        "--solvent_density", "-s"
            help = "The bulk solvent electron density in e/A^3"
            arg_type = Float64
            default = 0.335
        "--density_cutoff", "-c"
            help = "The electron density cutoff to exclude near-zero voxels"
            arg_type = Float64
            default = 0.001
        "--J", "-J"
            help = "Number of orientations to be averaged"
            arg_type = Int64
            default = 1200
        "--npr", "-n"
            help = "Number of parallel workers for computation"
            arg_type = Int64
            default = Sys.CPU_THREADS
        "qmin"
            help = "Starting q value, in (1/A)"
            arg_type = Float64
            required = true
        "qspacing"
            help = "q grid spacing, in (1/A)"
            arg_type = Float64
            required = true
        "qmax"
            help = "Ending q value, in (1/A)"
            arg_type = Float64
            required = true
        "--output", "-o"
            help = "Output file prefix for saving the .dat file"
            arg_type = String
            default = "output"

    end

    return parse_args(s)
end

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
    solvent_density = args["solvent_density"];
    density_cutoff = args["density_cutoff"];
    output = args["output"];

    @info("Please wait ... ");
    if !isnothing(density) && isnothing(pdb) && isnothing(solute) && isnothing(solvent)
        @info("Computing SWAXS (J=$J) using electron density file: $density, with sden=$solvent_density cutoff=$density_cutoff. ");
        m = mrc_reader(density);
        t = @elapsed intensity = DenSWAXS(m, q; density_cutoff=density_cutoff, J=J, sden=solvent_density);

    elseif isnothing(density) && !isnothing(pdb) && isnothing(solute) && isnothing(solvent)
        @info("Computing SWAXS (J=$J) using single PDB file: $pdb. ");
        t = @elapsed intensity = PDBSWAXS(pdb, q; J=J);

    elseif isnothing(density) && isnothing(pdb) && !isnothing(solute) && !isnothing(solvent)
        @info("Computing SWAXS (J=$J) using solute: $solute and solvent: $solvent. ");
        t = @elapsed intensity = PDBSWAXS(solute, solvent, q; J=J);

    else
        @warn("Observed conflicting options or didn't parse a pair of solute and solvent files. ");
        return nothing;
    end

    writedlm(output * ".dat", [q intensity]);
    @info("SWAXS program completed successfully: elapsed time = $t seconds with $(args["npr"]) cores.")
    return nothing;
end



# toplevel code
args = parse_commandline();

# check the options
NPROCS = args["npr"];
@info("SWAXS: Setting up parallel workers ...");
nprocs() > NPROCS ? rmprocs(workers()[end - (nprocs() - NPROCS) + 1: end]) : addprocs(NPROCS - nprocs());
@everywhere include("swaxs_module.jl");
@everywhere using .SWAXS

main();





# cd("G:\\My Drive\\18. Github Repo\\SWAXS.jl\\");
#
# using Distributed, DelimitedFiles
# NPROCS = Sys.CPU_THREADS;
# @info("SWAXS: Setting up parallel workers ...");
# nprocs() > NPROCS ? rmprocs(workers()[end - (nprocs() - NPROCS) + 1: end]) : addprocs(NPROCS - nprocs());
#
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
