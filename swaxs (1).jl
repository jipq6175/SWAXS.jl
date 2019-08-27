#swaxs.jl

# setting up user's environment
using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--density", "-d"
            help = "Electron density map file: .mrc or .map"
            arg_type = String
            default = ""
        "--pdb", "-p"
            help = "Single file containing atomic coordinates: .pdb"
            arg_type = String
            default = ""
        "--solute", "-t"
            help = "The solute.pdb file"
            arg_type = String
            default = ""
        "--solvent", "-v"
            help = "The solvent.pdb file"
            arg_type = String
            default = ""
        "--solvent_density", "-s"
            help = "The bulk solvent electron density in e/A^3, default = 0.335"
            arg_type = Float64
            default = 0.335
        "--density_cutoff", "-c"
            help = "The electron density cutoff to exclude near-zero voxels, default = 0.001"
            arg_type = Float64
            default = 0.001
        "--J", "-J"
            help = "Number of orientations to be averaged, default = 1200"
            arg_type = Int64
            default = 1200
    end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline();
    println("Parsed args:");
    for (arg,val) in parsed_args
        println("  $arg  =>  $val");
    end
    return nothing;
end

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
