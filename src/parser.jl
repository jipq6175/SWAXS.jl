# parser.jl


using ArgParse

function parse_commandline()
    s = ArgParseSettings()
    s.prog = "SWAXS";
    s.description = "Small and Wide-Angle X-ray Scattering calculator for (a) single atomic coordinates: .pdb, (b) solute-solvent PDB pair for buffer subtraction, (c) electron density: .mrc or .map with implicit density model, (d) voxelized 3D shape: .binvox with uniform excessive electron density and (f) DARE mode for accurate bulk solvent modeling. SWAXS implements the Debye formula, orientational average and theory of excessive electron density to account for buffer subtraction, solvent shell and excluded volumes. ";
    s.preformatted_description = false;
    s.epilog = "===================================================================\n===             Last Update: 08/10/20, Ithaca, NY.              ===\n===           Copyright (c) Yen-Lin Chen, 2018 - 2020           ===\n===                Academic Free License v. 3.0                 ===\n===                 Email: yc2253@cornell.edu                   ===\n===================================================================";
    s.preformatted_epilog = true;
    s.version = "0.3.1";
    s.add_version = true;
    s.add_help = true;

    @add_arg_table! s begin

        # densitySWAXS
        "--density", "-D"
            help = "CCP4 electron density map file: .mrc or .map"
            arg_type = String
        "--solvent_density", "-s"
            help = "The bulk solvent electron density in e/A^3"
            arg_type = Float64
            default = 0.335
        "--density_cutoff", "-c"
            help = "The electron density cutoff to exclude near-zero voxels"
            arg_type = Float64
            default = 0.001

        # single PDB
        "--pdb", "-P"
            help = "Single file containing atomic coordinates: .pdb"
            arg_type = String

        # solute and solvent pair
        "--solute", "-T"
            help = "The solute.pdb file containing atomic coordinates of molecules, ions and solvents"
            arg_type = String
        "--solvent", "-V"
            help = "The solvent.pdb file containing atomic coordinates of randomized bulk solvents"
            arg_type = String

        # md estimation of bulk solvent
        "--bulkdir", "-b"
            help = "The directory containing bulk frames; more than 50 frames are suggested or use --solute and --solvent"
            arg_type = String
        "--envelope", "-e"
            help = "Distance between the envelope and molecular surface, i.e. solvent layer width"
            arg_type = Float64
            default = 10.0;
        "--prefix", "-p"
            help = "The .pdb file prefix in the --bulkdir"
            arg_type = String
            default = "bulk"
        "--dare"
            help = "Make sure that you're actually doing it!! Computing cluster suggested."
            action = :store_true

        # ShapeSWAXS
        "--binvox", "-B"
            help = "The shape file of dummy voxels: .binvox"
            arg_type = String
        "--voxel_density", "-v"
            help = "The averaged electron density on dummy voxels in e/A^3"
            arg_type = Float64
            default = 0.5;

        # Common
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
