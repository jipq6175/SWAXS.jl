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
    bulkdir = args["bulkdir"];
    envelope = args["envelope"];
    prefix = args["prefix"];
    dare = args["dare"];
    binvox = args["binvox"];
    solvent_density = args["solvent_density"];
    voxel_density = args["voxel_density"];
    density_cutoff = args["density_cutoff"];
    output = args["output"];

    NPROCS = args["npr"];
    # @info("--- SWAXS: Setting up parallel workers ...");
    # nprocs() > NPROCS ? rmprocs(workers()[end - (nprocs() - NPROCS) + 1: end]) : addprocs(NPROCS - nprocs());
    #
    # @info("--- SWAXS: Please wait ... ");

    # Density SWAXS using .mrc
    # Command:
    # (1) julia -O 3 --color=yes swaxs_cmd.jl --density ".\\test\\malat1.mrc" -s 0.335 -c 0.001 -o test 0.0 0.1 1.0
    # (2) julia -O 3 --color=yes swaxs_cmd.jl --density ".\\test\\dna1.mrc" -s 0.335 -c 0.001 -o test 0.0 0.1 1.0
    if !isnothing(density) && isnothing(pdb) && isnothing(solute) && isnothing(solvent) && isnothing(binvox) && !dare
        @info("--- SWAXS: Computing SWAXS (J=$J) using electron density file: $density, with sden=$solvent_density cutoff=$density_cutoff. ");
        m = mrc_reader(density);
        t = @elapsed intensity = DenSWAXS(m, q; density_cutoff=density_cutoff, J=J, sden=solvent_density);

    # PDBSWAXS using 1 PDB
    # Command:
    # (1) julia -O 3 --color=yes swaxs_cmd.jl --pdb ".\\test\\rna.pdb" -o test 0.0 0.1 1.0
    elseif isnothing(density) && !isnothing(pdb) && isnothing(solute) && isnothing(solvent) && isnothing(binvox) && !dare
        @info("--- SWAXS: Computing SWAXS (J=$J) using single PDB file: $pdb. ");
        t = @elapsed intensity = PDBSWAXS(pdb, q; J=J);

    # PDBSWAXS using solute and solvent
    # Command:
    # (1) julia -O 3 --color=yes swaxs_cmd.jl --solute ".\\test\\solute.pdb" --solvent ".\\test\\solvent.pdb" -o test 0.0 0.1 1.0
    elseif isnothing(density) && isnothing(pdb) && !isnothing(solute) && !isnothing(solvent) && isnothing(binvox) && !dare
        @info("--- SWAXS: Computing SWAXS (J=$J) using solute: $solute and solvent: $solvent. ");
        @info("--- SWAXS: Processing solvent ... ");
        outputfn = "tmp.pdb";
        isfile(outputfn) ? @warn("--- SWAXS: $outputfn exists, parsing ... ") : make_solvent(solute, solvent, 2.0, outputfn);
        t = @elapsed intensity = PDBSWAXS(solute, outputfn, q; J=J);
        rm(outputfn);

    # dare mode
    # Command:
    # (1) julia -O 3 --color=yes swaxs_cmd.jl --dare --bulkdir ".\\test\\solvent" --prefix frame --solute ".\\test\\solvent\\solute.pdb" --envelope 10.5 -J 800 -o test 0.0 0.1 1.0
    elseif isnothing(density) && isnothing(pdb) && !isnothing(solute) && isnothing(solvent) && isnothing(binvox) && dare
        if !isnothing(bulkdir) && !isnothing(prefix)
            @info("--- SWAXS: DARE mode ... ");
            @info("--- SWAXS: Processing bulk solvents ...");
            # solventlist = joinpath.(bulkdir, readdir(bulkdir));
            solvent_batch(bulkdir, solute, prefix; solventprefix="tmp", d=envelope);
            solventlist = readdir(bulkdir);
            solventlist = joinpath.(bulkdir, solventlist[startswith.(solventlist, "tmp")]);

            @info("--- SWAXS: Computing SWAXS (J=$J) in DARE mode ... ");
            t = @elapsed intensity = PDBSWAXS(solute, solventlist, q; J=J);
            rm.(solventlist);

        else
            @warn("--- SWAXS: DARE mode found missing bulk directory or bulk prefix... ");
            return 1;
        end


    # ShapeSWAXS using .binvox
    # Command
    # (1) julia -O 3 --color=yes swaxs_cmd.jl --binvox ".\\test\\rabbit.binvox" -s 0.335 -v 1.0 -o test 0.0 0.1 1.0
    elseif isnothing(density) && isnothing(pdb) && isnothing(solute) && isnothing(solvent) && !isnothing(binvox) && !dare
        @info("--- SWAXS: Computing SWAXS (J=$J) using shape file: $binvox, with voxel_density=$voxel_density, sden=$solvent_density.");
        v = readvox(binvox);
        t = @elapsed intensity = ShapeSWAXS(v, voxel_density, q; J=J, sd=solvent_density);

    else
        @warn("--- SWAXS: Observed conflicting options or didn't parse a pair of solute and solvent files. ");
        return 1;
    end

    writedlm(output * ".dat", [q intensity]);
    @info("--- SWAXS: SWAXS program completed successfully: elapsed time = $(round(t, digits=2)) seconds with $(args["npr"]) cores.")
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
