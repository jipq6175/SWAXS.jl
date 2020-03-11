module SWAXS

using Distributed, DelimitedFiles, StatsBase, LinearAlgebra, Statistics

const AVec = AbstractVector;
const AMat = AbstractMatrix;
const AFloat = AbstractFloat;
const AS = AbstractString;

export mrc, mrc_reader, DenSWAXS
export PDBSWAXS
export Voxel, readvox, ShapeSWAXS


include("densityswaxs.jl");
include("pdbswaxs.jl");
include("voxel.jl");
include("parser.jl");


# for building executable
Base.@ccallable function julia_main()::Cint

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

    NPROCS = args["npr"];
    @info("--- SWAXS: Setting up parallel workers ...");
    nprocs() > NPROCS ? rmprocs(workers()[end - (nprocs() - NPROCS) + 1: end]) : addprocs(NPROCS - nprocs());

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
        return 1;
    end

    writedlm(output * ".dat", [q intensity]);
    @info("--- SWAXS: SWAXS program completed successfully: elapsed time = $t seconds with $(args["npr"]) cores.")
    @info("--- SWAXS: Removing parallel workers ...");
    rmprocs(workers());
    return 0;
end



end # module
