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


# for building executable



end # module
