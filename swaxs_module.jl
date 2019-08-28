# sawxs module

__precompile__();

module SWAXS

using Distributed, DelimitedFiles, StatsBase, LinearAlgebra, Statistics, Reexport

const AVec = AbstractVector;
const AMat = AbstractMatrix;
const AFloat = AbstractFloat;
const AS = AbstractString;

export DenSWAXS
export PDBSWAXS
export Voxel, readvox, ShapeSWAXS

include("mrc.jl");
@reexport using .MRC
include("densityswaxs.jl");
include("pdbswaxs.jl");
include("voxel.jl");



end  # module SWAXS
