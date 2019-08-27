# sawxs module


module SWAXS

using Distributed, DelimitedFiles, StatsBase, LinearAlgebra, Statistics, Reexport

const AVec = AbstractVector;
const AMat = AbstractMatrix;
const AFloat = AbstractFloat;
const AS = AbstractString;

export DenSWAXS
export PDBSWAXS

include("mrc.jl");
@reexport using .MRC
include("densityswaxs.jl");
include("pdbswaxs.jl");


end  # module SWAXS
