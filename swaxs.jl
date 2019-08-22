# sawxs module


module SWAXS

using Distributed, DelimitedFile, StatsBase, LinearAlgebra, Statistics, Reexport

const AVec = AbstractVector;
const AMat = AbstractMatrix;
const AFloat = AbstractFloat;









include("mrc.jl");
@reexport using .MRC
include("densityswaxs.jl");


include("voxel.jl");
include("dq.jl");







end  # module SWAXS
