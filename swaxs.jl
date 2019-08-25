# sawxs module


module SWAXS

using Distributed, DelimitedFiles, StatsBase, LinearAlgebra, Statistics, Reexport

const AVec = AbstractVector;
const AMat = AbstractMatrix;
const AFloat = AbstractFloat;



export DenSWAXS
export ShapeSWAXS, readvox, Voxel
export PDBSWAXS, SimplyPDB




include("mrc.jl");
@reexport using .MRC
include("densityswaxs.jl");


include("voxel.jl");
include("shapeswaxs.jl");


#include("scattering_form_factor.jl");






end  # module SWAXS
