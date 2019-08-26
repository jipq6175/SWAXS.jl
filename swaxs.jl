# sawxs module


module SWAXS

using Distributed, DelimitedFiles, StatsBase, LinearAlgebra, Statistics, Reexport

const AVec = AbstractVector;
const AMat = AbstractMatrix;
const AFloat = AbstractFloat;
const AS = AbstractString;



export DenSWAXS
export PDBSWAXS, SimplyPDB, AtomAFF




include("mrc.jl");
@reexport using .MRC
include("densityswaxs.jl");

include("pdbswaxs.jl");
#include("scattering_form_factor.jl");






end  # module SWAXS
