module FLOWFMM

import Base.:^
using LinearAlgebra
using StaticArrays
using ForwardDiff
using ReverseDiff
using ChainRulesCore
#using DelimitedFiles
using WriteVTK

const ONE_OVER_4PI = 1/4/pi

for file in ["containers", "derivatives", "element", "tree", "direct", "spherical", "fmm", "sortwrapper", "compatibility", "visualize", "estimate_cost", "grad_from_chainrules_extended", "rrules_definitions", "complex_track"]
    include(file*".jl")
end

end # module