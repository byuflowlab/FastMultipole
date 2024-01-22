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

for file in ["containers", "derivatives", "element", "tree", "direct", "spherical", "fmm", "sortwrapper", "compatibility", "b2m", "visualize", "estimate_cost", "rrules_definitions"]
    include(file*".jl")
end

end # module