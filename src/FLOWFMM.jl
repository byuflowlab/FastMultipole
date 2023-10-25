module FLOWFMM

import Base.:^
using LinearAlgebra
using StaticArrays
using ForwardDiff
using ReverseDiff
using ChainRulesCore

const ONE_OVER_4PI = 1/4/pi

for file in ["containers", "derivatives", "element", "tree", "direct", "spherical", "fmm", "sortwrapper", "compatibility", "complex_track", "rrules_definitions"]
    include(file*".jl")
end

end # module