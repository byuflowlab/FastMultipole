module FastMultipole

import Base.:^
using LinearAlgebra
using StaticArrays
using ForwardDiff
using ReverseDiff
using ChainRulesCore
#using DelimitedFiles
using WriteVTK

const ONE_OVER_4PI = 1/4/pi

# multithreading parameters
const MIN_NPT_B2M = 100
const MIN_NPT_M2M = 100
const MIN_NPT_M2L = 100
const MIN_NPT_L2L = 100
const MIN_NPT_L2B = 100
const MIN_NPT_NF = 100

for file in ["containers", "complex", "derivatives", "element", "tree", "direct", "spherical", "fmm", "sortwrapper", "compatibility", "b2m", "probe", "visualize","grad_from_chainrules_extended","rrules_definitions"]#, "estimate_cost"]
    include(file*".jl")
end

include("/home/eric/Research/VPM_derivatives/src/safe_accumulate.jl") # add for real after testing.

export fmm!, direct!, Tree, SortWrapper, ProbeSystem, add_line!, reset!, DerivativesSwitch

end # module
