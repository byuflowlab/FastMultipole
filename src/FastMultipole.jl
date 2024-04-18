module FastMultipole

import Base.:^
using LinearAlgebra, StaticArrays, WriteVTK#, DelimitedFiles

const ONE_OVER_4PI = 1/4/pi

# multithreading parameters
const MIN_NPT_B2M = 100
const MIN_NPT_M2M = 100
const MIN_NPT_M2L = 100
const MIN_NPT_L2L = 100
const MIN_NPT_L2B = 100
const MIN_NPT_NF = 100

for file in ["containers", "complex", "derivatives", "element", "tree", "direct", "spherical", "fmm", "sortwrapper", "compatibility", "b2m", "probe", "visualize"]#, "estimate_cost"]
    include(file*".jl")
end

export fmm!, direct!, Tree, SortWrapper, ProbeSystem, reset!, DerivativesSwitch

end # module
