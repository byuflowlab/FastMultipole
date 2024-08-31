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

# preallocate y-axis rotation matrices by π/2
const Hs_π2 = Float64[1.0]

for file in ["containers", "complex", "derivatives", "element", "tree", "direct", "rotate", "spherical", "fmm", "sortwrapper", "compatibility", "b2m", "probe", "visualize"]#, "estimate_cost"]
    include(file*".jl")
end

# preallocate y-axis rotation matrices up to 20th order
update_Hs_π2!(Hs_π2, 20)

export fmm!, direct!, Tree, SortWrapper, ProbeSystem, add_line!, reset!, DerivativesSwitch

end # module
