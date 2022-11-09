module FLOWFMM

import Base.:^
# import Statistics as S
using LinearAlgebra
using StaticArrays

# const THETA = 4
const ONE_OVER_4PI = 1/4/pi

# indices
const i_POSITION = 1:3
const i_STRENGTH = 4:7
const i_POTENTIAL = 1:4
const i_POTENTIAL_JACOBIAN = 5:16
const i_POTENTIAL_HESSIAN = 17:52
const i_VELOCITY = 1:3

for file in ["misc", "derivatives", "potential", "element", "tree", "direct", "spherical", "fmm"]
    include(file*".jl")
end

end # module
