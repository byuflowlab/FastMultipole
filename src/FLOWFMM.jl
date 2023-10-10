module FLOWFMM

import Base.:^
using LinearAlgebra
using StaticArrays

const ONE_OVER_4PI = 1/4/pi

for file in ["containers", "derivatives", "element", "tree", "direct", "spherical", "fmm", "sortwrapper", "compatibility"]
    include(file*".jl")
end

end # module
