module FLOWFMM

import Base.:^
using LinearAlgebra, StaticArrays, WriteVTK, DelimitedFiles

const ONE_OVER_4PI = 1/4/pi

for file in ["containers", "derivatives", "element", "tree", "direct", "spherical", "fmm", "sortwrapper", "compatibility", "visualize", "estimate_cost"]
    include(file*".jl")
end

end # module
