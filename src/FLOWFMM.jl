module FLOWFMM

import Base.:^
# import Statistics as S
# import LinearAlgebra as LA

const THETA = 4

for file in ["misc", "kernel", "element", "tree", "direct", "spherical", "fmm"]
    include(file*".jl")
end

end # module
