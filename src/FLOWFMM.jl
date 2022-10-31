module FLOWFMM

import Base.:^
# import Statistics as S
# import LinearAlgebra as LA

# const THETA = 4
const ONE_OVER_4PI = 1/4/pi

for file in ["misc", "potential", "element", "tree", "direct", "spherical", "fmm"]
    include(file*".jl")
end

end # module
