module FLOWFMM

import Base:^
import Statistics as S

const THETA = 4

for file in ["misc", "kernel", "element", "basis", "tree", "direct", "fmm"]
    include(file*".jl")
end

#= Options:

* n_divisions
* partition_algorithm
*

=#

end # module
