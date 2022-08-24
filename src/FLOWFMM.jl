module FLOWFMM

import Base:^
import Statistics
S = Statistics

for file in ["misc", "kernel", "element", "tree", "direct", "fmm", "gravitational"]
    include(file*".jl")
end

#= Options:

* n_divisions
* partition_algorithm
*

=#

end # module
