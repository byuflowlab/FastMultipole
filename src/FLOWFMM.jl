module FLOWFMM

import Base:iterate, length, ^, eltype
import Statistics
S = Statistics

for file in ["misc", "kernel", "element", "tree", "direct", "gravitational"]
    include(file*".jl")
end

#= Options:

* n_divisions
* partition_algorithm
*

=#

end # module
