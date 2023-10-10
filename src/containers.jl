#####
##### dispatch for common interface for external packages
#####
struct Body end
const BODY = Body()

struct Position end
const POSITION = Position()

struct Radius end
const RADIUS = Radius()

struct ScalarPotential end
const SCALAR_POTENTIAL = ScalarPotential()

struct VectorPotential end
const VECTOR_POTENTIAL = VectorPotential()

struct Velocity end
const VELOCITY = Velocity()

struct VelocityGradient end
const VELOCITY_GRADIENT = VelocityGradient()

struct Vertex1 end
const VERTEX1 = Vertex1()

struct Vertex2 end
const VERTEX2 = Vertex2()

struct Vertex3 end
const VERTEX3 = Vertex3()

struct Vertex4 end
const VERTEX4 = Vertex4()

struct ScalarStrength end
const SCALAR_STRENGTH = ScalarStrength()

struct VectorStrength end
const VECTOR_STRENGTH = VectorStrength()

##### 
##### dispatch convenience functions for multipole creation definition 
#####
abstract type AbstractKernel end

struct VortexPoint <: AbstractKernel end

struct VortexLine <: AbstractKernel end # not yet derived

struct VortexPanel <: AbstractKernel end # not yet derived

struct SourcePoint <: AbstractKernel end

struct SourcePanel <: AbstractKernel end # not yet implemented

struct DipolePanel <: AbstractKernel end # not yet implemented

#####
##### octree creation
#####
struct MultiBranch{TF,N}
    n_branches::Int8        # number of child branches
    n_bodies::SVector{N,Int32}         # number of descendent bodies
    first_branch::Int32     # index of the first branch
    first_body::SVector{N,Int32}       # index of the first element
    center::SVector{3,TF}   # center of the branch
    radius::TF              # side lengths of the cube encapsulating the branch
    multipole_expansion::NTuple{4,Vector{Complex{TF}}} # multipole expansion coefficients
    local_expansion::NTuple{4,Vector{Complex{TF}}}     # local expansion coefficients
    lock::ReentrantLock
    child_lock::ReentrantLock
end

struct SingleBranch{TF}
    n_branches::Int8        # number of child branches
    n_bodies::Int32         # number of descendent bodies
    first_branch::Int32     # index of the first branch
    first_body::Int32       # index of the first element
    center::SVector{3,TF}   # center of the branch
    radius::TF              # side lengths of the cube encapsulating the branch
    multipole_expansion::NTuple{4,Vector{Complex{TF}}} # multipole expansion coefficients
    local_expansion::NTuple{4,Vector{Complex{TF}}}     # local expansion coefficients
    lock::ReentrantLock
    child_lock::ReentrantLock
end

abstract type Tree end

"""
bodies[index_list] is the same sort operation as performed by the tree
sorted_bodies[inverse_index_list] undoes the sort operation performed by the tree
"""
struct MultiTree{TF,N} <: Tree
    branches::Vector{MultiBranch{TF,N}}        # a vector of `Branch` objects composing the tree
    expansion_order::Int16
    n_per_branch::Int32    # max number of bodies in a leaf
    index_list::NTuple{N,Vector{Int}}
    inverse_index_list::NTuple{N,Vector{Int}}
    leaf_index::Vector{Int}
    cumulative_count::Vector{Int} # starting with 0, a cumulative accounting of how many bodies are in leaf branches
end

struct SingleTree{TF} <: Tree
    branches::Vector{SingleBranch{TF}}        # a vector of `Branch` objects composing the tree
    expansion_order::Int16
    n_per_branch::Int32    # max number of bodies in a leaf
    index::Vector{Int}
    inverse_index::Vector{Int}
    leaf_index::Vector{Int}
    cumulative_count::Vector{Int} # starting with 0, a cumulative accounting of how many bodies are in leaf branches
end

#####
##### allow input systems to take any form when desired
#####
struct SortWrapper{TS}
    system::TS
    index::Vector{Int}
end

function SortWrapper(system)
    return SortWrapper(system,collect(1:length(system)))
end