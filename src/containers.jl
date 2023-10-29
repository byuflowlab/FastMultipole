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
abstract type Branch end

struct MultiBranch{TF,N} <: Branch
    bodies_index::SVector{N,UnitRange{Int64}}
    n_branches::Int64
    branch_index::UnitRange{Int64}
    center::SVector{3,TF}   # center of the branch
    radius::TF              # side lengths of the cube encapsulating the branch
    multipole_expansion::Array{Complex{TF},2} # multipole expansion coefficients
    local_expansion::Array{Complex{TF},2}     # local expansion coefficients
    lock::ReentrantLock
end

Base.eltype(::MultiBranch{TF}) where TF = TF

struct SingleBranch{TF} <: Branch
    bodies_index::UnitRange{Int64}
    n_branches::Int64
    branch_index::UnitRange{Int64}
    center::SVector{3,TF}   # center of the branch
    radius::TF              # side lengths of the cube encapsulating the branch
    multipole_expansion::Array{Complex{TF},2} # multipole expansion coefficients
    local_expansion::Array{Complex{TF},2}     # local expansion coefficients
    lock::ReentrantLock
end

Base.eltype(::SingleBranch{TF}) where TF = TF

abstract type Tree end

"""
bodies[index_list] is the same sort operation as performed by the tree
sorted_bodies[inverse_index_list] undoes the sort operation performed by the tree
"""
struct MultiTree{TF,N,TB} <: Tree
    branches::Vector{MultiBranch{TF,N}}        # a vector of `Branch` objects composing the tree
    levels_index::Vector{UnitRange{Int64}}
    sort_index_list::NTuple{N,Vector{Int}}
    inverse_sort_index_list::NTuple{N,Vector{Int}}
    buffers::TB
    expansion_order::Int64
    n_per_branch::Int64    # max number of bodies in a leaf
end

struct SingleTree{TF,TB} <: Tree
    branches::Vector{SingleBranch{TF}}        # a vector of `Branch` objects composing the tree
    levels_index::Vector{UnitRange{Int64}}
    sort_index::Vector{Int64}
    inverse_sort_index::Vector{Int64}
    buffer::TB
    expansion_order::Int64
    n_per_branch::Int64    # max number of bodies in a leaf
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