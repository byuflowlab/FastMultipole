#####
##### dispatch for common interface for external packages
#####
abstract type Indexable end

struct Body <: Indexable end
const BODY = Body()

struct Position <: Indexable end
const POSITION = Position()

struct Radius <: Indexable end
const RADIUS = Radius()

struct ScalarPotential <: Indexable end
const SCALAR_POTENTIAL = ScalarPotential()

struct VectorPotential <: Indexable end
const VECTOR_POTENTIAL = VectorPotential()

struct Velocity <: Indexable end
const VELOCITY = Velocity()

struct VelocityGradient <: Indexable end
const VELOCITY_GRADIENT = VelocityGradient()

struct Vertex <: Indexable end
const VERTEX = Vertex()

struct ScalarStrength <: Indexable end
const SCALAR_STRENGTH = ScalarStrength()

struct VectorStrength <: Indexable end
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
##### cost parameters
#####
# abstract type CostParameters end

# struct SingleCostParameters <: CostParameters
#     alloc_M2M_L2L::SVector{3,Float64}
#     tau_M2M_L2L::SVector{5,Float64}
#     alloc_M2L::SVector{3,Float64}
#     tau_M2L::SVector{5,Float64}
#     alloc_L2B::SVector{3,Float64}
#     tau_L2B::SVector{3,Float64}
#     C_nearfield::Float64
#     tau_B2M::SVector{3,Float64}
# end

# SingleCostParameters(;
#     alloc_M2M_L2L = ALLOC_M2M_L2L_DEFAULT, 
#     tau_M2M_L2L = TAU_M2M_DEFAULT, 
#     alloc_M2L = ALLOC_M2L_DEFAULT, 
#     tau_M2L = TAU_M2L_DEFAULT, 
#     tau_L2L = TAU_L2L_DEFAULT, 
#     alloc_L2B = ALLOC_L2B_DEFAULT, 
#     tau_L2B = TAU_L2B_DEFAULT, 
#     C_nearfield = C_NEARFIELD_DEFAULT,
#     tau_B2M = TAU_B2M_DEFAULT
# ) = SingleCostParameters(alloc_M2M_L2L, tau_M2M_L2L, alloc_M2L, tau_M2L, alloc_L2B, tau_L2B, C_nearfield, tau_B2M)

# struct MultiCostParameters{N} <: CostParameters
#     alloc_M2M_L2L::SVector{3,Float64}
#     tau_M2M_L2L::SVector{5,Float64}
#     alloc_M2L::SVector{3,Float64}
#     tau_M2L::SVector{5,Float64}
#     alloc_L2B::SVector{3,Float64}
#     tau_L2B::SVector{3,Float64}
#     C_nearfield::SVector{N,Float64}
#     tau_B2M::SVector{N,SVector{3,Float64}}
# end

# MultiCostParameters{N}(;
#     alloc_M2M_L2L = ALLOC_M2M_L2L_DEFAULT,
#     tau_M2M_L2L = TAU_M2M_L2L_DEFAULT,
#     alloc_M2L = ALLOC_M2L_DEFAULT,
#     tau_M2L = TAU_M2L_DEFAULT,
#     alloc_L2B = ALLOC_L2B_DEFAULT,
#     tau_L2B = TAU_L2B_DEFAULT,
#     C_nearfield = SVector{N,Float64}(C_NEARFIELD_DEFAULT for _ in 1:N),
#     tau_B2M = SVector{N,SVector{3,Float64}}(TAU_B2M_DEFAULT for _ in 1:N)
# ) where N = MultiCostParameters{N}(alloc_M2M_L2L, tau_M2M_L2L, alloc_M2L, tau_M2L, tau_L2L, alloc_L2B, tau_L2B, C_nearfield, tau_B2M)

# CostParameters(systems::Tuple) = MultiCostParameters()
# CostParameters(system) = SingleCostParameters()

# CostParameters(alloc_M2M_L2L, tau_B2M, alloc_M2L, tau_M2L, tau_L2L, alloc_L2B, tau_L2B, C_nearfield::Float64, tau_M2M_L2L) = 
#     SingleCostParameters(alloc_M2M_L2L, tau_B2M, alloc_M2L, tau_M2L, tau_L2L, alloc_L2B, tau_L2B, C_nearfield, tau_M2M_L2L)

# CostParameters(alloc_M2M_L2L, tau_B2M, alloc_M2L, tau_M2L, tau_L2L, alloc_L2B, tau_L2B, C_nearfield::SVector, tau_M2M_L2L) = 
#     MultiCostParameters(alloc_M2M_L2L, tau_B2M, alloc_M2L, tau_M2L, tau_L2L, alloc_L2B, tau_L2B, C_nearfield, tau_M2M_L2L)

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
    multipole_expansion::Array{TF,3} # multipole expansion coefficients
    local_expansion::Array{TF,3}     # local expansion coefficients
    harmonics::Array{TF,2}
    ML::MArray{Tuple{2,4},TF}
    lock::ReentrantLock
end

Base.eltype(::MultiBranch{TF}) where TF = TF

struct SingleBranch{TF} <: Branch
    bodies_index::UnitRange{Int64}
    n_branches::Int64
    branch_index::UnitRange{Int64}
    center::SVector{3,TF}   # center of the branch
    radius::TF              # side lengths of the cube encapsulating the branch
    multipole_expansion::Array{TF,3} # multipole expansion coefficients
    local_expansion::Array{TF,3}     # local expansion coefficients
    harmonics::Array{TF,2}
    ML::MArray{Tuple{2,4},TF}
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
    leaf_index::Vector{Int}
    sort_index_list::NTuple{N,Vector{Int}}
    inverse_sort_index_list::NTuple{N,Vector{Int}}
    buffers::TB
    expansion_order::Int64
    n_per_branch::Int64    # max number of bodies in a leaf
    # cost_parameters::MultiCostParameters{N}
    cost_parameters::SVector{N,Float64}
end

struct SingleTree{TF,TB} <: Tree
    branches::Vector{SingleBranch{TF}}        # a vector of `Branch` objects composing the tree
    levels_index::Vector{UnitRange{Int64}}
    leaf_index::Vector{Int}
    sort_index::Vector{Int64}
    inverse_sort_index::Vector{Int64}
    buffer::TB
    expansion_order::Int64
    n_per_branch::Int64    # max number of bodies in a leaf
    # cost_parameters::SingleCostParameters
    cost_parameters::Float64
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
