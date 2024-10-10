#------- dispatch for common interface for external packages -------#

abstract type Indexable end

struct Body <: Indexable end
const BODY = Body()

struct Position <: Indexable end
const POSITION = Position()

struct Radius <: Indexable end
const RADIUS = Radius()

struct ScalarPotential <: Indexable end
const SCALAR_POTENTIAL = ScalarPotential()

struct Velocity <: Indexable end
const VELOCITY = Velocity()

struct VelocityGradient <: Indexable end
const VELOCITY_GRADIENT = VelocityGradient()

struct Vertex <: Indexable end
const VERTEX = Vertex()

struct Normal <: Indexable end
const NORMAL = Normal()

struct Strength <: Indexable end
const STRENGTH = Strength()

#------- dispatch convenience functions for multipole creation definition -------#

abstract type AbstractKernel end

abstract type Vortex <: AbstractKernel end

abstract type Source <: AbstractKernel end

abstract type Dipole <: AbstractKernel end

abstract type SourceDipole <: AbstractKernel end

abstract type AbstractElement{TK<:AbstractKernel} end

abstract type Point{TK} <: AbstractElement{TK} end

abstract type Filament{TK} <: AbstractElement{TK} end

abstract type Panel{NS,TK} <: AbstractElement{TK} end

#####
##### dispatch convenience functions to determine which derivatives are desired
#####
"""
    DerivativesSwitch

Switch indicating whether the scalar potential, vector potential, velocity, and/or velocity gradient should be computed for a target system. Information is stored as type parameters, allowing the compiler to compile away if statements.
"""
struct DerivativesSwitch{PS,VS,GS} end

"""
    ExpansionSwitch

Switch indicating which expansions should be used:

1. scalar potential (`SP`)
2. vector potential via Lamb-Helmholtz decomposition (`VP`)
"""
struct ExpansionSwitch{SP,VP} end

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
abstract type Branch{TF} end

struct MultiBranch{TF,N} <: Branch{TF}
    bodies_index::SVector{N,UnitRange{Int64}}
    n_branches::Int64
    branch_index::UnitRange{Int64}
    i_parent::Int64
    i_leaf::Int64
    center::SVector{3,TF}   # center of the branch
    radius::TF              # side lengths of the cube encapsulating the branch
    multipole_expansion::Array{TF,3} # multipole expansion coefficients
    local_expansion::Array{TF,3}     # local expansion coefficients
    harmonics::Array{TF,3}
    ML::Matrix{TF}
    lock::ReentrantLock
end

Base.eltype(::MultiBranch{TF}) where TF = TF

struct SingleBranch{TF} <: Branch{TF}
    bodies_index::UnitRange{Int64}
    n_branches::Int64
    branch_index::UnitRange{Int64}
    i_parent::Int64
    i_leaf::Int64
    center::SVector{3,TF}   # center of the branch
    radius::TF              # side lengths of the cube encapsulating the branch
    multipole_expansion::Array{TF,3} # multipole expansion coefficients
    local_expansion::Array{TF,3}     # local expansion coefficients
    harmonics::Array{TF,3}
    #expansion_storage::Array{TF,3}
    ML::Matrix{TF}
    lock::ReentrantLock
end

Base.eltype(::SingleBranch{TF}) where TF = TF

"""
    abstract type Tree{TF,P} end

Supertype of all octree structures with `TF` the floating point type and `P` the expansion order.
"""
abstract type Tree{TF,P} end

"""
bodies[index_list] is the same sort operation as performed by the tree
sorted_bodies[inverse_index_list] undoes the sort operation performed by the tree
"""
struct MultiTree{TF,N,TB,P} <: Tree{TF,P}
    branches::Vector{MultiBranch{TF,N}}        # a vector of `Branch` objects composing the tree
    levels_index::Vector{UnitRange{Int64}}
    leaf_index::Vector{Int}
    sort_index_list::NTuple{N,Vector{Int}}
    inverse_sort_index_list::NTuple{N,Vector{Int}}
    buffers::TB
    expansion_order::Val{P}
    leaf_size::Int64    # max number of bodies in a leaf
    # cost_parameters::MultiCostParameters{N}
    # cost_parameters::SVector{N,Float64}
end

struct SingleTree{TF,TB,P} <: Tree{TF,P}
    branches::Vector{SingleBranch{TF}}        # a vector of `Branch` objects composing the tree
    levels_index::Vector{UnitRange{Int64}}
    leaf_index::Vector{Int}
    sort_index::Vector{Int64}
    inverse_sort_index::Vector{Int64}
    buffer::TB
    expansion_order::Val{P}
    leaf_size::Int64    # max number of bodies in a leaf
    # cost_parameters::SingleCostParameters
    # cost_parameters::Float64
end

struct InteractionList{TF}
    influence_matrices::Vector{Matrix{TF}}
    strengths::Vector{TF}
    influence::Vector{TF}
    direct_list::Vector{SVector{2,Int32}}
end

Base.length(list::InteractionList) = length(list.direct_list)

#####
##### allow input systems to take any form when desired
#####
struct SortWrapper{TS}
    system::TS
    index::Vector{Int}
end

#####
##### when we desire to evaluate the potential at locations not coincident with source centers
#####

"""
    ProbeSystem

Convenience system for defining locations at which the potential, velocity, or velocity gradient may be desired.
"""
struct ProbeSystem{TF,TSP,TVP,TV,TVG}
    position::Vector{SVector{3,TF}}
    scalar_potential::TSP
    vector_potential::TVP
    velocity::TV
    velocity_gradient::TVG
end
