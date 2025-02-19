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

abstract type SourceVortex <: AbstractKernel end

abstract type Dipole <: AbstractKernel end

abstract type SourceDipole <: AbstractKernel end

abstract type AbstractElement{TK<:AbstractKernel} end

abstract type Point{TK} <: AbstractElement{TK} end

abstract type Filament{TK} <: AbstractElement{TK} end

abstract type Panel{NS,TK} <: AbstractElement{TK} end

#------- dispatch convenience functions to determine which derivatives are desired -------#

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

#------- error predictors -------#

abstract type ErrorMethod end

struct UnequalSpheres <: ErrorMethod end

struct UnequalBoxes <: ErrorMethod end

struct UniformUnequalSpheres <: ErrorMethod end

struct UniformUnequalBoxes <: ErrorMethod end

struct RotatedCoefficients <: ErrorMethod end

#------- octree creation -------#

"""
    Branch{TF,N}

Branch object used to sort more than one system into an octree. Type parameters represent:

* `TF`: the floating point type (would be a dual number if using algorithmic differentiation)
* `N`: the number of systems represented

**Fields**

* `bodies_index::Vector{UnitRange}`: vector of unit ranges indicating the index of bodies in each represented system, respectively
* `n_branches::Int`: number of child branches corresponding to this branch
* `branch_index::UnitRange`: indices of this branch's child branches
* `i_parent::Int`: index of this branch's parent
* `i_leaf::Int`: if this branch is a leaf, what is its index in its parent `<:Tree`'s `leaf_index` field
* `center::Vector{TF}`: center of this branch at which its multipole and local expansions are centered
* `source_radius::TF`: if this branch is a leaf, distance from `center` to the outer edge of farthest body contained in this branch; otherwise, this is the distance from `center` to the corner of its `source_box`
* `target_radius::TF`: distance from `center` to the farthest body center contained in this branch
* `source_box::Vector{TF}`: vector of length 6 containing the distances from the center to faces of a rectangular prism completely enclosing all bodies with their finite radius in the negative x, positive x, negative y, positive y, negative z, and positive z directions, respectively
* `target_box::Vector{TF}`: vector of length 3 containing the distances from the center to faces of a rectangular prism completely enclosing all body centers in the x, y, and z direction, respectively
* multipole_expansion::Array{TF,3}`: array of size (2,2,) containing the multipole expansion coefficients; the first index indicates real or imaginary, the second index indicates scalar potential or the second component of the Lamb-Helmholtz decomposition, and the third index `k` indicates the expansion coefficient of degree \$n\$ and order \$m\$, as \$k = p(p+1)/2 + m + 1\$
* local_expansion::Array{TF,3}`: array of size `(2,2,(p+1)(p+2)/2)` containing the local expansion coefficients; the first index indicates real or imaginary, the second index indicates scalar potential or the second component of the Lamb-Helmholtz decomposition, and the third index `k` indicates the expansion coefficient of degree \$n\$ and order \$m\$, as \$k = p(p+1)/2 + m + 1\$
* `harmonics::Array{TF,3}`: array of size `(2,2,(p+1)(p+2)/2)` used as storage for regular harmonics, irregular harmonics, or whatever is needed, and is indexed as multipole and local expansions
* `lock::ReentrantLock`: lock used to avoid data race conditions when modifying this branch or its corresponding bodies

"""
struct Branch{TF,N}
    n_bodies::SVector{N,Int64}
    bodies_index::SVector{N,UnitRange{Int64}}
    n_branches::Int64
    branch_index::UnitRange{Int64}
    i_parent::Int64
    i_leaf::Int64
    source_center::SVector{3,TF}   # center of the branch
    target_center::SVector{3,TF}   # center of the branch
    source_radius::TF
    target_radius::TF
    source_box::SVector{3,TF} # x, y, and z half widths of the box encapsulating all sources
    target_box::SVector{3,TF} # x, y, and z half widths of the box encapsulating all sources
    # multipole_expansion::Array{TF,3} # multipole expansion coefficients
    # local_expansion::Array{TF,3}     # local expansion coefficients
    # harmonics::Array{TF,3}
    lock::ReentrantLock
end

function Branch(n_bodies, bodies_index, n_branches, branch_index, i_parent, i_leaf, source_center, target_center, source_radius::TF, target_radius, source_box, target_box, multipole_expansion, local_expansion, harmonics, lock) where TF
    Branch(n_bodies, bodies_index, n_branches, branch_index, i_parent, i_leaf, source_center, target_center, source_radius, target_radius, source_box, target_box, error, multipole_expansion, local_expansion, harmonics, lock)
end

Base.eltype(::Branch{TF,<:Any}) where TF = TF

"""
bodies[index_list] is the same sort operation as performed by the tree
sorted_bodies[inverse_index_list] undoes the sort operation performed by the tree
"""
struct Tree{TF,N,TB}
    branches::Vector{Branch{TF,N}}        # a vector of `Branch` objects composing the tree
    expansions::Array{TF,4}
    levels_index::Vector{UnitRange{Int64}}
    leaf_index::Vector{Int}
    sort_index_list::NTuple{N,Vector{Int}}
    inverse_sort_index_list::NTuple{N,Vector{Int}}
    buffers::TB
    expansion_order::Int64
    leaf_size::SVector{N,Int64}    # max number of bodies in a leaf
    # cost_parameters::MultiCostParameters{N}
    # cost_parameters::SVector{N,Float64}
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
