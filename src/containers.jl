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

struct Normal <: Indexable end
const NORMAL = Normal()

struct Strength <: Indexable end
const STRENGTH = Strength()

#####
##### dispatch convenience functions for multipole creation definition
#####
abstract type AbstractKernel{sign} end

struct VortexPoint{sign} <: AbstractKernel{sign} end
VortexPoint(sign=1) = VortexPoint{sign}()
struct VortexLine{sign} <: AbstractKernel{sign} end # not yet derived
VortexLine(sign=1) = VortexLine{sign}()
struct VortexPanel{sign} <: AbstractKernel{sign} end # not yet derived
VortexPanel(sign=1) = VortexPanel{sign}()
struct SourcePoint{sign} <: AbstractKernel{sign} end
SourcePoint(sign=1) = SourcePoint{sign}()
struct UniformSourcePanel{sign} <: AbstractKernel{sign} end
UniformSourcePanel(sign=1) = UniformSourcePanel{sign}()
struct UniformNormalDipolePanel{sign} <: AbstractKernel{sign} end
UniformNormalDipolePanel(sign=1) = UniformNormalDipolePanel{sign}()
struct UniformSourceNormalDipolePanel{sign} <: AbstractKernel{sign} end
UniformSourceNormalDipolePanel(sign=1) = UniformSourceNormalDipolePanel{sign}()

#####
##### dispatch convenience functions to determine which derivatives are desired
#####
"""
    DerivativesSwitch

Switch indicating whether the scalar potential, vector potential, velocity, and/or velocity gradient should be computed for a target system. Information is stored as type parameters, allowing the compiler to compile away if statements.
"""
struct DerivativesSwitch{PS,VPS,VS,GS} end

"""
    DerivativesSwitch(scalar_potential, vector_potential, velocity, velocity_gradient)

Constructs a tuple of [`DerivativesSwitch`](@ref) objects.

# Arguments

- `scalar_potential::Vector{Bool}`: a vector of `::Bool` indicating whether the scalar potential should be computed for each target system
- `vector_potential::Vector{Bool}`: a vector of `::Bool` indicating whether the vector potential should be computed for each target system
- `velocity::Vector{Bool}`: a vector of `::Bool` indicating whether the velocity should be computed for each target system
- `velocity_gradient::Vector{Bool}`: a vector of `::Bool` indicating whether the velocity gradient should be computed for each target system

"""
function DerivativesSwitch(scalar_potential, vector_potential, velocity, velocity_gradient)
    return Tuple(DerivativesSwitch(ps,vps,vs,gs) for (ps,vps,vs,gs) in zip(scalar_potential, vector_potential, velocity, velocity_gradient))
end

"""
    DerivativesSwitch(scalar_potential, vector_potential, velocity, velocity_gradient)

Constructs a single [`DerivativesSwitch`](@ref) object.

# Arguments

- `scalar_potential::Bool`: a `::Bool` indicating whether the scalar potential should be computed for the target system
- `vector_potential::Bool`: a `::Bool` indicating whether the vector potential should be computed for the target system
- `velocity::Bool`: a `::Bool` indicating whether the velocity should be computed for the target system
- `velocity_gradient::Bool`: a `::Bool` indicating whether the velocity gradient should be computed for the target system

"""
function DerivativesSwitch(scalar_potential::Bool, vector_potential::Bool, velocity::Bool, velocity_gradient::Bool)
    return DerivativesSwitch{scalar_potential, vector_potential, velocity, velocity_gradient}()
end

"""
    DerivativesSwitch(scalar_potential, vector_potential, velocity, velocity_gradient, target_systems)

Constructs a `::Tuple` of indentical [`DerivativesSwitch`](@ref) objects of the same length as `target_systems` (if it is a `::Tuple`), or a single [`DerivativesSwitch`](@ref) (if `target_system` is not a `::Tuple`)

# Arguments

- `scalar_potential::Bool`: a `::Bool` indicating whether the scalar potential should be computed for each target system
- `vector_potential::Bool`: a `::Bool` indicating whether the vector potential should be computed for each target system
- `velocity::Bool`: a `::Bool` indicating whether the velocity should be computed for each target system
- `velocity_gradient::Bool`: a `::Bool` indicating whether the velocity gradient should be computed for each target system

"""
function DerivativesSwitch(scalar_potential::Bool, vector_potential::Bool, velocity::Bool, velocity_gradient::Bool, target_systems::Tuple)
    return Tuple(DerivativesSwitch{scalar_potential, vector_potential, velocity, velocity_gradient}() for _ in target_systems)
end

function DerivativesSwitch(scalar_potential, vector_potential, velocity, velocity_gradient, target_systems::Tuple)
    @assert length(scalar_potential) == length(vector_potential) == length(velocity) == length(velocity_gradient) == length(target_systems) "length of inputs to DerivativesSwitch inconsistent"
    return Tuple(DerivativesSwitch{scalar_potential[i], vector_potential[i], velocity[i], velocity_gradient[i]}() for i in eachindex(target_systems))
end

function DerivativesSwitch(scalar_potential::Bool, vector_potential::Bool, velocity::Bool, velocity_gradient::Bool, target_system)
    return DerivativesSwitch{scalar_potential, vector_potential, velocity, velocity_gradient}()
end

DerivativesSwitch() = DerivativesSwitch{true, true, true, true}()

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
    center::SVector{3,TF}   # center of the branch
    radius::TF              # side lengths of the cube encapsulating the branch
    multipole_expansion::Array{TF,3} # multipole expansion coefficients
    local_expansion::Array{TF,3}     # local expansion coefficients
    harmonics::Array{TF,2}
    ML::Matrix{TF}
    lock::ReentrantLock
end

Base.eltype(::MultiBranch{TF}) where TF = TF

struct SingleBranch{TF} <: Branch{TF}
    bodies_index::UnitRange{Int64}
    n_branches::Int64
    branch_index::UnitRange{Int64}
    i_parent::Int64
    center::SVector{3,TF}   # center of the branch
    radius::TF              # side lengths of the cube encapsulating the branch
    multipole_expansion::Array{TF,3} # multipole expansion coefficients
    local_expansion::Array{TF,3}     # local expansion coefficients
    harmonics::Array{TF,2}
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

#####
##### allow input systems to take any form when desired
#####
struct SortWrapper{TS}
    system::TS
    index::Vector{Int}
end

"""
    SortWrapper(system)

Convenience wrapper for systems whose elements cannot be sorted in-place (e.g. structured grids). The resulting object is treated like any other `system`.
"""
function SortWrapper(system)
    return SortWrapper(system,collect(1:get_n_bodies(system)))
end

@inline wrap_duplicates(target_systems::Tuple, source_systems::Tuple) = Tuple(target_system in source_systems ? SortWrapper(target_system) : target_system for target_system in target_systems)

@inline wrap_duplicates(target_system, source_system) = target_system == source_system ? SortWrapper(target_system) : target_system

@inline wrap_duplicates(target_system, source_systems::Tuple) = target_system in source_systems ? SortWrapper(target_system) : target_system

@inline wrap_duplicates(target_systems::Tuple, source_system) = Tuple(target_system == source_system ? SortWrapper(target_system) : target_system for target_system in target_systems)

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
