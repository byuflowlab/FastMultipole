module FastMultipole

#------- IMPORTS -------#

import Base.:^
using LinearAlgebra
using StaticArrays
using WriteVTK

#------- CONSTANTS -------#

const ONE_OVER_4π = 1/4/pi
const ONE_THIRD = 1/3
const π_over_2 = pi/2
const DEBUG = Array{Bool,0}(undef)
DEBUG[] = false

const WARNING_FLAG_LEAF_SIZE = Array{Bool,0}(undef)
WARNING_FLAG_LEAF_SIZE[] = true

const WARNING_FLAG_PMAX = Array{Bool,0}(undef)
WARNING_FLAG_PMAX[] = true

# multithreading parameters
const MIN_NPT_B2M = 100
const MIN_NPT_M2M = 100
const MIN_NPT_M2L = 100
const MIN_NPT_L2L = 100
const MIN_NPT_L2B = 100
const MIN_NPT_NF = 100

# preallocate y-axis rotation matrices by π/2
const Hs_π2 = Float64[1.0]

# preallocate y-axis rotation Wigner matrix normalization
const ζs_mag = Float64[1.0]
const ηs_mag = Float64[1.0]

# preallocate error integrals
const ε_MAX_N = 20
const ε_NX = 100
const ε_Nθ = 40
const ε_Δθ = π_over_2 / ε_Nθ
const ε_Nϕ = 10
const ε_Δϕ = π / ε_Nϕ
const ε_Nω = 100
const ε_Δω = π_over_2 / ε_Nω
const ε_Nγ = 80
const ε_Δγ = π / ε_Nγ

include("error_preintegration.jl")

fpath = joinpath(@__DIR__, "multipole_integrals.csv")
const MULTIPOLE_INTEGRALS = read_write_multipole(fpath)
fpath = joinpath(@__DIR__, "local_integrals.csv")
const LOCAL_INTEGRALS = read_write_local(fpath)

#------- HEADERS AND EXPORTS -------#

include("containers.jl")
include("complex.jl")
include("derivatives.jl")
include("harmonics.jl")
include("rotate.jl")
include("translate.jl")
include("evaluate_expansions.jl")
include("tree.jl")

export Branch, SingleBranch, MultiBranch, Tree, SingleTree, MultiTree, initialize_expansion, initialize_harmonics
export unsort!, resort!, unsorted_index_2_sorted_index, sorted_index_2_unsorted_index

include("compatibility.jl")

export Body, Position, Radius, ScalarPotential, VectorPotential, Velocity, VelocityGradient, Vertex, Normal, Strength
export Vortex, Source, Dipole, SourceDipole, Point, Filament, Panel
export get_n_bodies, buffer_element, body_to_multipole!, direct!, direct_gpu!

include("bodytomultipole.jl")

export body_to_multipole!

include("direct.jl")

export direct!

include("derivativesswitch.jl")

export DerivativesSwitch

include("error.jl")

export multipole_error, local_error, error

include("dynamic_expansion_order.jl")

export get_P

include("sortwrapper.jl")

export SortWrapper

include("interaction_list.jl")

export EqualSpheres, UnequalSpheres, UnequalBoxes, UniformUnequalSpheres, UniformUnequalBoxes, UniformCubes, Dynamic, build_interaction_lists

include("fmm.jl")

export InteractionList, fmm!

include("probe.jl")

export ProbeSystem

include("visualize.jl")

export visualize

#------- PRECALCULATIONS -------#

# precompute y-axis rotation by π/2 matrices up to 20th order
update_Hs_π2!(Hs_π2, Val(20))

# precompute y-axis Wigner matrix normalization up to 20th order
update_ζs_mag!(ζs_mag, Val(20))
update_ηs_mag!(ηs_mag, Val(20))

end # module
