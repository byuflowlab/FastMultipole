module FastMultipole

#------- IMPORTS -------#

import Base.:^
using LinearAlgebra
using StaticArrays
using WriteVTK

#------- CONSTANTS -------#

const ONE_OVER_4π = 1/4/pi
const ONE_THIRD = 1/3
const DEBUG = Array{Bool,0}(undef)
DEBUG[] = false

# multithreading parameters
#const MIN_NPT_B2M = 1
#const MIN_NPT_M2M = 1
#const MIN_NPT_M2L = 1
#const MIN_NPT_L2L = 1
#const MIN_NPT_L2B = 1
#const MIN_NPT_NF = 1
const MIN_NPT_B2M = 100
const MIN_NPT_M2M = 100
const MIN_NPT_M2L = 100
const MIN_NPT_L2L = 100
const MIN_NPT_L2B = 100
const MIN_NPT_NF = 100

# preallocate y-axis rotation matrices by π/2
const Hs_π2 = Float64[1.0]

# preallocate y-axis Wigner matrix normalization
const ζs_mag = Float64[1.0]
const ηs_mag = Float64[1.0]

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

include("sortwrapper.jl")

export SortWrapper

include("interaction_list.jl")

export EqualSpheres, UnequalSpheres, UnequalBoxes, Dynamic, build_interaction_lists

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
