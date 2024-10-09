using LegendrePolynomials
using LinearAlgebra
using StaticArrays
using Test

const ONE_OVER_4π = 1/4/pi

# include headers
include("containers.jl")
include("derivativesswitch.jl")
include("harmonics.jl")
include("rotate.jl")
include("bodytomultipole.jl")
include("translate.jl")
include("tree.jl")
include("evaluate_expansions.jl")
include("../test/gravitational_noimport.jl")

# run tests
include("harmonics_test.jl")
include("rotate_test.jl")
include("b2m_test.jl")
include("m2m_test.jl")
include("m2l_test.Jl")
include("l2l_test.jl")
include("evaluate_expansions_test.jl")
