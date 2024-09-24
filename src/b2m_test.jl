using StaticArrays
using LegendrePolynomials
using Test

include("containers.jl")
include("compatibility.jl")
include("derivativesswitch.jl")
include("tree.jl")
include("harmonics.jl")
include("bodytomultipole.jl")
include("rotate.jl")
include("../test/gravitational_noimport.jl")

#--- point source ---#
@testset "body-to-multipole: point source" begin

x = SVector{3}(0.1,0.2,-0.3)
xs = x + SVector{3}(-0.2,0.07,-0.1)
bodies = [xs[1]; xs[2]; xs[3]; 0.0; 0.7;;]
masses = Gravitational(bodies)
expansion_order = 10
branch = Branch(1:1, 0, 1:0, 0, 1, x, 0.0, expansion_order)

body_to_multipole!(Point{Source{1.0}}, masses, branch, 1:1, branch.harmonics, Val(expansion_order))

Δx = SVector{3}(bodies[1:3]) - branch.center
ρ, θ, ϕ = cartesian_to_spherical(Δx...)

# check against known values
expansion_order_check = 10
for n in 0:expansion_order_check
    for m in -n:-1
        i = harmonic_index(n,-m)
        Mnm = (branch.multipole_expansion[1,1,i] - im*branch.multipole_expansion[2,1,i]) * (-1)^m
        local Rnm = (-1)^n * im^abs(m) * ρ^n * Plm(cos(θ),n,abs(m)) * exp(im*m*ϕ) / factorial(n+abs(m))
        Mnm_check = (-1)^(n+m) * bodies[5,1] * conj(Rnm)
        @test isapprox(Mnm, Mnm_check; atol=1e-12)
    end
    for m in 0:n
        i = harmonic_index(n,m)
        Mnm = branch.multipole_expansion[1,1,i] + im*branch.multipole_expansion[2,1,i]
        local Rnm = (-1)^n * im^abs(m) * ρ^n * Plm(cos(θ),n,abs(m)) * exp(im*m*ϕ) / factorial(n+abs(m))
        Mnm_check = (-1)^(n+m) * bodies[5,1] * conj(Rnm)
        @test isapprox(Mnm, Mnm_check; atol=1e-12)
    end
end

#--- evaluate multipole expansion ---#

include("m2b.jl")

# evaluate the multipole expansion
x_target = SVector{3}(2.3,-4.1, 0.4)
r = x_target - xs
rnorm = sqrt(r'*r)
ϕ_analytic = masses.bodies[1].strength/rnorm
@show branch.multipole_expansion
ϕ_m2b, v_m2b, g_m2b = evaluate_multipole(x_target, branch.center, branch.multipole_expansion, DerivativesSwitch(), Val(expansion_order))

@test isapprox(ϕ_m2b, ϕ_analytic; atol=1e-12)

end

