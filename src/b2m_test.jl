using StaticArrays
using LegendrePolynomials
using Test

include("containers.jl")
include("compatibility.jl")
include("tree.jl")
include("harmonics.jl")
include("b2m.jl")
include("rotate.jl")
include("../test/gravitational_noimport.jl")

#--- point source ---#

x = SVector{3}(0.1,0.2,-0.3)
xs = x + SVector{3}(-0.2,0.07,-0.1)
bodies = [xs[1]; xs[2]; xs[3]; 0.0; 0.7;;]
masses = Gravitational(bodies)
expansion_order = 10
branch = Branch(1:1, 0, 1:0, 0, 1, x, 0.0, expansion_order)

B2M!(Point{Source{1.0}}, masses, branch, 1:1, branch.harmonics, Val(expansion_order))

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
ϕ_m2b, v_m2b, g_m2b = evaluate_multipole(x_target, branch.center, branch.multipole_expansion, DerivativesSwitch(), Val(expansion_order))

@test isapprox(ϕ_m2b, ϕ_analytic; atol=1e-12)

#--- translate multipole ---#

include("translate.jl")

# preallocate containers
Hs_π2 = [1.0]
update_Hs_π2!(Hs_π2, expansion_order)
Ts = zeros(length_Ts(expansion_order))
eimϕs = zeros(2, expansion_order+1)
weights_tmp_1 = initialize_expansion(expansion_order, eltype(Ts))
weights_tmp_2 = initialize_expansion(expansion_order, eltype(Ts))
translated_weights = initialize_expansion(expansion_order, eltype(Ts))

# normalization
global ζs_mag = zeros(length_ζs(expansion_order))
update_ζs_mag!(ζs_mag, 0, expansion_order)

# next multipole branch
branch_2 = Branch(2:2, 0, 1:0, 0, 1, x + SVector{3}(0.1, 0.2, 0.14), 0.0, expansion_order)

multipole_to_multipole!(branch_2, branch, weights_tmp_1, weights_tmp_2, Ts, eimϕs, Hs_π2, expansion_order)


