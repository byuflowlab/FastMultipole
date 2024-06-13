using Test
import Statistics
S = Statistics

import FastMultipole

using FastMultipole.LinearAlgebra
using FastMultipole.StaticArrays
using ForwardDiff
using LegendrePolynomials
using Random
using SpecialFunctions

test_dir = @__DIR__

#####
##### define gravitational kernel and mass elements
#####

include(joinpath(test_dir, "gravitational.jl"))

# allow expansion_order::Int to be used while testing
FastMultipole.M2M!(branch, child, harmonics, M, expansion_order::Int) = FastMultipole.M2M!(branch, child, harmonics, M, Val(expansion_order))
FastMultipole.M2L!(target_branch, source_branch, harmonics, L, expansion_order::Int) = FastMultipole.M2L!(target_branch, source_branch, harmonics, L, Val(expansion_order))
FastMultipole.upward_pass_singlethread!(branches, systems, expansion_order::Int) = FastMultipole.upward_pass_singlethread!(branches, systems, Val(expansion_order))
FastMultipole.M2L!(target_branch, source_branch, expansion_order::Int) = FastMultipole.M2L!(target_branch, source_branch, Val(expansion_order))

@testset "complex" begin
    z1 = rand(Complex{Float64})
    z2 = rand(Complex{Float64})
    z3 = rand(Complex{Float64})

    # addition
    @test real(z1+z2) ≈ FastMultipole.complex_add(real(z1), imag(z1), real(z2), imag(z2))[1]
    @test imag(z1+z2) ≈ FastMultipole.complex_add(real(z1), imag(z1), real(z2), imag(z2))[2]

    # subtraction
    @test real(z1-z2) ≈ FastMultipole.complex_subtract(real(z1), imag(z1), real(z2), imag(z2))[1]
    @test imag(z1-z2) ≈ FastMultipole.complex_subtract(real(z1), imag(z1), real(z2), imag(z2))[2]

    # multiplication
    @test real(z1*z2) ≈ FastMultipole.complex_multiply(real(z1), imag(z1), real(z2), imag(z2))[1]
    @test imag(z1*z2) ≈ FastMultipole.complex_multiply(real(z1), imag(z1), real(z2), imag(z2))[2]

    @test real(z1*z2*z3) ≈ FastMultipole.complex_multiply(real(z1), imag(z1), real(z2), imag(z2), real(z3), imag(z3))[1]
    @test imag(z1*z2*z3) ≈ FastMultipole.complex_multiply(real(z1), imag(z1), real(z2), imag(z2), real(z3), imag(z3))[2]

    # division
    @test real(z1/z2) ≈ FastMultipole.complex_divide(real(z1), imag(z1), real(z2), imag(z2))[1]
    @test imag(z1/z2) ≈ FastMultipole.complex_divide(real(z1), imag(z1), real(z2), imag(z2))[2]
    @test real(z1/z2) ≈ FastMultipole.complex_divide_real(real(z1), imag(z1), real(z2), imag(z2))
    @test imag(z1/z2) ≈ FastMultipole.complex_divide_imag(real(z1), imag(z1), real(z2), imag(z2))

	# cross product
	z1_vec = SVector{3}(rand(Complex{Float64}) for _ in 1:3)
	z2_vec = SVector{3}(rand(Complex{Float64}) for _ in 1:3)
	@test real.(cross(z1_vec,z2_vec)) ≈ FastMultipole.complex_cross_real(real(z1_vec[1]), imag(z1_vec[1]), real(z1_vec[2]), imag(z1_vec[2]), real(z1_vec[3]), imag(z1_vec[3]), real(z2_vec[1]), imag(z2_vec[1]), real(z2_vec[2]), imag(z2_vec[2]), real(z2_vec[3]), imag(z2_vec[3]))
end

@testset "direct" begin

    function V(xi, xj, mj; G=1)
        Rho_ij = xi - xj
        rho_ij = sqrt(Rho_ij' * Rho_ij)
        vij = G * mj / rho_ij / 4 / pi
        if isinf(vij); vij = 0.0; end
        return Rho_ij, rho_ij, vij
    end

    x = [
        -5.0 -4.8 -5.1 4.9 5.2 5.1 5.3;
        0.4 0.2 0.1 0.3 0.2 -0.1 0.1;
        -0.1 0.2 0.1 0.1 0.0 -0.1 -0.2
    ]

    m = [
        1.4,
        2.0,
        2.3,
        4.1,
        1.1,
        3.4,
        4.5
    ]

    Rho_ijs = zeros(3,length(m),length(m))
    rho_ijs = zeros(length(m),length(m))
    V_ijs = zeros(length(m),length(m))

    for i in 1:length(m)
        for j in 1:length(m)
            Rho_ij, rho_ij, vij = V(x[:,i], x[:,j], m[j])
            Rho_ijs[:,i,j] .= Rho_ij
            rho_ijs[i,j] = rho_ij
            V_ijs[i,j] = vij
        end
    end

    V_tots = zeros(length(m))
    for i = 1:length(m)
        V_tots[i] = sum(V_ijs[i,:])
    end

    bodies = vcat(x,rand(1,length(m)),m',zeros(3,length(m)))
    mass = Gravitational(bodies)

    FastMultipole.direct!(mass)
    V_tots_direct = mass.potential[1,:]

    for i in 1:length(V_tots)
        @test isapprox(V_tots[i], V_tots_direct[i]; atol=1e-4)
    end
end

@testset "tree" begin

    # build list of elements to sort
    xs = [
        1.2 1.1 0.8;
        0.8 0.9 0.2;
        0.1 0.2 0.9;
        0.1 0.3 0.2;
        0.2 0.25 0.4
    ]
    radii = zeros(Float64,size(xs)[1])
    ms = rand(size(xs)[1])
    bodies = vcat(xs',radii',ms',zeros(3,length(ms)))
    elements = Gravitational(bodies)

    # test center_radius function
    center, radius = FastMultipole.center_radius((elements,); scale_radius = 1.00001)
    test_center = [0.65, 0.65, 0.55]
    test_radius = 0.5500055

    for i in 1:3
        @test isapprox(center[i], test_center[i]; atol=1e-4)
    end
    @test isapprox(radius, test_radius; atol=1e-4)

    # test branch! function
    expansion_order, leaf_size, multipole_threshold = 2, 1, 0.5
    tree = FastMultipole.Tree((elements,); expansion_order, leaf_size, shrink_recenter=false)

    test_branches = [
        5 0.65 0.65 0.55 0.5500055;
        2 0.37499725 0.37499725 0.27499725 0.27500275;
        1 0.92500275 0.92500275 0.27499725 0.27500275;
        1 0.37499725 0.37499725 0.82500275 0.27500275;
        1 0.92500275 0.92500275 0.82500275 0.27500275;
        1 0.237495875 0.237495875 0.137495875 0.137501375;
        1 0.237495875 0.237495875 0.412498625 0.137501375;
    ]

    @test length(tree.branches) == size(test_branches)[1]

    for i_branch in 1:length(tree.branches)
        @test isapprox(length(tree.branches[i_branch].bodies_index[1]), test_branches[i_branch,1]; atol=1e-8)
        for i in 1:3
            @test isapprox(tree.branches[i_branch].center[i], test_branches[i_branch,1+i]; atol=1e-7)
        end
        @test isapprox(tree.branches[i_branch].radius, test_branches[i_branch,5]; atol=1e-7)
    end
end

@testset "cartesian to spherical" begin
    # cartesian to spherical
    rho = 1.0
    theta = pi/4
    phi = pi/2
    that = [rho, theta, phi]
    x = rho * sin(theta) * cos(phi)
    y = rho * sin(theta) * sin(phi)
    z = rho * cos(theta)
    this = [x,y,z]
    FastMultipole.cartesian_2_spherical!(this)
    for i in 1:3
        @test isapprox(this[i], that[i]; atol=1e-10)
    end
end


new_order = [4,5,2,3,1]

const new_order_index = [5,3,4,1,2]

# get the new index of mass_i as new_order_index[mass_i]

@testset "spherical B2M" begin

xs = [
    1.2 1.1 0.8;
    0.8 0.9 0.2;
    0.1 0.2 0.9;
    0.1 0.3 0.2;
    0.2 0.25 0.4
]

ms = [
    0.8,
    1.1,
    2.2,
    0.5,
    1.9
]

bodies = zeros(8,length(ms))
for i in 1:length(ms)
    bodies[1:3,i] .= xs[i,1:3]
    bodies[4,i] = 0.0
    bodies[5,i] = ms[i]
end
system = Gravitational(bodies)

expansion_order = 2
leaf_size = 1
multipole_threshold = 0.25
tree = FastMultipole.Tree((system,); expansion_order, leaf_size, shrink_recenter=false)

i_mass = 1
i_branch = 5 # use the first mass

harmonics = zeros(Float64,2,(expansion_order+1)^2)

FastMultipole.B2M!(system, tree.branches[i_branch], new_order_index[i_mass]:new_order_index[i_mass], harmonics, tree.expansion_order)

center = tree.branches[i_branch].center

x_target = [10.1,-7.3,8.6]
target_potential = zeros(4)
FastMultipole.M2B!(target_potential, x_target, i_branch, tree)

u_fmm = target_potential[1]

dx = x_target - xs[1,:]
u_check = ms[1] / sqrt(dx' * dx) * ONE_OVER_4PI

function Ylm(theta, phi, l, m)
    ylm = sqrt(factorial(big(l-abs(m)))/ factorial(big(l+abs(m)))) * Plm(cos(theta), l, abs(m)) * exp(im * m * phi)
end

function evaluate_biot_savart(x_source, x_target, q_source, P)
    v = 0.0
    i = 1

    for l in 0:P
        for m in -l:l

            v += q_source * x_source[1]^l / x_target[1]^(l+1) * real(Ylm(x_target[2], x_target[3], l, m) * conj(Ylm(x_source[2], x_source[3], l, m))) * ONE_OVER_4PI
            i += 1
        end
    end
    return v
end

x_source_sph = FastMultipole.cartesian_2_spherical(xs[1,:] - center)
x_target_sph = FastMultipole.cartesian_2_spherical(x_target - center)
u_check_man = evaluate_biot_savart(x_source_sph, x_target_sph, ms[1], expansion_order);
@test isapprox(u_check_man, u_check; atol=1e-6)
@test isapprox(u_check, u_fmm; atol=1e-6)

end

@testset "spherical M2M" begin
xs = [
    1.2 1.1 0.8;
    0.8 0.9 0.2;
    0.1 0.2 0.9;
    0.1 0.3 0.2;
    0.2 0.25 0.4
]

ms = [
    0.8,
    1.1,
    2.2,
    0.5,
    1.9
]

bodies = zeros(8,length(ms))
for i in 1:length(ms)
    bodies[1:3,i] .= xs[i,1:3]
    bodies[5,i] = ms[i]
end
elements = Gravitational(bodies)

expansion_order = 3
tree = FastMultipole.Tree((elements,); expansion_order, leaf_size=1, shrink_recenter=false)

i_branch = 2 # contains 4th and 5th elements
i_branch_4 = 6 # use the fourth mass
# i_branch_5 = 7 # use the fifth mass
harmonics = zeros(Float64,2,(expansion_order+1)^2)
M = zeros(Float64,2,4)
# using only the 4th mass: (note it has been reordered)
FastMultipole.B2M!(elements, tree.branches[i_branch_4], new_order_index[4]:new_order_index[4], harmonics, tree.expansion_order) # evaluate multipole coefficients
for child_branch in view(tree.branches, tree.branches[i_branch].branch_index) # translate coefficients to the center of branch 2
    FastMultipole.M2M!(tree.branches[i_branch], child_branch, harmonics, M, expansion_order)
end

x_target = [8.3,1.4,-4.2]
target_potential = zeros(4)
target = x_target
FastMultipole.M2B!(target_potential, target, i_branch, tree)
u_fmm = target_potential[1]

target_potential .*= 0
FastMultipole.M2B!(target_potential, target, i_branch_4, tree)
u_fmm_no_x = target_potential[1]

dx = x_target - xs[4,:]
u_check = ms[4] / sqrt(dx'*dx) * ONE_OVER_4PI

@test isapprox(u_fmm, u_fmm_no_x; atol=1e-5)
@test isapprox(u_fmm, u_check; atol=1e-5)

end

@testset "spherical L2P" begin
xs = [
    1.2 1.1 0.8;
    0.8 0.9 0.2;
    0.1 0.2 0.9;
    0.1 0.3 0.2;
    0.2 0.25 0.4
]

ms = [
    0.8,
    1.1,
    2.2,
    0.5,
    1.9
]

bodies = zeros(8,length(ms))
for i in 1:length(ms)
    bodies[1:3,i] .= xs[i,1:3]
    bodies[5,i] = ms[i]
end
elements = Gravitational(bodies)

expansion_order = 20
tree = FastMultipole.Tree((elements,); expansion_order, leaf_size=1, shrink_recenter=false)

branch_i = 2 # contains two elements; 4 and 5
target_i = new_order_index[4]
source_i = new_order_index[1] # just needs to be farther away than the target to ensure convergence

dx_source = FastMultipole.cartesian_2_spherical(elements[source_i,FastMultipole.POSITION] - tree.branches[branch_i].center)
dx_target = FastMultipole.cartesian_2_spherical(elements[target_i,FastMultipole.POSITION] - tree.branches[branch_i].center)

local_coefficients_theta = zeros(Float64, 2, ((expansion_order+1)*(expansion_order+2))>>1)
local_coefficients_expanded = zeros(Float64, 2, (expansion_order+1)^2)
local_coefficients_expanded_theta = zeros(Float64, 2, (expansion_order+1)^2)
FastMultipole.irregular_harmonic!(local_coefficients_expanded, dx_source..., expansion_order)
local_coefficients_expanded .*= ms[1]
regular_harmonics_expanded = zeros(Float64, 2, (expansion_order+1)^2)
FastMultipole.regular_harmonic!(regular_harmonics_expanded, dx_target..., expansion_order)

FastMultipole.B2L!(tree, branch_i, elements[source_i,FastMultipole.POSITION], elements.bodies[source_i].strength)

FastMultipole.L2B!(elements, target_i:target_i, tree.branches[branch_i].local_expansion, FastMultipole.DerivativesSwitch(), Val(expansion_order), tree.branches[branch_i].center)
u_fmm = elements.potential[1,target_i] * ONE_OVER_4PI

dx_direct = xs[4,:] - xs[1,:]
u_check = 1 / sqrt(dx_direct' * dx_direct) * ONE_OVER_4PI
u_check *= ms[1]

u_man = sum(regular_harmonics_expanded[1,:].*local_coefficients_expanded[1,:] .+ regular_harmonics_expanded[2,:].*local_coefficients_expanded[2,:]) * ONE_OVER_4PI
#u_man = real(sum(regular_harmonics_expanded' * local_coefficients_expanded)) # appears to work

@test isapprox(u_check, u_fmm; atol=1e-12)
@test isapprox(u_check, u_man; atol=1e-12)

end

@testset "spherical L2L" begin
xs = [
    1.2 1.1 0.8;
    0.8 0.9 0.2;
    0.1 0.2 0.9;
    0.1 0.3 0.2;
    0.2 0.25 0.4
]

ms = [
    0.8,
    1.1,
    2.2,
    0.5,
    1.9
]

bodies = zeros(8,length(ms))
for i in 1:length(ms)
    bodies[1:3,i] .= xs[i,1:3]
    bodies[5,i] = ms[i]
end
elements = Gravitational(bodies)

expansion_order = 20
tree = FastMultipole.Tree((elements,); expansion_order, leaf_size=1, shrink_recenter=false)

# local coefficient at branch 2 due to mass 1
FastMultipole.B2L!(tree, 2, elements[new_order_index[1],FastMultipole.POSITION], elements.bodies[new_order_index[1]].strength)
# local_2 = deepcopy(tree.branches[2].local_expansion)

# check L2P now:
harmonics = zeros(Float64,2,(expansion_order+1)^2)
harmonics_theta = zeros(Float64,2,(expansion_order+1)^2)
harmonics_theta_2 = zeros(Float64,2,(expansion_order+1)^2)
workspace = zeros(3,4)
spherical_potential = zeros(52)
FastMultipole.L2B!(elements, new_order_index[5]:new_order_index[5], tree.branches[2].local_expansion, FastMultipole.DerivativesSwitch(), Val(expansion_order), tree.branches[2].center)
u_fmm_no_x = elements.potential[1,new_order_index[5]] * ONE_OVER_4PI
elements.potential[1,new_order_index[5]] *= 0

# translate local expansion to branch 7 (mass 5)
FastMultipole.L2L!(tree.branches[2], tree.branches[7], harmonics, zeros(eltype(tree.branches[1].multipole_expansion),2,4), tree.expansion_order)

local_coefficients_check = zeros(Float64,2, (expansion_order+1)^2)
dx_check, dy_check, dz_check = FastMultipole.cartesian_2_spherical(elements[new_order_index[1],FastMultipole.POSITION] - tree.branches[7].center)
FastMultipole.irregular_harmonic!(local_coefficients_check, dx_check, dy_check, dz_check, expansion_order)
local_coefficients_check .*= ms[1]

# evaluate local expansion at mass 5
FastMultipole.L2B!((elements,), tree.branches[7], (FastMultipole.DerivativesSwitch(),), Val(expansion_order))
u_fmm = elements.potential[1,new_order_index[5]] * ONE_OVER_4PI

dx_direct = elements[new_order_index[5],FastMultipole.POSITION] - elements[new_order_index[1],FastMultipole.POSITION]
u_check = ms[1] / sqrt(dx_direct' * dx_direct) * ONE_OVER_4PI

regular_harmonics = zeros(Float64,2, (expansion_order+1)^2)
dx_target = FastMultipole.cartesian_2_spherical(elements[new_order_index[5],FastMultipole.POSITION] - tree.branches[7].center)
FastMultipole.regular_harmonic!(regular_harmonics, dx_target..., expansion_order)
#u_man = real(sum(regular_harmonics' * local_coefficients_check))
u_man = sum(regular_harmonics[1,:].*local_coefficients_check[1,:] .+ regular_harmonics[2,:].*local_coefficients_check[2,:]) * ONE_OVER_4PI

@test isapprox(u_check, u_man; atol=1e-12)
@test isapprox(u_check, u_fmm_no_x; atol=1e-12)
@test isapprox(u_check, u_fmm; atol=1e-12)

end

@testset "spherical M2L" begin
xs = [
    1.2 1.1 0.8;
    0.8 0.9 0.2;
    0.1 0.2 0.9;
    0.1 0.3 0.2;
    0.2 0.25 0.4
]

ms = [
    0.8,
    1.1,
    2.2,
    0.5,
    1.9
]

bodies = zeros(8,length(ms))
for i in 1:length(ms)
    bodies[1:3,i] .= xs[i,1:3]
    bodies[5,i] = ms[i]
end
elements = Gravitational(bodies)

expansion_order = 30
tree = FastMultipole.Tree((elements,); expansion_order, leaf_size=1, shrink_recenter=false)

i_branch_multipole = 7 # mass 5
i_branch_local = 5 # mass 1
harmonics = zeros(Float64,2, (expansion_order+1)^2)
harmonics_theta = zeros(Float64,2, (expansion_order+1)^2)
harmonics_theta_2 = zeros(Float64,2, (expansion_order+1)^2)
workspace = zeros(3,4)

FastMultipole.B2M!(elements, tree.branches[i_branch_multipole], new_order_index[5]:new_order_index[5], harmonics, tree.expansion_order)

# test Multipole # checks out
dx_mp = xs[5,:] - tree.branches[i_branch_multipole].center
dx_mp = FastMultipole.cartesian_2_spherical(dx_mp)
FastMultipole.regular_harmonic!(harmonics, dx_mp..., expansion_order)

these_harmonics = harmonics[1,:] + im .* harmonics[2,:]
multipole_check =  these_harmonics * ms[5] * ONE_OVER_4PI
dx_mp = xs[1,:] - tree.branches[i_branch_multipole].center
dx_mp = FastMultipole.cartesian_2_spherical(dx_mp)
FastMultipole.irregular_harmonic!(harmonics, dx_mp..., expansion_order)
those_harmonics = harmonics[1,:] + im .* harmonics[2,:]
u_check_mp = real(sum(those_harmonics' * multipole_check))
# reformat_expansion = tree.branches[i_branch_multipole].multipole_expansion[1,1,:] .+ im .* tree.branches[i_branch_multipole].multipole_expansion[2,1,:]
# @show tree.branches[i_branch_multipole].multipole_expansion[:,1,1:10] multipole_check[1:10]
# @show size(reformat_expansion) size(multipole_check)
# for i in eachindex(multipole_check)
# 	@test isapprox(reformat_expansion[i], multipole_check[i]; atol=1e-12)
# end

m2l_harmonics = zeros(eltype(tree.branches[1].multipole_expansion), 2, (expansion_order<<1 + 1)*(expansion_order<<1 + 1))
L = zeros(eltype(tree.branches[1].local_expansion), 2, 4)
FastMultipole.M2L!(tree.branches[i_branch_local], tree.branches[i_branch_multipole], m2l_harmonics, L, expansion_order)
FastMultipole.L2B!(elements, new_order_index[1]:new_order_index[1], tree.branches[i_branch_local].local_expansion, FastMultipole.DerivativesSwitch(), Val(expansion_order), tree.branches[i_branch_local].center)
u_fmm = elements.potential[1,new_order_index[1]]

local_exp = tree.branches[i_branch_local].local_expansion[1]

dx_direct = elements[new_order_index[1],FastMultipole.POSITION] - elements[new_order_index[5], FastMultipole.POSITION]
u_direct = elements.bodies[new_order_index[5]].strength[1] / sqrt(dx_direct' * dx_direct) * ONE_OVER_4PI

@test isapprox(u_direct, u_check_mp; atol=1e-12)
@test isapprox(u_fmm, u_direct; atol=1e-12)

end

@testset "fmm" begin

xs = [
    1.2 1.1 0.8;
    0.8 0.9 0.2;
    0.1 0.2 0.9;
    0.1 0.3 0.2;
    0.2 0.25 0.4
]

ms = [
    0.8,
    1.1,
    2.2,
    0.5,
    1.9
]

bodies = zeros(8,length(ms))
for i in 1:length(ms)
    bodies[1:3,i] .= xs[i,1:3]
    bodies[5,i] = ms[i]
end
elements = Gravitational(bodies)

expansion_order = 24
multipole_threshold = 0.5
tree = FastMultipole.Tree((elements,); expansion_order, leaf_size=1, shrink_recenter=false)

# perform upward pass
FastMultipole.upward_pass_singlethread!(tree.branches, (elements,), expansion_order)

# m6 = tree.branches[6].multipole_expansion
target = [4.1,2.2,3.4]
dx_direct_6 = target - elements[new_order_index[4],FastMultipole.POSITION]
u_direct_6 = ms[4] / sqrt(dx_direct_6' * dx_direct_6) * ONE_OVER_4PI

mass_target_potential = zeros(4)
mass_target = target
FastMultipole.M2B!(mass_target_potential, mass_target, 6, tree)
u_fmm_6 = mass_target_potential[1]

# add branches 6 and 7
dx_direct_7 = target - elements[new_order_index[5],FastMultipole.POSITION]
u_direct_67 = u_direct_6 + ms[5] / sqrt(dx_direct_7' * dx_direct_7) * ONE_OVER_4PI

# reset target potential
mass_target_potential *= 0

# use summed multipole expansion from branches 6 and 7 (summed at 2)
FastMultipole.M2B!(mass_target_potential, mass_target, 2, tree)
u_fmm_67 = mass_target_potential[1]

# perform horizontal pass
m2l_list, direct_list = FastMultipole.build_interaction_lists(tree.branches, tree.branches, multipole_threshold, true, true, true)
FastMultipole.nearfield_singlethread!((elements,), tree.branches, (FastMultipole.DerivativesSwitch(),), (elements,), tree.branches, direct_list)
FastMultipole.horizontal_pass_singlethread!(tree.branches, tree.branches, m2l_list, expansion_order)

# consider the effect on branch 3 (mass 2)
elements.potential[i_POTENTIAL,new_order_index[2]] .*= 0 # reset potential at mass 2
# P2P is performed from branches 3 (mass 2), 4 (mass 3), and 5 (mass 1) to branch 3
# elements.direct!(elements.potential[i_POTENTIAL,new_order_index[2]], elements.bodies[i_POSITION,new_order_index[2]], elements.bodies[:,new_order_index[1]])
# elements.direct!(elements.potential[i_POTENTIAL,new_order_index[2]], elements.bodies[i_POSITION,new_order_index[2]], elements.bodies[:,new_order_index[2]])
# elements.direct!(elements.potential[i_POTENTIAL,new_order_index[2]], elements.bodies[i_POSITION,new_order_index[2]], elements.bodies[:,new_order_index[3]])
FastMultipole.P2P!((elements,), tree.branches[3], (FastMultipole.DerivativesSwitch(),), (elements,), tree.branches[3])
FastMultipole.P2P!((elements,), tree.branches[3], (FastMultipole.DerivativesSwitch(),), (elements,), tree.branches[4])
FastMultipole.P2P!((elements,), tree.branches[3], (FastMultipole.DerivativesSwitch(),), (elements,), tree.branches[5])
u_fmm_123 = elements.potential[i_POTENTIAL[1],new_order_index[2]]

dx_12 = elements[new_order_index[2],FastMultipole.POSITION] - elements[new_order_index[1],FastMultipole.POSITION]
u_direct_12 = elements.bodies[new_order_index[1]].strength[1] / sqrt(dx_12' * dx_12) * ONE_OVER_4PI
u_direct_22 = 0.0
dx_32 = elements[new_order_index[2],FastMultipole.POSITION] - elements[new_order_index[3],FastMultipole.POSITION]
u_direct_32 = elements.bodies[new_order_index[3]].strength[1] / sqrt(dx_32' * dx_32) * ONE_OVER_4PI

u_direct_123 = u_direct_12 + u_direct_22 + u_direct_32

# M2L is performed from branches 6, 7 to branch 3 (containing mass 2)
harmonics = zeros(Float64,2, (expansion_order+1)^2)
harmonics_theta = zeros(Float64,2, (expansion_order+1)^2)
harmonics_theta_2 = zeros(Float64,2, (expansion_order+1)^2)
workspace = zeros(3,4)
FastMultipole.L2B!((elements,), tree.branches[3], (FastMultipole.DerivativesSwitch(),), Val(expansion_order))
u_fmm_12345 = elements.potential[i_POTENTIAL[1],new_order_index[2]]

dx_42 = elements[new_order_index[4],FastMultipole.POSITION] - elements[new_order_index[2],FastMultipole.POSITION]
u_direct_42 = elements.bodies[new_order_index[4]].strength[1] / sqrt(dx_42' * dx_42) * ONE_OVER_4PI
dx_52 = elements[new_order_index[5],FastMultipole.POSITION] - elements[new_order_index[2],FastMultipole.POSITION]
u_direct_52 = elements.bodies[new_order_index[5]].strength[1] / sqrt(dx_52' * dx_52) * ONE_OVER_4PI

u_direct_12345 = u_direct_123 + u_direct_42 + u_direct_52

@test isapprox(u_direct_123, u_fmm_123; atol=1e-12)

@test isapprox(u_direct_12345, u_fmm_12345; atol=1e-12)

# reset potentials
elements.potential .*= 0

# run FastMultipole.(reset potentials with reset_tree flag)
FastMultipole.fmm!(tree, (elements,); multipole_threshold=multipole_threshold, reset_tree=true)
u_fmm = deepcopy(elements.potential[1,:])

elements.potential .= 0.0

FastMultipole.direct!((elements,))

u_direct = deepcopy(elements.potential[1,:])

for i in 1:FastMultipole.get_n_bodies(elements)
    @test isapprox(u_fmm[i], u_direct[i]; atol=1e-12)
end

end

#####
##### vector potential
#####

include("vortex.jl")

@testset "derivatives" begin
"""
dr_k/dx_idx_j
"""
function d2rdx2(r, theta, phi)
    derivatives = zeros(3,3,3)
    derivatives[:,:,1] .= [
        (1-cos(phi)^2 * sin(theta)^2)/r -sin(theta)^2*cos(phi)*sin(phi)/r -sin(theta)*cos(phi)*cos(theta)/r;
        (-sin(theta)^2*cos(phi)*sin(phi))/r (1-sin(theta)^2*sin(phi)^2)/r -sin(theta)*sin(phi)*cos(theta)/r;
        -sin(theta)*cos(phi)*cos(theta)/r -sin(theta)*sin(phi)*cos(theta)/r sin(theta)^2/r
    ]
    derivatives[:,:,2] .= [
        cos(theta)/sin(theta)*(1-cos(phi)^2*(1+2*sin(theta)^2))/r^2 -cos(theta)/sin(theta)*sin(phi)*cos(phi)*(1+2*sin(theta)^2)/r^2 cos(phi)*(1-2*cos(theta)^2)/r^2;
        -cos(theta)/sin(theta)*sin(phi)*cos(phi)*(1+2*sin(theta)^2)/r^2 cos(theta)/sin(theta)*(1-sin(phi)^2*(1+2*sin(theta)^2))/r^2 (2*sin(theta)^2-1)/r^2*sin(phi);
        cos(phi)*(1-2*cos(theta)^2)/r^2 (2*sin(theta)^2-1)/r^2*sin(phi) 2*sin(theta)*cos(theta)/r^2
    ]
    derivatives[:,:,3] .= [
        2*cos(phi)*sin(phi)/r^2/sin(theta)^2 (2*sin(phi)^2-1)/r^2/sin(theta)^2 0;
        (2*sin(phi)^2-1)/r^2/sin(theta)^2 -2*sin(phi)*cos(phi)/r^2/sin(theta)^2 0;
        0 0 0
    ]
    return derivatives
end

function cartesian_2_spherical(x,y,z; vec=false)
    r = sqrt(x^2+y^2+z^2)
    theta = acos(z/r)
    phi = atan(y,x)
    vec && return [r,theta,phi]
    return r, theta, phi
end

function d2rdx2_cart(x,y,z)
    r, theta, phi = cartesian_2_spherical(x,y,z)
    return d2rdx2(r,theta,phi)
end

"""
dr_j/dx_i
"""
function drdx(r,theta,phi)
    derivatives = [
        sin(theta)*cos(phi) cos(theta)*cos(phi)/r -sin(phi)/r/sin(theta);
        sin(theta)*sin(phi) cos(theta)*sin(phi)/r cos(phi)/r/sin(theta);
        cos(theta) -sin(theta)/r 0
    ]
    # derivatives = [
    #     sin(theta)*cos(phi) sin(theta)*sin(phi) cos(theta);
    #     cos(theta)*cos(phi)/r cos(theta)*sin(phi)/r -sin(theta)/r;
    #     -sin(phi)/r/sin(theta) cos(phi)/r/sin(theta) 0
    # ]
    return derivatives
end

function drdx_cart(x,y,z)
    r, theta, phi = cartesian_2_spherical(x,y,z)
    return drdx(r, theta, phi)
end

function d2rdx2_fd(x,y,z)
    derivatives = zeros(3,3,3)
    derivatives[1,:,:] .= (drdx_cart(x+1e-6,y,z) - drdx_cart(x,y,z))/1e-6
    derivatives[2,:,:] .= (drdx_cart(x,y+1e-6,z) - drdx_cart(x,y,z))/1e-6
    derivatives[3,:,:] .= (drdx_cart(x,y,z+1e-6) - drdx_cart(x,y,z))/1e-6
    return derivatives
end

x,y,z = rand(3)

function drdx_fd(x,y,z)
    derivatives = zeros(3,3)
    derivatives[1,:] .= (cartesian_2_spherical(x+1e-6,y,z;vec=true) - cartesian_2_spherical(x,y,z;vec=true))/1e-6
    derivatives[2,:] .= (cartesian_2_spherical(x,y+1e-6,z;vec=true) - cartesian_2_spherical(x,y,z;vec=true))/1e-6
    derivatives[3,:] .= (cartesian_2_spherical(x,y,z+1e-6;vec=true) - cartesian_2_spherical(x,y,z;vec=true))/1e-6
    return derivatives
end

fd_1 = drdx_fd(x,y,z)
anal_1 = drdx_cart(x,y,z)

for i in 1:length(fd_1)
    @test isapprox(fd_1[i], anal_1[i]; atol=1e-4)
end

fd = d2rdx2_fd(x,y,z)
anal = d2rdx2_cart(x,y,z)

for i in 1:length(fd)
    @test isapprox(fd[i], anal[i]; atol=1e-3)
end

expansion_order = 9
Random.seed!(123)
local_expansion = rand(2,4,55)
target = SVector{3}([0.0, 0.5, 0.0])
center = SVector{3}([0.01, 0.52, -0.03])

phi_fd(x) = FastMultipole.L2B(x, center, local_expansion, FastMultipole.DerivativesSwitch(), Val(expansion_order))[1]
psi_fd(x) = FastMultipole.L2B(x, center, local_expansion, FastMultipole.DerivativesSwitch(), Val(expansion_order))[2]

function vector_velocity_fd(x)
	j_psi = ForwardDiff.jacobian(psi_fd, x)
	v_psi = SVector{3}(j_psi[3,2] - j_psi[2,3],
		j_psi[1,3] - j_psi[3,1],
	 	j_psi[2,1] - j_psi[1,2])
	return v_psi
end

function scalar_velocity_fd(x)
	v_phi = -ForwardDiff.gradient(phi_fd, x)
	return v_phi
end

scalar_velocity_gradient_fd(x) = ForwardDiff.jacobian(scalar_velocity_fd, x)
vector_velocity_gradient_fd(x) = ForwardDiff.jacobian(vector_velocity_fd, x)

_, _, scalar_velocity_check_fmm, scalar_velocity_gradient_check_fmm = FastMultipole.L2B(target, center, local_expansion, FastMultipole.DerivativesSwitch(true, false, true, true), Val(expansion_order))
_, _, vector_velocity_check_fmm, vector_velocity_gradient_check_fmm = FastMultipole.L2B(target, center, local_expansion, FastMultipole.DerivativesSwitch(false, true, true, true), Val(expansion_order))

scalar_velocity_check_fd = scalar_velocity_fd(target)
vector_velocity_check_fd = vector_velocity_fd(target)
scalar_velocity_gradient_check_fd = scalar_velocity_gradient_fd(target)
vector_velocity_gradient_check_fd = vector_velocity_gradient_fd(target)

@test isapprox(scalar_velocity_check_fmm, scalar_velocity_check_fd; atol=1e-12)
@test isapprox(vector_velocity_check_fmm, vector_velocity_check_fd; atol=1e-12)
@test isapprox(scalar_velocity_gradient_check_fmm, scalar_velocity_gradient_check_fd; atol=1e-12)
@test isapprox(vector_velocity_gradient_check_fmm, vector_velocity_gradient_check_fd; atol=1e-12)

end

@testset "2D vortex particles" begin

# 2-D vortex ring
xs = [
    0 0
    0.5 -0.5
    0 0
]

Gammas = [
    0 0;
    0 0;
    1 -1.0
]

vortexparticles = VortexParticles(xs, Gammas)

# using direct method
FastMultipole.direct!((vortexparticles,))
update_velocity_stretching!(vortexparticles)

@test isapprox(vortexparticles.velocity_stretching[1,1], 1/4/pi; atol=1e-10)
@test isapprox(vortexparticles.velocity_stretching[1,2], 1/4/pi; atol=1e-10)
@test isapprox(vortexparticles.velocity_stretching[2,1], 0; atol=1e-10)
@test isapprox(vortexparticles.velocity_stretching[2,2], 0; atol=1e-10)
@test isapprox(vortexparticles.velocity_stretching[3,1], 0; atol=1e-10)
@test isapprox(vortexparticles.velocity_stretching[3,2], 0; atol=1e-10)

vorton_potential_check = zeros(4,2)
vorton_potential_check[i_POTENTIAL_VECTOR,:] = deepcopy(vortexparticles.potential[i_POTENTIAL_VECTOR,1:2])

vorton_velocity_check = deepcopy(vortexparticles.velocity_stretching[i_VELOCITY_vortex,:])

# reset vortons
vortexparticles.potential .*= 0
vortexparticles.velocity_stretching .*= 0

# manually build tree for testing
# Branch(n_branches, n_bodies, first_branch, first_body, center, radius, multipole_expansion, local_expansion, lock, child_lock)
expansion_order = 9
leaf_size = 1
x_branch_1 = FastMultipole.SVector{3}([0.0,0,0])
branch_1 = FastMultipole.MultiBranch(SVector{1}([1:2]), 2, 2:3, 0, x_branch_1, 1/8, FastMultipole.initialize_expansion(expansion_order), FastMultipole.initialize_expansion(expansion_order), FastMultipole.initialize_harmonics(expansion_order, Float64), FastMultipole.initialize_ML(expansion_order, Float64), ReentrantLock())
x_branch_2 = FastMultipole.SVector{3}(xs[:,1] .+ [0.01, 0.02, -0.03])
branch_2 = FastMultipole.MultiBranch(SVector{1}([1:1]), 0, 3:2, 1, x_branch_2, 1/8, FastMultipole.initialize_expansion(expansion_order), FastMultipole.initialize_expansion(expansion_order), FastMultipole.initialize_harmonics(expansion_order, Float64), FastMultipole.initialize_ML(expansion_order, Float64), ReentrantLock())
x_branch_3 = FastMultipole.SVector{3}(xs[:,2] .+ [0.02, -0.04, 0.01])
branch_3 = FastMultipole.MultiBranch(SVector{1}([2:2]), 0, 3:2, 1, x_branch_3, 1/8, FastMultipole.initialize_expansion(expansion_order), FastMultipole.initialize_expansion(expansion_order), FastMultipole.initialize_harmonics(expansion_order, Float64), FastMultipole.initialize_ML(expansion_order, Float64), ReentrantLock())

# using FMM
# tree = FastMultipole.Tree(branches, [expansion_order], leaf_size, B2M!, P2P!)
dummy_index = (zeros(Int64,FastMultipole.get_n_bodies(vortexparticles)),)
dummy_leaf_index = collect(1:3)
# dummy_cost_parameter = FastMultipole.dummy_direct_cost_estimate((vortexparticles,), leaf_size)
tree = FastMultipole.MultiTree([branch_1, branch_2, branch_3], [1:1,2:3], dummy_leaf_index, dummy_index, dummy_index, (deepcopy(vortexparticles),), Val(expansion_order), leaf_size)#, dummy_cost_parameter)
harmonics = zeros(Complex{Float64},2,(expansion_order+1)^2)
FastMultipole.B2M!(branch_2, (vortexparticles,), harmonics, Val(expansion_order))
FastMultipole.B2M!(branch_3, (vortexparticles,), harmonics, Val(expansion_order))
# @show tree.branches[2].multipole_expansion[2:4,:] # checks out

m2l_harmonics = zeros(eltype(tree.branches[1].multipole_expansion), 2, (expansion_order<<1 + 1)*(expansion_order<<1 + 1))
L = zeros(eltype(tree.branches[1].local_expansion), 2, 4)
FastMultipole.M2L!(tree.branches[2], tree.branches[3], m2l_harmonics, L, Val(expansion_order))
FastMultipole.M2L!(tree.branches[3], tree.branches[2], m2l_harmonics, L, Val(expansion_order))
# @show tree.branches[2].local_expansion
harmonics = zeros(Float64,2, (expansion_order+1)^2)
harmonics_theta = zeros(Float64,2, (expansion_order+1)^2)
harmonics_theta_2 = zeros(Float64,2, (expansion_order+1)^2)
workspace = zeros(3,4)
FastMultipole.L2B!((vortexparticles,), tree.branches[2], (FastMultipole.DerivativesSwitch(),), Val(expansion_order))
FastMultipole.L2B!((vortexparticles,), tree.branches[3], (FastMultipole.DerivativesSwitch(),), Val(expansion_order))
update_velocity_stretching!(vortexparticles)

for i in 1:2
    for ind in 1:4
        @test isapprox(vortexparticles.potential[ind,i], vorton_potential_check[ind,i]; atol=1e-12)
    end
    for dim in 1:3
        @test isapprox(vortexparticles.velocity_stretching[dim,i], vorton_velocity_check[dim,i]; atol=1e-12)
    end
end

end

function psi(target_x, source_x, source_gamma)
    dx = target_x - source_x
    dx_norm = sqrt(dx' * dx)
    return source_gamma ./ dx_norm * ONE_OVER_4PI
end

function dpsidx(target_x, source_x, source_gamma)
    dx = target_x - source_x
    dx_norm = sqrt(dx' * dx)
    x, y, z = dx
    jacobian = [
        -x*source_gamma[1] -x*source_gamma[2] -x*source_gamma[3];
        -y*source_gamma[1] -y*source_gamma[2] -y*source_gamma[3];
        -z*source_gamma[1] -z*source_gamma[2] -z*source_gamma[3];
    ] ./ dx_norm^3 * ONE_OVER_4PI
    return jacobian
end

function d2psidx2(target_x, source_x, source_gamma)
    dx = target_x - source_x
    dx_norm = sqrt(dx' * dx)
    x, y, z = dx
    hessian = zeros(3,3,3)
    d2dr2 = [
        2x^2-y^2-z^2 3x*y 3x*z;
        3x*y 2y^2-x^2-z^2 3y*z;
        3x*z 3y*z 2z^2-x^2-y^2
    ] / dx_norm^5
    hessian[:,:,1] = d2dr2 * source_gamma[1] * ONE_OVER_4PI
    hessian[:,:,2] = d2dr2 * source_gamma[2] * ONE_OVER_4PI
    hessian[:,:,3] = d2dr2 * source_gamma[3] * ONE_OVER_4PI
    return hessian
end

function u(target_x, source_x, source_gamma)
    dx = target_x  - source_x
    dx_norm = sqrt(dx' * dx)
    return 1/4/pi/dx_norm^3 * [
        -dx[2]*source_gamma[3] + dx[3]*source_gamma[2],
        -dx[3]*source_gamma[1] + dx[1]*source_gamma[3],
        -dx[1]*source_gamma[2] + dx[2]*source_gamma[1]
    ]
end

function duidxj_fd_fun(target_x, source_x, source_gamma; h=1e-8)
    duidx = (u(target_x+[h,0,0], source_x, source_gamma) - u(target_x,source_x,source_gamma))/h
    duidy = (u(target_x+[0,h,0], source_x, source_gamma) - u(target_x,source_x,source_gamma))/h
    duidz = (u(target_x+[0,0,h], source_x, source_gamma) - u(target_x,source_x,source_gamma))/h
    duidxj_res = hcat(duidx, duidy, duidz) .* 4 * pi
    return duidxj_res
end

function stretching(target_x, source_x, target_gamma, source_gamma)
    dx = target_x - source_x
    x, y, z = dx
    xy = x*y
    yz = y*z
    xz = x*z
    gx, gy, gz = source_gamma
    dx_norm = sqrt(dx' * dx)
    duidxj = [
        (3xy*gz-3xz*gy) ((2y^2-x^2-z^2)*gz-3yz*gy) (3yz*gz-(2z^2-x^2-y^2)*gy);
        (3xz*gx-(2x^2-y^2-z^2)*gz) (3yz*gx-3xy*gz) ((2z^2-x^2-y^2)*gx-3xz*gz);
        ((2x^2-y^2-z^2)*gy-3xy*gx) (3xy*gy-(2y^2-x^2-z^2)*gx) (3xz*gy-3yz*gx)
    ]/dx_norm^5
    stretch = 1/4/pi*duidxj*target_gamma
    return stretch
end

@testset "two 3D vortex particles" begin

bodies = [
    0.4 0.1
    0.1 -0.5
    -0.3 0.2
    1/8 1/8
    0.3 -0.4
    -0.1 -0.2
    0.08 0.5
]

vortex_particles = VortexParticles(bodies)

#####
##### obtain psi, u, and stretching analytically
#####

psis = zeros(3,2)
psis[:,1] = psi(bodies[1:3,1], bodies[1:3,2], bodies[5:7,2])
psis[:,2] = psi(bodies[1:3,2], bodies[1:3,1], bodies[5:7,1])
hessians = zeros(3,3,3,2)
hessians[:,:,:,1] = d2psidx2(bodies[1:3,1], bodies[1:3,2], bodies[5:7,2])
hessians[:,:,:,2] = d2psidx2(bodies[1:3,2], bodies[1:3,1], bodies[5:7,1])
us = zeros(3,2)
us[:,1] = u(bodies[1:3,1], bodies[1:3,2], bodies[5:7,2])
us[:,2] = u(bodies[1:3,2], bodies[1:3,1], bodies[5:7,1])
ss = zeros(3,2)
ss[:,1] = stretching(bodies[1:3,1], bodies[1:3,2], bodies[5:7,1], bodies[5:7,2])
ss[:,2] = stretching(bodies[1:3,2], bodies[1:3,1], bodies[5:7,2], bodies[5:7,1])

#####
##### use direct method
#####
FastMultipole.direct!((vortex_particles,))
update_velocity_stretching!(vortex_particles)

psis_direct = deepcopy(vortex_particles.potential[2:4,:])
for i in 1:length(psis_direct)
    @test isapprox(psis_direct[i], psis[i]; atol=1e-10)
end
# hessians_direct = deepcopy(reshape(vortex_particles.potential[i_POTENTIAL_HESSIAN[10:end],:],3,3,3,2))
# for i in 1:length(hessians)
#     @test isapprox(hessians_direct[i], hessians[i]; atol=1e-10)
# end
us_direct = deepcopy(vortex_particles.velocity_stretching[1:3,:])
for i in 1:length(us)
    @test isapprox(us_direct, us;atol=1e-10)
end
ss_direct = deepcopy(vortex_particles.velocity_stretching[4:6,:])
for i in 1:length(ss)
    @test isapprox(ss_direct[i], ss[i];atol=1e-10)
end

#####
##### use fmm
#####
# reset potential
vortex_particles.potential .*= 0
vortex_particles.velocity_stretching .*= 0

# branch = Branch(n_branches, n_bodies, i_child, i_start, center, radius, multipole_expansion, local_expansion)
expansion_order = 9
leaf_size = 1
x_branch_1 = SVector{3}((bodies[1:3,1] + bodies[1:3,2])/2)
branch_1 = FastMultipole.MultiBranch(SVector{1}([1:2]), 2, 2:3, 0, x_branch_1, 1/8, FastMultipole.initialize_expansion(expansion_order), FastMultipole.initialize_expansion(expansion_order), FastMultipole.initialize_harmonics(expansion_order, Float64), FastMultipole.initialize_ML(expansion_order, Float64), ReentrantLock())
x_branch_2 = SVector{3}(bodies[1:3,1] .+ [0.01, 0.02, -0.03])
branch_2 = FastMultipole.MultiBranch(SVector{1}([1:1]), 0, 3:2, 1, x_branch_2, 1/8, FastMultipole.initialize_expansion(expansion_order), FastMultipole.initialize_expansion(expansion_order), FastMultipole.initialize_harmonics(expansion_order, Float64), FastMultipole.initialize_ML(expansion_order, Float64), ReentrantLock())
x_branch_3 = SVector{3}(bodies[1:3,2] .+ [0.02, -0.04, 0.01])
branch_3 = FastMultipole.MultiBranch(SVector{1}([2:2]), 0, 3:2, 1, x_branch_2, 1/8, FastMultipole.initialize_expansion(expansion_order), FastMultipole.initialize_expansion(expansion_order), FastMultipole.initialize_harmonics(expansion_order, Float64), FastMultipole.initialize_ML(expansion_order, Float64), ReentrantLock())

dummy_index = (zeros(Int,length(vortex_particles.bodies)),)
dummy_leaf_index = collect(1:3)
# dummy_cost_parameter = FastMultipole.dummy_direct_cost_estimate((vortex_particles,), leaf_size)
tree = FastMultipole.MultiTree([branch_1, branch_2, branch_3], [1:1,2:3], dummy_leaf_index, dummy_index, dummy_index, (deepcopy(vortex_particles),), Val(expansion_order), leaf_size)#, dummy_cost_parameter)
# FastMultipole.B2M!(tree, vortex_particles, 2)
# FastMultipole.B2M!(tree, vortex_particles, 3)

# FastMultipole.M2L!(tree, 2, 3)
# FastMultipole.M2L!(tree, 3, 2)
# FastMultipole.L2B!(tree, vortex_particles, 2)
# FastMultipole.L2B!(tree, vortex_particles, 3)
FastMultipole.fmm!(tree, (vortex_particles,); multipole_threshold=0.5, unsort_bodies=false)
update_velocity_stretching!(vortex_particles)

psis_fmm = deepcopy(vortex_particles.potential[2:4,:])
for i in 1:length(psis_fmm)
    @test isapprox(psis_fmm[i], psis[i]; atol=1e-10)
end
# hessians_fmm.= deepcopy(reshape(vortex_particles.potential[i_POTENTIAL_HESSIAN[10:end],:],3,3,3,2))
# for i in 1:length(hessians)
#     @test isapprox(hessians_fmm.i], hessians[i]; atol=1e-8)
# end
us_fmm = deepcopy(vortex_particles.velocity_stretching[1:3,:])
for i in 1:length(us)
    @test isapprox(us_fmm, us;atol=1e-8)
end
ss_fmm = deepcopy(vortex_particles.velocity_stretching[4:6,:])
for i in 1:length(ss)
    @test isapprox(ss_fmm[i], ss[i];atol=1e-8)
end

end



@testset "three 3D vortex particles" begin

bodies = [
    0.4 0.1 -0.1
    0.1 -0.5 0.25
    -0.3 0.2 0.1
    1/8 1/8 1/8
    0.3 -0.4 0.2
    -0.1 -0.2 0.5
    0.08 0.5 1.1
]

vortex_particles = VortexParticles(bodies)

psis = zeros(3,3)
psis[:,1] = psi(bodies[1:3,1], bodies[1:3,2], bodies[5:7,2]) + psi(bodies[1:3,1], bodies[1:3,3], bodies[5:7,3])
psis[:,2] = psi(bodies[1:3,2], bodies[1:3,1], bodies[5:7,1]) + psi(bodies[1:3,2], bodies[1:3,3], bodies[5:7,3])
psis[:,3] = psi(bodies[1:3,3], bodies[1:3,1], bodies[5:7,1]) + psi(bodies[1:3,3], bodies[1:3,2], bodies[5:7,2])
jacobians = zeros(3,3,3)
jacobians[:,:,1] = dpsidx(bodies[1:3,1], bodies[1:3,2], bodies[5:7,2]) + dpsidx(bodies[1:3,1], bodies[1:3,3], bodies[5:7,3])
jacobians[:,:,2] = dpsidx(bodies[1:3,2], bodies[1:3,1], bodies[5:7,1]) + dpsidx(bodies[1:3,2], bodies[1:3,3], bodies[5:7,3])
jacobians[:,:,3] = dpsidx(bodies[1:3,3], bodies[1:3,1], bodies[5:7,1]) + dpsidx(bodies[1:3,3], bodies[1:3,2], bodies[5:7,2])
hessians = zeros(3,3,3,3)
hessians[:,:,:,1] = d2psidx2(bodies[1:3,1], bodies[1:3,2], bodies[5:7,2]) + d2psidx2(bodies[1:3,1], bodies[1:3,3], bodies[5:7,3])
hessians[:,:,:,2] = d2psidx2(bodies[1:3,2], bodies[1:3,1], bodies[5:7,1]) + d2psidx2(bodies[1:3,2], bodies[1:3,3], bodies[5:7,3])
hessians[:,:,:,3] = d2psidx2(bodies[1:3,3], bodies[1:3,1], bodies[5:7,1]) + d2psidx2(bodies[1:3,3], bodies[1:3,2], bodies[5:7,2])
us = zeros(3,3)
us[:,1] = u(bodies[1:3,1], bodies[1:3,2], bodies[5:7,2]) + u(bodies[1:3,1], bodies[1:3,3], bodies[5:7,3])
us[:,2] = u(bodies[1:3,2], bodies[1:3,1], bodies[5:7,1]) + u(bodies[1:3,2], bodies[1:3,3], bodies[5:7,3])
us[:,3] = u(bodies[1:3,3], bodies[1:3,1], bodies[5:7,1]) + u(bodies[1:3,3], bodies[1:3,2], bodies[5:7,2])
ss = zeros(3,3)
ss[:,1] = stretching(bodies[1:3,1], bodies[1:3,2], bodies[5:7,1], bodies[5:7,2]) + stretching(bodies[1:3,1], bodies[1:3,3], bodies[5:7,1], bodies[5:7,3])
ss[:,2] = stretching(bodies[1:3,2], bodies[1:3,1], bodies[5:7,2], bodies[5:7,1]) + stretching(bodies[1:3,2], bodies[1:3,3], bodies[5:7,2], bodies[5:7,3])
ss[:,3] = stretching(bodies[1:3,3], bodies[1:3,1], bodies[5:7,3], bodies[5:7,1]) + stretching(bodies[1:3,3], bodies[1:3,2], bodies[5:7,3], bodies[5:7,2])

#####
##### use direct method
#####
FastMultipole.direct!((vortex_particles,))
update_velocity_stretching!(vortex_particles)

psis_direct = deepcopy(vortex_particles.potential[2:4,:])
for i in 1:length(psis_direct)
    @test isapprox(psis_direct[i], psis[i]; atol=1e-10)
end
# jacobians_direct = deepcopy(reshape(vortex_particles.potential[i_POTENTIAL_JACOBIAN[4:end],:],3,3,3))
# for i in 1:length(jacobians_direct)
#     @test isapprox(jacobians_direct[i], jacobians[i]; atol=1e-10)
# end
# hessians_direct = deepcopy(reshape(vortex_particles.potential[i_POTENTIAL_HESSIAN[10:end],:],3,3,3,3))
# for i in 1:length(hessians)
#     @test isapprox(hessians_direct[i], hessians[i]; atol=1e-10)
# end
us_direct = deepcopy(vortex_particles.velocity_stretching[1:3,:])
for i in 1:length(us)
    @test isapprox(us_direct[i], us[i];atol=1e-10)
end
ss_direct = deepcopy(vortex_particles.velocity_stretching[4:6,:])
for i in 1:length(ss)
    @test isapprox(ss_direct[i], ss[i];atol=1e-10)
end

#####
##### use fmm
#####
# reset potential
vortex_particles.potential .*= 0
vortex_particles.velocity_stretching .*= 0

expansion_order = 32
leaf_size = 1

tree = FastMultipole.Tree((vortex_particles,); expansion_order, leaf_size, shrink_recenter=false)

FastMultipole.fmm!(tree, (vortex_particles,); multipole_threshold=0.5)

update_velocity_stretching!(vortex_particles)

psis_fmm = deepcopy(vortex_particles.potential[2:4,:])
for i in 1:length(psis_fmm)
    @test isapprox(psis_fmm[i], psis[i]; rtol=1e-12)
end
# hessians_fmm.= deepcopy(reshape(vortex_particles.potential[i_POTENTIAL_HESSIAN[10:end],:],3,3,3,3))
# for i in 1:length(hessians)
#     @test isapprox(hessians_fmm.i], hessians[i]; rtol=1e-12)
# end
us_fmm = deepcopy(vortex_particles.velocity_stretching[1:3,:])
for i in 1:length(us)
    @test isapprox(us_fmm, us;rtol=1e-12)
end
ss_fmm = deepcopy(vortex_particles.velocity_stretching[4:6,:])
for i in 1:length(ss)
    @test isapprox(ss_fmm[i], ss[i];rtol=1e-12)
end

end

function generate_gravitational(seed, n_bodies; radius_factor=0.1)
    Random.seed!(123)
    bodies = rand(8,n_bodies)
    bodies[4,:] ./= (n_bodies^(1/3)*2)
    bodies[4,:] .*= radius_factor
    system = Gravitational(bodies)
end

@testset "single/dual tree, single/multi branch" begin

expansion_order, leaf_size, multipole_threshold = 14, 100, 0.31
n_bodies = 5000
shrink_recenter = true
seed = 123
validation_system = generate_gravitational(seed, n_bodies; radius_factor=0.1)
FastMultipole.direct!(validation_system)
validation_potential = validation_system.potential[1,:]

system3 = generate_gravitational(seed, n_bodies; radius_factor=0.1)
FastMultipole.fmm!(system3; expansion_order=expansion_order, leaf_size=leaf_size, multipole_threshold=multipole_threshold, nearfield=true, farfield=true, unsort_bodies=true)
potential3 = system3.potential[1,:]
@test isapprox(maximum(abs.(potential3 - validation_potential)), 0.0; atol=1e-9)

system4 = generate_gravitational(seed, n_bodies; radius_factor=0.1)
FastMultipole.fmm!((system4,); expansion_order=expansion_order, leaf_size=leaf_size, multipole_threshold=multipole_threshold, nearfield=true, farfield=true, unsort_bodies=true)
potential4 = system4.potential[1,:]
@test isapprox(maximum(abs.(potential4 - validation_potential)), 0.0; atol=1e-9)

system5 = generate_gravitational(seed, n_bodies; radius_factor=0.1)
FastMultipole.fmm!(system5, system5; expansion_order=expansion_order, leaf_size_source=leaf_size, leaf_size_target=leaf_size, multipole_threshold=multipole_threshold, nearfield=true, farfield=true, unsort_source_bodies=true, unsort_target_bodies=true)
potential5 = system5.potential[1,:]
@test isapprox(maximum(abs.(potential5 - validation_potential)), 0.0; atol=1e-9)

system6 = generate_gravitational(seed, n_bodies; radius_factor=0.1)
FastMultipole.fmm!((system6,), system6; expansion_order=expansion_order, leaf_size_source=leaf_size, leaf_size_target=leaf_size, multipole_threshold=multipole_threshold, nearfield=true, farfield=true, unsort_source_bodies=true, unsort_target_bodies=true)
potential6 = system6.potential[1,:]
@test isapprox(maximum(abs.(potential6 - validation_potential)), 0.0; atol=1e-9)

system7 = generate_gravitational(seed, n_bodies; radius_factor=0.1)
FastMultipole.fmm!((system7,), (system7,); expansion_order=expansion_order, leaf_size_source=leaf_size, leaf_size_target=leaf_size, multipole_threshold=multipole_threshold, nearfield=true, farfield=true, unsort_source_bodies=true, unsort_target_bodies=true)
potential7 = system7.potential[1,:]
@test isapprox(maximum(abs.(potential7 - validation_potential)), 0.0; atol=1e-9)

end

@testset "sortwrapper" begin

expansion_order, leaf_size, multipole_threshold = 14, 100, 0.31
n_bodies = 5000
shrink_recenter = true
seed = 123
validation_system = generate_gravitational(seed, n_bodies; radius_factor=0.1)
FastMultipole.direct!(validation_system)
validation_potential = validation_system.potential[1,:]

system8 = FastMultipole.SortWrapper(generate_gravitational(seed, n_bodies; radius_factor=0.1))
FastMultipole.fmm!(system8; expansion_order=expansion_order, leaf_size=leaf_size, multipole_threshold=multipole_threshold, nearfield=true, farfield=true, unsort_bodies=true)
potential8 = system8.system.potential[1,:]
@test isapprox(maximum(abs.(potential8 - validation_potential)), 0.0; atol=1e-9)

system9 = FastMultipole.SortWrapper(generate_gravitational(seed, n_bodies; radius_factor=0.1))
FastMultipole.fmm!((system9,); expansion_order=expansion_order, leaf_size=leaf_size, multipole_threshold=multipole_threshold, nearfield=true, farfield=true, unsort_bodies=true)
potential9 = system9.system.potential[1,:]
@test isapprox(maximum(abs.(potential9 - validation_potential)), 0.0; atol=1e-9)

system10 = FastMultipole.SortWrapper(generate_gravitational(seed, n_bodies; radius_factor=0.1))
FastMultipole.fmm!(system10, system10; expansion_order=expansion_order, leaf_size_source=leaf_size, leaf_size_target=leaf_size, multipole_threshold=multipole_threshold, nearfield=true, farfield=true, unsort_source_bodies=true, unsort_target_bodies=true)
potential10 = system10.system.potential[1,:]
@test isapprox(maximum(abs.(potential10 - validation_potential)), 0.0; atol=1e-9)

system11 = FastMultipole.SortWrapper(generate_gravitational(seed, n_bodies; radius_factor=0.1))
FastMultipole.fmm!((system11,), system11; expansion_order=expansion_order, leaf_size_source=leaf_size, leaf_size_target=leaf_size, multipole_threshold=multipole_threshold, nearfield=true, farfield=true, unsort_source_bodies=true, unsort_target_bodies=true)
potential11 = system11.system.potential[1,:]
@test isapprox(maximum(abs.(potential11 - validation_potential)), 0.0; atol=1e-9)

system12 = FastMultipole.SortWrapper(generate_gravitational(seed, n_bodies; radius_factor=0.1))
FastMultipole.fmm!((system12,), (system12,); expansion_order=expansion_order, leaf_size_source=leaf_size, leaf_size_target=leaf_size, multipole_threshold=multipole_threshold, nearfield=true, farfield=true, unsort_source_bodies=true, unsort_target_bodies=true)
potential12 = system12.system.potential[1,:]
@test isapprox(maximum(abs.(potential12 - validation_potential)), 0.0; atol=1e-9)

end

@testset "b2m source and dipole panels" begin

# source panel expansion
strength = 0.1428909901797533 / 4 / pi
R0_global = [-0.11530793537122011, 0.1997192025788218, -0.9730448705798238]
R0 = [-0.48862125790260025, -0.1735941199525583, -1.8441092898197107]
Ru = [0.04707023255275322, -0.10579806212499904, -0.020193487162119217]
Rv = [-0.039004202054645, -0.02833821156361363, 0.0]
normal = cross(Rv,Ru)
normal /= norm(normal)
qnm_prev = [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]
jnm_prev = [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]
inm_prev = [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]
this_expansion_order = 14
n_coefficients = ((this_expansion_order+1) * (this_expansion_order+2)) >> 1
multipole_expansion = zeros(2,4,n_coefficients)
multipole_expansion_test = [3.15347369905358e-5 0.0 ; 5.836576687179473e-5 0.0 ; -7.661878069585717e-6 3.442114476513879e-6 ; 5.176954094092041e-5 0.0 ; -1.4179736287285782e-5 6.372425519286581e-6 ; 7.42180196746523e-7 -8.340091929224401e-7 ; 2.9171240607705508e-5 0.0 ; -1.2848310937891748e-5 5.7763762881869375e-6 ; 1.3731914469928547e-6 -1.5438750448677817e-6 ; -2.981931723472997e-8 9.412368214502325e-8 ; 1.1616879034380539e-5 0.0 ; -7.589495677738776e-6 3.4136975727862394e-6 ; 1.2526909809559565e-6 -1.409195559129566e-6 ; -5.512043533657875e-8 1.742126360703301e-7 ; -7.340115650494341e-10 -6.495176725450135e-9 ; 3.411340562667307e-6 0.0 ; -3.2811221779399064e-6 1.4766336464937334e-6 ; 7.508033462390918e-7 -8.45140417586639e-7 ; -5.0408680959505405e-8 1.5954795989679934e-7 ; -1.3643415621015169e-9 -1.2019285292447762e-8 ; 1.7490047653047897e-10 2.9802226104568586e-10 ; 7.320829843849779e-7 0.0 ; -1.1040063516878706e-6 4.97169838052928e-7 ; 3.3233488404786064e-7 -3.7435741595385594e-7 ; -3.0399512277559916e-8 9.636699936286549e-8 ; -1.257805313390412e-9 -1.1028077206194902e-8 ; 3.239977186634296e-10 5.512947043790685e-10 ; -1.2339557354154617e-11 -8.868558661043951e-12 ; 1.0121462827844036e-7 0.0 ; -2.9977386626355873e-7 1.3510437997734334e-7 ; 1.1575725282596027e-7 -1.3049796958370626e-7 ; -1.3594389585018197e-8 4.3167709924881135e-8 ; -7.667019220226146e-10 -6.6880743655575325e-9 ; 2.9803138008504807e-10 5.0634784455077e-10 ; -2.284557373740596e-11 -1.6393720070527073e-11 ; 5.595126561276888e-13 1.1787154792703014e-13 ; 1.4751810727352213e-9 0.0 ; -6.708468877007811e-8 3.026408896811465e-8 ; 3.300201956443528e-8 -3.723787837566326e-8 ; -4.805896384168386e-9 1.5289010248988546e-8 ; -3.4754378686830237e-10 -3.0152402096102556e-9 ; 1.81479446950944e-10 3.07835507484142e-10 ; -2.1022722355702963e-11 -1.506079745721181e-11 ; 1.035467726604547e-12 2.1721150343128253e-13 ; -1.8395652399672212e-14 3.852107985386056e-15 ; -4.0609549751469545e-9 0.0 ; -1.2482117506512014e-8 5.638433267215509e-9 ; 7.905670401846284e-9 -8.929464291174913e-9 ; -1.3981037159203388e-9 4.456832611637532e-9 ; -1.2492898592786975e-10 -1.0775898987149912e-9 ; 8.228418174866392e-11 1.3933781542187157e-10 ; -1.2818959531322538e-11 -9.167593119553919e-12 ; 9.531511322909608e-13 1.9904649104822634e-13 ; -3.402893923417914e-14 7.160468249308245e-15 ; 4.465711619580417e-16 -3.1804698481402137e-16 ; -1.4824011156256148e-9 0.0 ; -1.919833813101314e-9 8.688650303478487e-10 ; 1.6200773449979007e-9 -1.8320295085795508e-9 ; -3.4396303878189783e-10 1.0989242195318225e-9 ; -3.708192715671552e-11 -3.1787656903648365e-10 ; 2.962574958499155e-11 5.007695694535446e-11 ; -5.8262901006951636e-12 -4.159070625703628e-12 ; 5.818200312063137e-13 1.2092681969926203e-13 ; -3.1327310064348024e-14 6.625396053471999e-15 ; 8.255548969437712e-16 -5.893192421106608e-16 ; -7.403943579266824e-18 1.2317599122655843e-17 ; -3.515934968301032e-10 0.0 ; -2.3594925434518903e-10 1.07127957503847e-10 ; 2.8743136385190034e-10 -3.2548499126874616e-10 ; -7.297407969269781e-11 2.3371958216102407e-10 ; -9.344496321234191e-12 -7.957394809536844e-11 ; 8.820872426341902e-12 1.4881530759394605e-11 ; -2.10510183229474e-12 -1.4998051808823597e-12 ; 2.6492865541168975e-13 5.4789022475871284e-14 ; -1.9135933669548904e-14 4.0683844309145415e-15 ; 7.598586435383344e-16 -5.43729257580773e-16 ; -1.366948657228694e-17 2.2803585253454954e-17 ; 4.513564991747822e-20 -3.389937360059917e-19 ; -6.573572046881621e-11 0.0 ; -2.043863801299287e-11 9.352606219154242e-12 ; 4.445562910300353e-11 -5.042557742977247e-11 ; -1.3543695917757064e-11 4.3496602156693426e-11 ; -2.0397436768532417e-12 -1.7246583947075334e-11 ; 2.2333444964799568e-12 3.76016740753104e-12 ; -6.297270638036358e-13 -4.477390964671143e-13 ; 9.597644376139124e-14 1.974437557076672e-14 ; -8.724659175859616e-15 1.8650669605364297e-15 ; 4.642657322793347e-16 -3.3304590676846863e-16 ; -1.2569421365745768e-17 2.1028191987394226e-17 ; 8.266548146039412e-20 -6.271851481459042e-19 ; 2.019255448918255e-21 7.17077959778311e-21 ; -1.0278419931152499e-11 0.0 ; -4.203242166211679e-13 2.1066117376944255e-13 ; 6.002693212564717e-12 -6.823256754648399e-12 ; -2.222360876797338e-12 7.159343863048163e-12 ; -3.915530871566449e-13 -3.2854874061920314e-12 ; 4.906975535591016e-13 8.243739866537346e-13 ; -1.603914412140437e-13 -1.1379290070808765e-13 ; 2.881187837425885e-14 5.8944206681181154e-15 ; -3.166737553935407e-15 6.808073260383793e-16 ; 2.11823342218566e-16 -1.5235001948604837e-16 ; -7.674869280589808e-18 1.2877827358552214e-17 ; 7.540003064291612e-20 -5.781607461451349e-19 ; 3.753862216640856e-21 1.3258129373615239e-20 ; -9.398251768197244e-23 -1.164083006095856e-22 ; -1.3762401781072233e-12 0.0 ; 2.8087332169094867e-13 -1.2390042397442847e-13 ; 7.032457492620986e-13 -8.016843386759864e-13 ; -3.248637046058368e-13 1.0502423099069615e-12 ; -6.687546341485134e-14 -5.565382595152538e-13 ; 9.500740818317414e-14 1.5924469307571103e-13 ; -3.5499221407458265e-14 -2.5128180166288814e-14 ; 7.370959401207499e-15 1.4991782388637415e-15 ; -9.53078284611239e-16 2.0611156438331795e-16 ; 7.697575427474568e-17 -5.551354761137937e-17 ; -3.500656875141304e-18 5.891922242973979e-18 ; 4.5664032899842155e-20 -3.540569052057126e-19 ; 3.4784857373427885e-21 1.2216677736179948e-20 ; -1.7414905033630244e-22 -2.1502849742163935e-22 ; 2.460117981125405e-24 1.3094830233058996e-24]
multipole_expansion_test2 = zeros(2,4,n_coefficients)

for i in 1:n_coefficients
    multipole_expansion_test2[1,1,i] = multipole_expansion_test[i,1]
    multipole_expansion_test2[2,1,i] = multipole_expansion_test[i,2]
end

expansion_order = Val{14}()
coefficients = view(multipole_expansion, :, 1, :)
FastMultipole._B2M!_panel(coefficients, qnm_prev, jnm_prev, inm_prev, R0, Ru, Rv, strength, normal, expansion_order, FastMultipole.Panel{3,FastMultipole.ConstantSource{1}})
for i in eachindex(multipole_expansion)
    @test isapprox(multipole_expansion[i], multipole_expansion_test2[i]; atol=1e-12)
end

# single tri panel
displ = rand(SVector{3,Float64})
vertices = [
    SVector{3}(0.0,0.0,-0.2) + displ,
    SVector{3}(2,1,0.3) + displ,
    SVector{3}(-1,2,0.1) + displ
]
target = [5.5, -4.2, 5.3] + displ
centroid = (vertices[1] + vertices[2] + vertices[3])/3
normal = cross(vertices[2]-vertices[1],vertices[3]-vertices[1])
area = norm(normal)/2
normal /= area*2
strength = 1.0 / 4 / pi
P = 20
expansion_order = Val(P)
branch = FastMultipole.SingleBranch(
    1:1,
    0,
    1:0,
    0,
    SVector{3}(0.0,0.0,0.0) + displ,
    0.4,
    FastMultipole.initialize_expansion(P),
    FastMultipole.initialize_expansion(P),
    zeros(2, (P+1)*(P+1)),
    zeros(2,4),
    ReentrantLock(),
)

# build expansion
qnm_prev = zeros(2, P+2)
jnm_prev = zeros(2, P+2)
inm_prev = zeros(2, P+2)
R0_global = vertices[1]
R0 = R0_global - branch.center
Ru = vertices[2] - R0_global
Rv = vertices[3] - R0_global
panel_type = FastMultipole.Panel{3, FastMultipole.ConstantNormalDipole{1}}
coefficients = view(branch.multipole_expansion, :, 1, :)
FastMultipole._B2M!_panel(coefficients, qnm_prev, jnm_prev, inm_prev, R0, Ru, Rv, strength, normal, expansion_order, panel_type)

# evaluate
target_potential_dipole = zeros(4)
irregular_harmonics_dipole = FastMultipole.initialize_harmonics(P,Float64)
multipole_expansion_dipole = branch.multipole_expansion
FastMultipole.M2B!(target_potential_dipole, target, displ, irregular_harmonics_dipole, multipole_expansion_dipole, P)

this_potential = 0.0015577438029241901

@test isapprox(target_potential_dipole[1], this_potential; atol=1e-12)

# reset and repeat for source panel
branch.multipole_expansion .= 0.0

panel_type = FastMultipole.Panel{3, FastMultipole.ConstantSource{1}}
coefficients = view(branch.multipole_expansion, :, 1, :)
FastMultipole._B2M!_panel(coefficients, qnm_prev, jnm_prev, inm_prev, R0, Ru, Rv, strength, normal, expansion_order, panel_type)

target_potential = zeros(4)
irregular_harmonics = FastMultipole.initialize_harmonics(P,Float64)
multipole_expansion = branch.multipole_expansion
FastMultipole.M2B!(target_potential, target, displ, irregular_harmonics, multipole_expansion, P)

@test isapprox(target_potential[1], 0.022849841377239076; atol=1e-12)

end

@testset "probes" begin

n_bodies = 101
bodies = rand(8,n_bodies)
bodies[5:8,2:n_bodies] .= 0.0
probes = FastMultipole.ProbeSystem(bodies[1:3,2:end]; scalar_potential=true, vector_potential=true, velocity=true, velocity_gradient=true)
mass = Gravitational(bodies)
FastMultipole.fmm!(probes, mass; expansion_order=5, leaf_size_source=50, leaf_size_target=50, multipole_threshold=0.4)
FastMultipole.fmm!(mass; expansion_order=5, leaf_size=50, multipole_threshold=0.4)

for i_mass in 2:n_bodies
    @test isapprox(mass.potential[1,i_mass], probes.scalar_potential[i_mass-1]; atol=1e-12)
end

end

