using Test
import Statistics
S = Statistics

import FastMultipole
fmm = FastMultipole

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

fmm.M2M!(branch, child, harmonics, M, expansion_order::Int) = fmm.M2M!(branch, child, harmonics, M, Val(expansion_order))
fmm.L2B_loop!(vector_potential, potential_jacobian, potential_hessian, body_position, expansion_center, local_expansion, harmonics, harmonics_theta, harmonics_theta_2, expansion_order::Int, workspace) = 
    fmm.L2B_loop!(vector_potential, potential_jacobian, potential_hessian, body_position, expansion_center, local_expansion, harmonics, harmonics_theta, harmonics_theta_2, Val(expansion_order), workspace)
fmm.M2L!(target_branch, source_branch, harmonics, L, expansion_order::Int) = fmm.M2L!(target_branch, source_branch, harmonics, L, Val(expansion_order))
fmm.upward_pass_singlethread!(branches, systems, expansion_order::Int) = fmm.upward_pass_singlethread!(branches, systems, Val(expansion_order))
fmm.M2L!(target_branch, source_branch, expansion_order::Int) = fmm.M2L!(target_branch, source_branch, Val(expansion_order))

@testset "complex" begin
    z1 = rand(Complex{Float64})
    z2 = rand(Complex{Float64})
    z3 = rand(Complex{Float64})
    
    # addition
    @test real(z1+z2) ≈ fmm.complex_add(real(z1), imag(z1), real(z2), imag(z2))[1]
    @test imag(z1+z2) ≈ fmm.complex_add(real(z1), imag(z1), real(z2), imag(z2))[2]
    
    # subtraction
    @test real(z1-z2) ≈ fmm.complex_subtract(real(z1), imag(z1), real(z2), imag(z2))[1]
    @test imag(z1-z2) ≈ fmm.complex_subtract(real(z1), imag(z1), real(z2), imag(z2))[2]
    
    # multiplication
    @test real(z1*z2) ≈ fmm.complex_multiply(real(z1), imag(z1), real(z2), imag(z2))[1]
    @test imag(z1*z2) ≈ fmm.complex_multiply(real(z1), imag(z1), real(z2), imag(z2))[2]
    
    @test real(z1*z2*z3) ≈ fmm.complex_multiply(real(z1), imag(z1), real(z2), imag(z2), real(z3), imag(z3))[1]
    @test imag(z1*z2*z3) ≈ fmm.complex_multiply(real(z1), imag(z1), real(z2), imag(z2), real(z3), imag(z3))[2]
    
    # division
    @test real(z1/z2) ≈ fmm.complex_divide(real(z1), imag(z1), real(z2), imag(z2))[1]
    @test imag(z1/z2) ≈ fmm.complex_divide(real(z1), imag(z1), real(z2), imag(z2))[2]
    @test real(z1/z2) ≈ fmm.complex_divide_real(real(z1), imag(z1), real(z2), imag(z2))
    @test imag(z1/z2) ≈ fmm.complex_divide_imag(real(z1), imag(z1), real(z2), imag(z2))
end

@testset "direct" begin

    function V(xi, xj, mj; G=1)
        Rho_ij = xi - xj
        rho_ij = sqrt(Rho_ij' * Rho_ij)
        vij = G * mj / rho_ij
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

    fmm.direct!(mass)
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
    center, radius = fmm.center_radius((elements,); scale_radius = 1.00001)
    test_center = [0.65, 0.65, 0.55]
    test_radius = 0.5500055

    for i in 1:3
        @test isapprox(center[i], test_center[i]; atol=1e-4)
    end
    @test isapprox(radius, test_radius; atol=1e-4)

    # test branch! function
    expansion_order, n_per_branch, multipole_acceptance_criterion = 2, 1, 0.5
    tree = fmm.Tree((elements,); expansion_order, n_per_branch, shrink_recenter=false)

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
    fmm.cartesian_2_spherical!(this)
    for i in 1:3
        @test isapprox(this[i], that[i]; atol=1e-10)
    end
end

#= Renormalized, so these tests don't work anymore
# However, I may need to revert the normalization
# once I begin using multiple kernels
@testset "solid harmonics" begin

function Ylm(theta, phi, l, m)
    lm! = sqrt(factorial(big(l-abs(m)))/ factorial(big(l+abs(m))))
    plm = Plm(cos(theta), l, abs(m))
    eim = exp(im * m * phi)
    ylm = lm! * plm * eim
end

function regular_harmonic_manual(rho, theta, phi, p)
    reg_harmonics = Vector{Complex{Float64}}(undef,(p+1)^2)
    i = 1
    for l in 0:p
        for m in -l:l
            reg_harmonics[i] = Ylm(theta, phi, l, m) * rho^l
            i+=1
        end
    end
    return reg_harmonics
end

function irregular_harmonic_manual(rho, theta, phi, p)
    reg_harmonics = Vector{Complex{Float64}}(undef,(p+1)^2)
    i = 1
    for l in 0:p
        for m in -l:l
            reg_harmonics[i] = Ylm(theta, phi, l, m) / rho^(l+1)
            i+=1
        end
    end
    return reg_harmonics
end

rho = 1.2
alpha = pi/4 * 1.4
beta = pi/6 * 0.9
P = 3
# rh_exa = regular_harmonic(rho, alpha, beta, P+1)

rh_man = regular_harmonic_manual(rho, alpha, beta, P)
rh_fmm = zeros(Complex{Float64},length(rh_man))
rh_fmm_theta = zeros(Complex{Float64},length(rh_man))
fmm.regular_harmonic!(rh_fmm, rh_fmm_theta, rho, alpha, beta, P)

for i in 1:length(rh_man)
    @test isapprox(rh_man[i], rh_fmm[i]; atol=1e-11)
end

ih_man = irregular_harmonic_manual(rho, alpha, beta, P)
ih_fmm = zeros(Complex{Float64},length(rh_man))
fmm.irregular_harmonic!(ih_fmm, rho, alpha, beta, P)

for i in 1:length(ih_man)
    @test isapprox(ih_man[i], ih_fmm[i]; atol=1e-11)
end
end
=#

new_order = [4,5,2,3,1]

const new_order_index = [5,3,4,1,2]

# get the new index of mass_i as new_order_index[mass_i]

@testset "spherical P2M" begin
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
n_per_branch = 1
multipole_acceptance_criterion = 0.25
tree = fmm.Tree((system,); expansion_order, n_per_branch, shrink_recenter=false)

i_mass = 1
i_branch = 5 # use the first mass

harmonics = zeros(Float64,2,(expansion_order+1)^2)
fmm.B2M!(system, tree.branches[i_branch], new_order_index[i_mass]:new_order_index[i_mass], harmonics, tree.expansion_order)

center = tree.branches[i_branch].center

x_target = [10.1,-7.3,8.6]
target_potential = zeros(4)
fmm.M2B!(target_potential, x_target, i_branch, tree)

u_fmm = target_potential[1]

dx = x_target - xs[1,:]
u_check = ms[1] / sqrt(dx' * dx)

function Ylm(theta, phi, l, m)
    ylm = sqrt(factorial(big(l-abs(m)))/ factorial(big(l+abs(m)))) * Plm(cos(theta), l, abs(m)) * exp(im * m * phi)
end

function evaluate_biot_savart(x_source, x_target, q_source, P)
    v = 0.0
    i = 1

    for l in 0:P
        for m in -l:l

            v += q_source * x_source[1]^l / x_target[1]^(l+1) * real(Ylm(x_target[2], x_target[3], l, m) * conj(Ylm(x_source[2], x_source[3], l, m)))
            i += 1
        end
    end
    return v
end

x_source_sph = fmm.cartesian_2_spherical(xs[1,:] - center)
x_target_sph = fmm.cartesian_2_spherical(x_target - center)
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
tree = fmm.Tree((elements,); expansion_order, n_per_branch=1, shrink_recenter=false)

i_branch = 2 # contains 4th and 5th elements
i_branch_4 = 6 # use the fourth mass
# i_branch_5 = 7 # use the fifth mass
harmonics = zeros(Float64,2,(expansion_order+1)^2)
M = zeros(Float64,2,4)
# using only the 4th mass: (note it has been reordered)
fmm.B2M!(elements, tree.branches[i_branch_4], new_order_index[4]:new_order_index[4], harmonics, tree.expansion_order) # evaluate multipole coefficients
for child_branch in view(tree.branches, tree.branches[i_branch].branch_index) # translate coefficients to the center of branch 2
    fmm.M2M!(tree.branches[i_branch], child_branch, harmonics, M, expansion_order)
end

x_target = [8.3,1.4,-4.2]
target_potential = zeros(4)
target = x_target
fmm.M2B!(target_potential, target, i_branch, tree)
u_fmm = target_potential[1]

target_potential .*= 0
fmm.M2B!(target_potential, target, i_branch_4, tree)
u_fmm_no_x = target_potential[1]

dx = x_target - xs[4,:]
u_check = ms[4] / sqrt(dx'*dx)

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
tree = fmm.Tree((elements,); expansion_order, n_per_branch=1, shrink_recenter=false)

branch_i = 2 # contains two elements; 4 and 5
target_i = new_order_index[4]
source_i = new_order_index[1] # just needs to be farther away than the target to ensure convergence

dx_source = fmm.cartesian_2_spherical(elements[source_i,fmm.POSITION] - tree.branches[branch_i].center)
dx_target = fmm.cartesian_2_spherical(elements[target_i,fmm.POSITION] - tree.branches[branch_i].center)

local_coefficients_theta = zeros(Float64, 2, ((expansion_order+1)*(expansion_order+2))>>1)
local_coefficients_expanded = zeros(Float64, 2, (expansion_order+1)^2)
local_coefficients_expanded_theta = zeros(Float64, 2, (expansion_order+1)^2)
fmm.irregular_harmonic!(local_coefficients_expanded, dx_source..., expansion_order)
local_coefficients_expanded .*= ms[1]
regular_harmonics_expanded = zeros(Float64, 2, (expansion_order+1)^2)
fmm.regular_harmonic!(regular_harmonics_expanded, dx_target..., expansion_order)

fmm.B2L!(tree, branch_i, elements[source_i,fmm.POSITION], elements.bodies[source_i].strength)

harmonics = zeros(Float64, 2, (expansion_order+1)^2)
harmonics_theta = zeros(Float64, 2, (expansion_order+1)^2)
harmonics_theta_2 = zeros(Float64, 2, (expansion_order+1)^2)
workspace = zeros(3,4)
fmm.L2B!(elements, target_i:target_i, tree.branches[branch_i].local_expansion, fmm.DerivativesSwitch(), Val(expansion_order), tree.branches[branch_i].center, zeros(3), zeros(3,4), zeros(3,3,4), harmonics, harmonics_theta, harmonics_theta_2, workspace)
u_fmm = elements.potential[1,target_i]

dx_direct = xs[4,:] - xs[1,:]
u_check = 1 / sqrt(dx_direct' * dx_direct)
u_check *= ms[1]

u_man = sum(regular_harmonics_expanded[1,:].*local_coefficients_expanded[1,:] .+ regular_harmonics_expanded[2,:].*local_coefficients_expanded[2,:])
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
tree = fmm.Tree((elements,); expansion_order, n_per_branch=1, shrink_recenter=false)

# local coefficient at branch 2 due to mass 1
fmm.B2L!(tree, 2, elements[new_order_index[1],fmm.POSITION], elements.bodies[new_order_index[1]].strength)
# local_2 = deepcopy(tree.branches[2].local_expansion)

# check L2P now:
harmonics = zeros(Float64,2,(expansion_order+1)^2)
harmonics_theta = zeros(Float64,2,(expansion_order+1)^2)
harmonics_theta_2 = zeros(Float64,2,(expansion_order+1)^2)
workspace = zeros(3,4)
spherical_potential = zeros(52)
fmm.L2B!(elements, new_order_index[5]:new_order_index[5], tree.branches[2].local_expansion, fmm.DerivativesSwitch(), Val(expansion_order), tree.branches[2].center, zeros(3), zeros(3,4), zeros(3,3,4), harmonics, harmonics_theta, harmonics_theta_2, workspace)
u_fmm_no_x = elements.potential[1,new_order_index[5]]
elements.potential[1,new_order_index[5]] *= 0

# translate local expansion to branch 7 (mass 5)
fmm.L2L!(tree.branches[2], tree.branches[7], harmonics, zeros(eltype(tree.branches[1].multipole_expansion),2,4), tree.expansion_order)

local_coefficients_check = zeros(Float64,2, (expansion_order+1)^2)
dx_check, dy_check, dz_check = fmm.cartesian_2_spherical(elements[new_order_index[1],fmm.POSITION] - tree.branches[7].center)
fmm.irregular_harmonic!(local_coefficients_check, dx_check, dy_check, dz_check, expansion_order)
local_coefficients_check .*= ms[1]

# evaluate local expansion at mass 5
fmm.L2B!((elements,), tree.branches[7], (fmm.DerivativesSwitch(),), Val(expansion_order), zeros(3), zeros(3,4), zeros(3,3,4), harmonics, harmonics_theta, harmonics_theta_2, workspace)
u_fmm = elements.potential[1,new_order_index[5]]

dx_direct = elements[new_order_index[5],fmm.POSITION] - elements[new_order_index[1],fmm.POSITION]
u_check = ms[1] / sqrt(dx_direct' * dx_direct)

regular_harmonics = zeros(Float64,2, (expansion_order+1)^2)
dx_target = fmm.cartesian_2_spherical(elements[new_order_index[5],fmm.POSITION] - tree.branches[7].center)
fmm.regular_harmonic!(regular_harmonics, dx_target..., expansion_order)
#u_man = real(sum(regular_harmonics' * local_coefficients_check))
u_man = sum(regular_harmonics[1,:].*local_coefficients_check[1,:] .+ regular_harmonics[2,:].*local_coefficients_check[2,:])

@test isapprox(u_check, u_man; atol=1e-12)
@test isapprox(u_check, u_fmm_no_x; atol=1e-12)
@test isapprox(u_check, u_fmm; atol=1e-12)

end

@testset "spherical: M2L" begin
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
tree = fmm.Tree((elements,); expansion_order, n_per_branch=1, shrink_recenter=false)

i_branch_multipole = 7 # mass 5
i_branch_local = 5 # mass 1
harmonics = zeros(Float64,2, (expansion_order+1)^2)
harmonics_theta = zeros(Float64,2, (expansion_order+1)^2)
harmonics_theta_2 = zeros(Float64,2, (expansion_order+1)^2)
workspace = zeros(3,4)

fmm.B2M!(elements, tree.branches[i_branch_multipole], new_order_index[5]:new_order_index[5], harmonics, tree.expansion_order)

# # test Multipole # checks out
# dx_mp = xs[5,:] - tree.branches[i_branch_multipole].center
# fmm.cartesian_2_spherical!(dx_mp)
# fmm.regular_harmonic!(harmonics, dx_mp..., expansion_order)
# multipole_check =  harmonics * ms[5]
# dx_mp = xs[1,:] - tree.branches[i_branch_multipole].center
# fmm.cartesian_2_spherical!(dx_mp)
# fmm.irregular_harmonic!(harmonics, dx_mp..., expansion_order)
# u_check_mp = real(sum(harmonics' * multipole_check))

# @show tree.branches[i_branch_multipole].multipole_expansion[1][1:10] multipole_check[1:10]
# ###
m2l_harmonics = zeros(eltype(tree.branches[1].multipole_expansion), 2, (expansion_order<<1 + 1)*(expansion_order<<1 + 1))
L = zeros(eltype(tree.branches[1].local_expansion), 2, 4)
fmm.M2L!(tree.branches[i_branch_local], tree.branches[i_branch_multipole], m2l_harmonics, L, expansion_order)
fmm.L2B!(elements, new_order_index[1]:new_order_index[1], tree.branches[i_branch_local].local_expansion, fmm.DerivativesSwitch(), Val(expansion_order), tree.branches[i_branch_local].center, zeros(3), zeros(3,4), zeros(3,3,4), harmonics, harmonics_theta, harmonics_theta_2, workspace)
u_fmm = elements.potential[1,new_order_index[1]]

local_exp = tree.branches[i_branch_local].local_expansion[1]

# test local
# dx_l = xs[1,:] - tree.branches[i_branch_local].center
# fmm.regular_harmonic!(harmonics, dx_l..., expansion_order)
# u_check_local = real(sum(harmonics' * local_exp))

dx_direct = elements[new_order_index[1],fmm.POSITION] - elements[new_order_index[5], fmm.POSITION]
u_direct = elements.bodies[new_order_index[5]].strength[1] / sqrt(dx_direct' * dx_direct)

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
multipole_acceptance_criterion = 0.5
tree = fmm.Tree((elements,); expansion_order, n_per_branch=1, shrink_recenter=false)

# perform upward pass
fmm.upward_pass_singlethread!(tree.branches, (elements,), expansion_order)

# m6 = tree.branches[6].multipole_expansion
target = [4.1,2.2,3.4]
dx_direct_6 = target - elements[new_order_index[4],fmm.POSITION]
u_direct_6 = ms[4] / sqrt(dx_direct_6' * dx_direct_6)

mass_target_potential = zeros(4)
mass_target = target
fmm.M2B!(mass_target_potential, mass_target, 6, tree)
u_fmm_6 = mass_target_potential[1]

# add branches 6 and 7
dx_direct_7 = target - elements[new_order_index[5],fmm.POSITION]
u_direct_67 = u_direct_6 + ms[5] / sqrt(dx_direct_7' * dx_direct_7)

# reset target potential
mass_target_potential *= 0

# use summed multipole expansion from branches 6 and 7 (summed at 2)
fmm.M2B!(mass_target_potential, mass_target, 2, tree)
u_fmm_67 = mass_target_potential[1]

# perform horizontal pass
m2l_list, direct_list = fmm.build_interaction_lists(tree.branches, tree.branches, multipole_acceptance_criterion, true, true, true)
fmm.nearfield_singlethread!((elements,), tree.branches, (fmm.DerivativesSwitch(),), (elements,), tree.branches, direct_list)
fmm.horizontal_pass_singlethread!(tree.branches, tree.branches, m2l_list, expansion_order)

# consider the effect on branch 3 (mass 2)
elements.potential[i_POTENTIAL,new_order_index[2]] .*= 0 # reset potential at mass 2
# P2P is performed from branches 3 (mass 2), 4 (mass 3), and 5 (mass 1) to branch 3
# elements.direct!(elements.potential[i_POTENTIAL,new_order_index[2]], elements.bodies[i_POSITION,new_order_index[2]], elements.bodies[:,new_order_index[1]])
# elements.direct!(elements.potential[i_POTENTIAL,new_order_index[2]], elements.bodies[i_POSITION,new_order_index[2]], elements.bodies[:,new_order_index[2]])
# elements.direct!(elements.potential[i_POTENTIAL,new_order_index[2]], elements.bodies[i_POSITION,new_order_index[2]], elements.bodies[:,new_order_index[3]])
fmm.P2P!((elements,), tree.branches[3], (fmm.DerivativesSwitch(),), (elements,), tree.branches[3])
fmm.P2P!((elements,), tree.branches[3], (fmm.DerivativesSwitch(),), (elements,), tree.branches[4])
fmm.P2P!((elements,), tree.branches[3], (fmm.DerivativesSwitch(),), (elements,), tree.branches[5])
u_fmm_123 = elements.potential[i_POTENTIAL[1],new_order_index[2]]

dx_12 = elements[new_order_index[2],fmm.POSITION] - elements[new_order_index[1],fmm.POSITION]
u_direct_12 = elements.bodies[new_order_index[1]].strength[1] / sqrt(dx_12' * dx_12)
u_direct_22 = 0.0
dx_32 = elements[new_order_index[2],fmm.POSITION] - elements[new_order_index[3],fmm.POSITION]
u_direct_32 = elements.bodies[new_order_index[3]].strength[1] / sqrt(dx_32' * dx_32)

u_direct_123 = u_direct_12 + u_direct_22 + u_direct_32

# M2L is performed from branches 6, 7 to branch 3 (containing mass 2)
harmonics = zeros(Float64,2, (expansion_order+1)^2)
harmonics_theta = zeros(Float64,2, (expansion_order+1)^2)
harmonics_theta_2 = zeros(Float64,2, (expansion_order+1)^2)
workspace = zeros(3,4)
fmm.L2B!((elements,), tree.branches[3], (fmm.DerivativesSwitch(),), Val(expansion_order), zeros(3), zeros(3,4), zeros(3,3,4), harmonics, harmonics_theta, harmonics_theta_2, workspace)
u_fmm_12345 = elements.potential[i_POTENTIAL[1],new_order_index[2]]

dx_42 = elements[new_order_index[4],fmm.POSITION] - elements[new_order_index[2],fmm.POSITION]
u_direct_42 = elements.bodies[new_order_index[4]].strength[1] / sqrt(dx_42' * dx_42)
dx_52 = elements[new_order_index[5],fmm.POSITION] - elements[new_order_index[2],fmm.POSITION]
u_direct_52 = elements.bodies[new_order_index[5]].strength[1] / sqrt(dx_52' * dx_52)

u_direct_12345 = u_direct_123 + u_direct_42 + u_direct_52

@test isapprox(u_direct_123, u_fmm_123; atol=1e-12)

@test isapprox(u_direct_12345, u_fmm_12345; atol=1e-12)

# reset potentials
elements.potential .*= 0

# run fmm (reset potentials with reset_tree flag)
fmm.fmm!(tree, (elements,); multipole_acceptance_criterion=multipole_acceptance_criterion, reset_tree=true)
u_fmm = deepcopy(elements.potential[1,:])

elements.potential .= 0.0

fmm.direct!((elements,))

u_direct = deepcopy(elements.potential[1,:])

for i in 1:fmm.get_n_bodies(elements)
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
end

@testset "chain rule" begin

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

# """
# dr_i/dx_j
# """
# function drdx(r,theta,phi)
#     # derivatives = [
#     #     sin(theta)*cos(phi) cos(theta)*cos(phi)/r -sin(phi)/r/sin(theta);
#     #     sin(theta)*sin(phi) cos(theta)*sin(theta)/r cos(phi)/r/sin(theta);
#     #     cos(theta) -sin(theta)/r 0
#     # ]
#     derivatives = [
#         sin(theta)*cos(phi) sin(theta)*sin(phi) cos(theta);
#         cos(theta)*cos(phi)/r cos(theta)*sin(phi)/r -sin(theta)/r;
#         -sin(phi)/r/sin(theta) cos(phi)/r/sin(theta) 0
#     ]
#     return derivatives
# end

function derivatives_2_cartesian(first_derivatives, second_derivatives, r, theta, phi)
    @assert size(spherical_derivatives) == (3,3)
    d2_unit_vec = d2rdx2(r, theta, phi)
    d_unit_vec = drdx(r,theta,phi)
    cartesian_derivatives = zeros(3,3)
    for k in 1:3
        cartesian_derivatives .+= d2_unit_vec[:,:,k] .* first_derivatives[k]
    end
    cartesian_derivatives .+= d_unit_vec * second_derivatives * transpose(d_unit_vec)
    return cartesian_derivatives
end

function cs_derivative(func, i, args; step=1e-25)
    args[i] += step*im
    return imag(func(args))/step
end

function simple(X)
    r = X[1]
    theta = X[2]
    phi = X[3]
    return 3*r^4*cos(theta) + 2*theta - sin(phi)
end

function cartesian_2_spherical(x,y,z; vec=false)
    r = sqrt(x^2+y^2+z^2)
    theta = acos(z/r)
    phi = atan(y/x)
    vec && return [r,theta,phi]
    return r, theta, phi
end

function spherical_2_cartesian(r,theta,phi; vec=false)
    x = r*sin(theta)*cos(phi)
    y = r*sin(theta)*sin(phi)
    z = r*cos(theta)
    vec && return [x,y,z]
    return x,y,z
end

function simple_cart(X)
    args = cartesian_2_spherical(X...; vec=true)
    return simple(args)
end

r0, theta0, phi0 = rand(3)

dsdr(r,theta,phi) = 12 * r^3 * cos(theta)
dsdt(r,theta,phi) = 3*r^4*-sin(theta) + 2
dsdp(r,theta,phi) = -cos(phi)
d2sdr2(r,theta,phi) = 36 * r^2 * cos(theta)
d2sdrdt(r,theta,phi) = -12*r^3*sin(theta)
d2sdrdp(r,theta,phi) = 0.0
d2sdt2(r,theta,phi) = -3*r^4*cos(theta)
d2sdtdp(r,theta,phi) = 0.0
d2sdp2(r,theta,phi) = sin(phi)

testd = cs_derivative(simple,1,Complex{Float64}[r0,theta0,phi0])
@test isapprox(testd, dsdr(r0,theta0,phi0);atol=1e-12)

function second_derivative(func,i,j,args; step=1e-8)
    first_func(args) = cs_derivative(func,i,args)
    mask = zeros(length(args))
    mask[j] += step
    sec = (first_func(args + mask) - first_func(args))/step
    return sec
end

function fd_hessian(func, args)
    simple_jacobian = zeros(3,3)
    for j in 1:3
        for i in 1:3
            simple_jacobian[j,i] = second_derivative(func,i,j,convert(Vector{Complex{Float64}},args))
        end
    end
    return simple_jacobian
end
simple_jacobian = fd_hessian(simple,[r0,theta0,phi0])
simple_jacobian_anal = [
    d2sdr2(r0,theta0,phi0) d2sdrdt(r0,theta0,phi0) d2sdrdp(r0,theta0,phi0);
    d2sdrdt(r0,theta0,phi0) d2sdt2(r0,theta0,phi0) d2sdtdp(r0,theta0,phi0);
    d2sdrdp(r0,theta0,phi0) d2sdtdp(r0,theta0,phi0) d2sdp2(r0,theta0,phi0);
]

for i in 1:length(simple_jacobian)
    @test isapprox(simple_jacobian[i], simple_jacobian_anal[i]; atol=1e-6)
end

# test chain rule
args = [r0,theta0,phi0]
spherical_grad = [cs_derivative(simple,i,convert(Vector{Complex{Float64}},args)) for i in 1:3]
spherical_hessian = fd_hessian(simple, args)
d2rkdxidxj = d2rdx2(args...)
drjdxi = drdx(args...)

cartesian_hessian = zeros(3,3)
for k in 1:3
    cartesian_hessian .+= d2rkdxidxj[:,:,k] * spherical_grad[k]
end

cartesian_hessian .+= drjdxi * spherical_hessian * transpose(drjdxi)

cartesian_hessian_fd = fd_hessian(simple_cart, spherical_2_cartesian(args...;vec=true))

for i in 1:length(cartesian_hessian)
    @test isapprox(cartesian_hessian[i], cartesian_hessian_fd[i]; rtol=1e-3)
end

# now test FMM function
potential_hessian = zeros(3,3,4)
potential_jacobian = zeros(3,4)
for i in 1:4
    potential_hessian[:,:,i] .= spherical_hessian
    potential_jacobian[:,i] .= spherical_grad
end
workspace = zeros(3,4)

fmm.spherical_2_cartesian!(potential_jacobian, potential_hessian, workspace, r0, theta0, phi0, fmm.DerivativesSwitch())

for i in 1:3
    for j in 1:3
        @test isapprox(cartesian_hessian[i,j], potential_hessian[i,j,1]; rtol=1e-14)
    end
end

#--- first and second derivatives of the regular solid harmonics in theta ---#

function test_rh_i(i; P=7, r=rand(), theta_test = rand()*pi, phi = rand()*2*pi)
    P = 7
    r = rand()
    r < 1e-2 && (r = 0.4)
    theta_test = rand() * pi
    isapprox(theta_test, 0.0; atol=1e-1) && (theta_test = 1e-1)
    isapprox(theta_test, pi; atol=1e-1) && (theta_test = pi - 1e-1)
    phi = rand()*2*pi

    function get_regular_harmonic_test(theta::TF; r=r, phi=phi, P=P) where TF
        _, _, _, derivative_harmonics, derivative_harmonics_theta, derivative_harmonics_theta_2, _ = fmm.preallocate_l2b(TF, TF, Val(P))
        fmm.regular_harmonic!(derivative_harmonics, derivative_harmonics_theta, derivative_harmonics_theta_2, r, theta, phi, P)
        return derivative_harmonics[:,i]
    end

    function get_regular_harmonic_theta_test(theta::TF; r=r, phi=phi, P=P) where TF
        _, _, _, derivative_harmonics, derivative_harmonics_theta, derivative_harmonics_theta_2, _ = fmm.preallocate_l2b(TF, TF, Val(P))
        fmm.regular_harmonic!(derivative_harmonics, derivative_harmonics_theta, derivative_harmonics_theta_2, r, theta, phi, P)
        return derivative_harmonics_theta[:,i]
    end

    function get_regular_harmonic_theta_2_test(theta::TF; r=r, phi=phi, P=P) where TF
        _, _, _, derivative_harmonics, derivative_harmonics_theta, derivative_harmonics_theta_2, _ = fmm.preallocate_l2b(TF, TF, Val(P))
        fmm.regular_harmonic!(derivative_harmonics, derivative_harmonics_theta, derivative_harmonics_theta_2, r, theta, phi, P)
        return derivative_harmonics_theta_2[:,i]
    end

    function check_regular_harmonic_theta(theta)
        return ForwardDiff.derivative(get_regular_harmonic_test, theta)
    end

    function check_regular_harmonic_theta_2(theta)
        return ForwardDiff.derivative((t) -> ForwardDiff.derivative(get_regular_harmonic_test, t), theta)
    end

    rh = get_regular_harmonic_test(theta_test)
    rh_p = get_regular_harmonic_theta_test(theta_test)
    rh_pp = get_regular_harmonic_theta_2_test(theta_test)
    rh_p_check = check_regular_harmonic_theta(theta_test)
    rh_pp_check = check_regular_harmonic_theta_2(theta_test)

    @test isapprox(rh_p[1], rh_p_check[1]; atol=1e-12)
    @test isapprox(rh_p[2], rh_p_check[2]; atol=1e-12)
    @test isapprox(rh_pp[1], rh_pp_check[1]; atol=1e-12)
    @test isapprox(rh_pp[2], rh_pp_check[2]; atol=1e-12)

    return nothing
end

P = 7
r = rand()
theta_test = rand() * pi
phi = rand()*2*pi

for i in 1:36
    test_rh_i(i; P, r, theta_test, phi)
end

#--- hessian of spherical coordinates w.r.t. cartesian coordinates ---#

# r coordinate
function get_dr2dxidxj(xvec)
    x, y, z = xvec
    rho, theta, phi = fmm.cartesian_2_spherical(x, y, z)
    s_theta, c_theta = sincos(theta)
    s_phi, c_phi = sincos(phi)
    dr2dxidxj = @SMatrix [
        (1-c_phi^2 * s_theta^2)/rho -s_theta^2*c_phi*s_phi/rho -s_theta*c_phi*c_theta/rho;
        (-s_theta^2*c_phi*s_phi)/rho (1-s_theta^2*s_phi^2)/rho -s_theta*s_phi*c_theta/rho;
        -s_theta*c_phi*c_theta/rho -s_theta*s_phi*c_theta/rho s_theta^2/rho
    ]
    return dr2dxidxj
end

function check_dr2dxidxj(xvec)
    return ForwardDiff.jacobian((x) -> ForwardDiff.gradient((xx) -> sqrt(xx'*xx), x), xvec)
end

function test_dr2dxidxj(xvec)
    dr2dxidxj = get_dr2dxidxj(xvec)
    dr2dxidxj_check = check_dr2dxidxj(xvec)
    for i in eachindex(dr2dxidxj)
        @test isapprox(dr2dxidxj[i], dr2dxidxj_check[i]; atol=1e-12)
    end
end

for _ in 1:10
    test_dr2dxidxj(rand(3))
    test_dr2dxidxj(rand(3)*10)
end

# theta coordinate
function get_dtheta2dxidxj(xvec)
    x, y, z = xvec
    rho, theta, phi = fmm.cartesian_2_spherical(x, y, z)
    s_theta, c_theta = sincos(theta)
    s_phi, c_phi = sincos(phi)
    dtheta2dxidxj = @SMatrix [ # this works? much better at least
        c_theta/s_theta*(1-c_phi^2*(1+2*s_theta^2))/rho^2 -c_theta/s_theta*s_phi*c_phi*(1+2*s_theta^2)/rho^2 c_phi*(1-2*c_theta^2)/rho^2;
        -c_theta/s_theta*s_phi*c_phi*(1+2*s_theta^2)/rho^2 c_theta/s_theta*(1-s_phi^2*(1+2*s_theta^2))/rho^2 (2*s_theta^2-1)/rho^2*s_phi;
        c_phi*(1-2*c_theta^2)/rho^2 (2*s_theta^2-1)/rho^2*s_phi 2*s_theta*c_theta/rho^2
    ]
    return dtheta2dxidxj
end

function check_dtheta2dxidxj(xvec)
    return ForwardDiff.jacobian((x) -> ForwardDiff.gradient((xx) -> acos(xx[3]/sqrt(xx'*xx)), x), xvec)
end

function test_dtheta2dxidxj(xvec)
    dtheta2dxidxj = get_dtheta2dxidxj(xvec)
    dtheta2dxidxj_check = check_dtheta2dxidxj(xvec)
    for i in eachindex(dtheta2dxidxj)
        @test isapprox(dtheta2dxidxj[i], dtheta2dxidxj_check[i]; atol=1e-12)
    end
end

for _ in 1:10
    test_dtheta2dxidxj(rand(3))
    test_dtheta2dxidxj(rand(3)*10)
end

# phi coordinate
function get_dphi2dxidxj(xvec)
    x, y, z = xvec
    rho, theta, phi = fmm.cartesian_2_spherical(x, y, z)
    s_theta, c_theta = sincos(theta)
    s_phi, c_phi = sincos(phi)
    dphi2dxidxj = @SMatrix [
        2*c_phi*s_phi/rho^2/s_theta^2 (2*s_phi^2-1)/rho^2/s_theta^2 0;
        (2*s_phi^2-1)/rho^2/s_theta^2 -2*s_phi*c_phi/rho^2/s_theta^2 0;
        0 0 0
    ]
    return dphi2dxidxj
end

function check_dphi2dxidxj(xvec)
    return ForwardDiff.jacobian((x) -> ForwardDiff.gradient((xx) -> atan(xx[2]/xx[1]), x), xvec)
end

function test_dphi2dxidxj(xvec)
    dphi2dxidxj = get_dphi2dxidxj(xvec)
    dphi2dxidxj_check = check_dphi2dxidxj(xvec)
    for i in eachindex(dphi2dxidxj)
        @test isapprox(dphi2dxidxj[i], dphi2dxidxj_check[i]; atol=1e-12)
    end
end

for _ in 1:10
    test_dphi2dxidxj(rand(3))
    test_dphi2dxidxj(rand(3)*10)
end

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
fmm.direct!((vortexparticles,))
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
n_per_branch = 1
x_branch_1 = fmm.SVector{3}([0.0,0,0])
branch_1 = fmm.MultiBranch(SVector{1}([1:2]), 2, 2:3, 0, x_branch_1, 1/8, fmm.initialize_expansion(expansion_order), fmm.initialize_expansion(expansion_order), fmm.initialize_harmonics(expansion_order, Float64), fmm.initialize_ML(expansion_order, Float64), ReentrantLock())
x_branch_2 = fmm.SVector{3}(xs[:,1] .+ [0.01, 0.02, -0.03])
branch_2 = fmm.MultiBranch(SVector{1}([1:1]), 0, 3:2, 1, x_branch_2, 1/8, fmm.initialize_expansion(expansion_order), fmm.initialize_expansion(expansion_order), fmm.initialize_harmonics(expansion_order, Float64), fmm.initialize_ML(expansion_order, Float64), ReentrantLock())
x_branch_3 = fmm.SVector{3}(xs[:,2] .+ [0.02, -0.04, 0.01])
branch_3 = fmm.MultiBranch(SVector{1}([2:2]), 0, 3:2, 1, x_branch_3, 1/8, fmm.initialize_expansion(expansion_order), fmm.initialize_expansion(expansion_order), fmm.initialize_harmonics(expansion_order, Float64), fmm.initialize_ML(expansion_order, Float64), ReentrantLock())

# using FMM
# tree = fmm.Tree(branches, [expansion_order], n_per_branch, B2M!, P2P!)
dummy_index = (zeros(Int64,fmm.get_n_bodies(vortexparticles)),)
dummy_leaf_index = collect(1:3)
# dummy_cost_parameter = fmm.dummy_direct_cost_estimate((vortexparticles,), n_per_branch)
tree = fmm.MultiTree([branch_1, branch_2, branch_3], [1:1,2:3], dummy_leaf_index, dummy_index, dummy_index, (deepcopy(vortexparticles),), Val(expansion_order), n_per_branch)#, dummy_cost_parameter)
harmonics = zeros(Complex{Float64},2,(expansion_order+1)^2)
fmm.B2M!(branch_2, (vortexparticles,), harmonics, Val(expansion_order))
fmm.B2M!(branch_3, (vortexparticles,), harmonics, Val(expansion_order))
# @show tree.branches[2].multipole_expansion[2:4,:] # checks out

m2l_harmonics = zeros(eltype(tree.branches[1].multipole_expansion), 2, (expansion_order<<1 + 1)*(expansion_order<<1 + 1))
L = zeros(eltype(tree.branches[1].local_expansion), 2, 4)
fmm.M2L!(tree.branches[2], tree.branches[3], m2l_harmonics, L, Val(expansion_order))
fmm.M2L!(tree.branches[3], tree.branches[2], m2l_harmonics, L, Val(expansion_order))
# @show tree.branches[2].local_expansion
harmonics = zeros(Float64,2, (expansion_order+1)^2)
harmonics_theta = zeros(Float64,2, (expansion_order+1)^2)
harmonics_theta_2 = zeros(Float64,2, (expansion_order+1)^2)
workspace = zeros(3,4)
fmm.L2B!((vortexparticles,), tree.branches[2], (fmm.DerivativesSwitch(),), Val(expansion_order), zeros(3), zeros(3,4), zeros(3,3,4), harmonics, harmonics_theta, harmonics_theta_2, workspace)
fmm.L2B!((vortexparticles,), tree.branches[3], (fmm.DerivativesSwitch(),), Val(expansion_order), zeros(3), zeros(3,4), zeros(3,3,4), harmonics, harmonics_theta, harmonics_theta_2, workspace)
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
fmm.direct!((vortex_particles,))
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
n_per_branch = 1
x_branch_1 = SVector{3}((bodies[1:3,1] + bodies[1:3,2])/2)
branch_1 = fmm.MultiBranch(SVector{1}([1:2]), 2, 2:3, 0, x_branch_1, 1/8, fmm.initialize_expansion(expansion_order), fmm.initialize_expansion(expansion_order), fmm.initialize_harmonics(expansion_order, Float64), fmm.initialize_ML(expansion_order, Float64), ReentrantLock())
x_branch_2 = SVector{3}(bodies[1:3,1] .+ [0.01, 0.02, -0.03])
branch_2 = fmm.MultiBranch(SVector{1}([1:1]), 0, 3:2, 1, x_branch_2, 1/8, fmm.initialize_expansion(expansion_order), fmm.initialize_expansion(expansion_order), fmm.initialize_harmonics(expansion_order, Float64), fmm.initialize_ML(expansion_order, Float64), ReentrantLock())
x_branch_3 = SVector{3}(bodies[1:3,2] .+ [0.02, -0.04, 0.01])
branch_3 = fmm.MultiBranch(SVector{1}([2:2]), 0, 3:2, 1, x_branch_2, 1/8, fmm.initialize_expansion(expansion_order), fmm.initialize_expansion(expansion_order), fmm.initialize_harmonics(expansion_order, Float64), fmm.initialize_ML(expansion_order, Float64), ReentrantLock())

dummy_index = (zeros(Int,length(vortex_particles.bodies)),)
dummy_leaf_index = collect(1:3)
# dummy_cost_parameter = fmm.dummy_direct_cost_estimate((vortex_particles,), n_per_branch)
tree = fmm.MultiTree([branch_1, branch_2, branch_3], [1:1,2:3], dummy_leaf_index, dummy_index, dummy_index, (deepcopy(vortex_particles),), Val(expansion_order), n_per_branch)#, dummy_cost_parameter)
# fmm.B2M!(tree, vortex_particles, 2)
# fmm.B2M!(tree, vortex_particles, 3)

# fmm.M2L!(tree, 2, 3)
# fmm.M2L!(tree, 3, 2)
# fmm.L2B!(tree, vortex_particles, 2)
# fmm.L2B!(tree, vortex_particles, 3)
fmm.fmm!(tree, (vortex_particles,); multipole_acceptance_criterion=0.5, unsort_bodies=false)
update_velocity_stretching!(vortex_particles)

psis_fmm = deepcopy(vortex_particles.potential[2:4,:])
for i in 1:length(psis_fmm)
    @test isapprox(psis_fmm[i], psis[i]; atol=1e-10)
end
# hessians_fmm = deepcopy(reshape(vortex_particles.potential[i_POTENTIAL_HESSIAN[10:end],:],3,3,3,2))
# for i in 1:length(hessians)
#     @test isapprox(hessians_fmm[i], hessians[i]; atol=1e-8)
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
fmm.direct!((vortex_particles,))
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
n_per_branch = 1

tree = fmm.Tree((vortex_particles,); expansion_order, n_per_branch, shrink_recenter=false)

fmm.fmm!(tree, (vortex_particles,); multipole_acceptance_criterion=0.5)

update_velocity_stretching!(vortex_particles)

psis_fmm = deepcopy(vortex_particles.potential[2:4,:])
for i in 1:length(psis_fmm)
    @test isapprox(psis_fmm[i], psis[i]; rtol=1e-12)
end
# hessians_fmm = deepcopy(reshape(vortex_particles.potential[i_POTENTIAL_HESSIAN[10:end],:],3,3,3,3))
# for i in 1:length(hessians)
#     @test isapprox(hessians_fmm[i], hessians[i]; rtol=1e-12)
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

expansion_order, n_per_branch, multipole_acceptance_criterion = 14, 100, 0.31
n_bodies = 5000
shrink_recenter = true
seed = 123
validation_system = generate_gravitational(seed, n_bodies; radius_factor=0.1)
fmm.direct!(validation_system)
validation_potential = validation_system.potential[1,:]

system3 = generate_gravitational(seed, n_bodies; radius_factor=0.1)
fmm.fmm!(system3; expansion_order=expansion_order, n_per_branch=n_per_branch, multipole_acceptance_criterion=multipole_acceptance_criterion, nearfield=true, farfield=true, unsort_bodies=true)
potential3 = system3.potential[1,:]
@test isapprox(maximum(abs.(potential3 - validation_potential)), 0.0; atol=1e-9)

system4 = generate_gravitational(seed, n_bodies; radius_factor=0.1)
fmm.fmm!((system4,); expansion_order=expansion_order, n_per_branch=n_per_branch, multipole_acceptance_criterion=multipole_acceptance_criterion, nearfield=true, farfield=true, unsort_bodies=true)
potential4 = system4.potential[1,:]
@test isapprox(maximum(abs.(potential4 - validation_potential)), 0.0; atol=1e-9)

system5 = generate_gravitational(seed, n_bodies; radius_factor=0.1)
fmm.fmm!(system5, system5; expansion_order=expansion_order, n_per_branch_source=n_per_branch, n_per_branch_target=n_per_branch, multipole_acceptance_criterion=multipole_acceptance_criterion, nearfield=true, farfield=true, unsort_source_bodies=true, unsort_target_bodies=true)
potential5 = system5.potential[1,:]
@test isapprox(maximum(abs.(potential5 - validation_potential)), 0.0; atol=1e-9)

system6 = generate_gravitational(seed, n_bodies; radius_factor=0.1)
fmm.fmm!((system6,), system6; expansion_order=expansion_order, n_per_branch_source=n_per_branch, n_per_branch_target=n_per_branch, multipole_acceptance_criterion=multipole_acceptance_criterion, nearfield=true, farfield=true, unsort_source_bodies=true, unsort_target_bodies=true)
potential6 = system6.potential[1,:]
@test isapprox(maximum(abs.(potential6 - validation_potential)), 0.0; atol=1e-9)

system7 = generate_gravitational(seed, n_bodies; radius_factor=0.1)
fmm.fmm!((system7,), (system7,); expansion_order=expansion_order, n_per_branch_source=n_per_branch, n_per_branch_target=n_per_branch, multipole_acceptance_criterion=multipole_acceptance_criterion, nearfield=true, farfield=true, unsort_source_bodies=true, unsort_target_bodies=true)
potential7 = system7.potential[1,:]
@test isapprox(maximum(abs.(potential7 - validation_potential)), 0.0; atol=1e-9)

end

@testset "sortwrapper" begin

expansion_order, n_per_branch, multipole_acceptance_criterion = 14, 100, 0.31
n_bodies = 5000
shrink_recenter = true
seed = 123
validation_system = generate_gravitational(seed, n_bodies; radius_factor=0.1)
fmm.direct!(validation_system)
validation_potential = validation_system.potential[1,:]

system8 = fmm.SortWrapper(generate_gravitational(seed, n_bodies; radius_factor=0.1))
fmm.fmm!(system8; expansion_order=expansion_order, n_per_branch=n_per_branch, multipole_acceptance_criterion=multipole_acceptance_criterion, nearfield=true, farfield=true, unsort_bodies=true)
potential8 = system8.system.potential[1,:]
@test isapprox(maximum(abs.(potential8 - validation_potential)), 0.0; atol=1e-9)

system9 = fmm.SortWrapper(generate_gravitational(seed, n_bodies; radius_factor=0.1))
fmm.fmm!((system9,); expansion_order=expansion_order, n_per_branch=n_per_branch, multipole_acceptance_criterion=multipole_acceptance_criterion, nearfield=true, farfield=true, unsort_bodies=true)
potential9 = system9.system.potential[1,:]
@test isapprox(maximum(abs.(potential9 - validation_potential)), 0.0; atol=1e-9)

system10 = fmm.SortWrapper(generate_gravitational(seed, n_bodies; radius_factor=0.1))
fmm.fmm!(system10, system10; expansion_order=expansion_order, n_per_branch_source=n_per_branch, n_per_branch_target=n_per_branch, multipole_acceptance_criterion=multipole_acceptance_criterion, nearfield=true, farfield=true, unsort_source_bodies=true, unsort_target_bodies=true)
potential10 = system10.system.potential[1,:]
@test isapprox(maximum(abs.(potential10 - validation_potential)), 0.0; atol=1e-9)

system11 = fmm.SortWrapper(generate_gravitational(seed, n_bodies; radius_factor=0.1))
fmm.fmm!((system11,), system11; expansion_order=expansion_order, n_per_branch_source=n_per_branch, n_per_branch_target=n_per_branch, multipole_acceptance_criterion=multipole_acceptance_criterion, nearfield=true, farfield=true, unsort_source_bodies=true, unsort_target_bodies=true)
potential11 = system11.system.potential[1,:]
@test isapprox(maximum(abs.(potential11 - validation_potential)), 0.0; atol=1e-9)

system12 = fmm.SortWrapper(generate_gravitational(seed, n_bodies; radius_factor=0.1))
fmm.fmm!((system12,), (system12,); expansion_order=expansion_order, n_per_branch_source=n_per_branch, n_per_branch_target=n_per_branch, multipole_acceptance_criterion=multipole_acceptance_criterion, nearfield=true, farfield=true, unsort_source_bodies=true, unsort_target_bodies=true)
potential12 = system12.system.potential[1,:]
@test isapprox(maximum(abs.(potential12 - validation_potential)), 0.0; atol=1e-9)

end

@testset "b2m source and dipole panels" begin

# source panel expansion
strength = 0.1428909901797533
R0_global = [-0.11530793537122011, 0.1997192025788218, -0.9730448705798238]
R0 = [-0.48862125790260025, -0.1735941199525583, -1.8441092898197107]
Ru = [0.04707023255275322, -0.10579806212499904, -0.020193487162119217]
Rv = [-0.039004202054645, -0.02833821156361363, 0.0]
normal = rand(SVector{3,Float64})
qnm_prev = [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]
jnm_prev = [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]
inm_prev = [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]
multipole_expansion = zeros(2,4,120)
multipole_expansion_test = [3.15347369905358e-5 0.0 0.0 0.0; 0.0 0.0 0.0 0.0;;; 5.836576687179473e-5 0.0 0.0 0.0; 0.0 0.0 0.0 0.0;;; -7.661878069585717e-6 0.0 0.0 0.0; 3.442114476513879e-6 0.0 0.0 0.0;;; 5.176954094092041e-5 0.0 0.0 0.0; 0.0 0.0 0.0 0.0;;; -1.4179736287285782e-5 0.0 0.0 0.0; 6.372425519286581e-6 0.0 0.0 0.0;;; 7.42180196746523e-7 0.0 0.0 0.0; -8.340091929224401e-7 0.0 0.0 0.0;;; 2.9171240607705508e-5 0.0 0.0 0.0; 0.0 0.0 0.0 0.0;;; -1.2848310937891748e-5 0.0 0.0 0.0; 5.7763762881869375e-6 0.0 0.0 0.0;;; 1.3731914469928547e-6 0.0 0.0 0.0; -1.5438750448677817e-6 0.0 0.0 0.0;;; -2.981931723472997e-8 0.0 0.0 0.0; 9.412368214502325e-8 0.0 0.0 0.0;;; 1.1616879034380539e-5 0.0 0.0 0.0; 0.0 0.0 0.0 0.0;;; -7.589495677738776e-6 0.0 0.0 0.0; 3.4136975727862394e-6 0.0 0.0 0.0;;; 1.2526909809559565e-6 0.0 0.0 0.0; -1.409195559129566e-6 0.0 0.0 0.0;;; -5.512043533657875e-8 0.0 0.0 0.0; 1.742126360703301e-7 0.0 0.0 0.0;;; -7.340115650494341e-10 0.0 0.0 0.0; -6.495176725450135e-9 0.0 0.0 0.0;;; 3.411340562667307e-6 0.0 0.0 0.0; 0.0 0.0 0.0 0.0;;; -3.2811221779399064e-6 0.0 0.0 0.0; 1.4766336464937334e-6 0.0 0.0 0.0;;; 7.508033462390918e-7 0.0 0.0 0.0; -8.45140417586639e-7 0.0 0.0 0.0;;; -5.0408680959505405e-8 0.0 0.0 0.0; 1.5954795989679934e-7 0.0 0.0 0.0;;; -1.3643415621015169e-9 0.0 0.0 0.0; -1.2019285292447762e-8 0.0 0.0 0.0;;; 1.7490047653047897e-10 0.0 0.0 0.0; 2.9802226104568586e-10 0.0 0.0 0.0;;; 7.320829843849779e-7 0.0 0.0 0.0; 0.0 0.0 0.0 0.0;;; -1.1040063516878706e-6 0.0 0.0 0.0; 4.97169838052928e-7 0.0 0.0 0.0;;; 3.3233488404786064e-7 0.0 0.0 0.0; -3.7435741595385594e-7 0.0 0.0 0.0;;; -3.0399512277559916e-8 0.0 0.0 0.0; 9.636699936286549e-8 0.0 0.0 0.0;;; -1.257805313390412e-9 0.0 0.0 0.0; -1.1028077206194902e-8 0.0 0.0 0.0;;; 3.239977186634296e-10 0.0 0.0 0.0; 5.512947043790685e-10 0.0 0.0 0.0;;; -1.2339557354154617e-11 0.0 0.0 0.0; -8.868558661043951e-12 0.0 0.0 0.0;;; 1.0121462827844036e-7 0.0 0.0 0.0; 0.0 0.0 0.0 0.0;;; -2.9977386626355873e-7 0.0 0.0 0.0; 1.3510437997734334e-7 0.0 0.0 0.0;;; 1.1575725282596027e-7 0.0 0.0 0.0; -1.3049796958370626e-7 0.0 0.0 0.0;;; -1.3594389585018197e-8 0.0 0.0 0.0; 4.3167709924881135e-8 0.0 0.0 0.0;;; -7.667019220226146e-10 0.0 0.0 0.0; -6.6880743655575325e-9 0.0 0.0 0.0;;; 2.9803138008504807e-10 0.0 0.0 0.0; 5.0634784455077e-10 0.0 0.0 0.0;;; -2.284557373740596e-11 0.0 0.0 0.0; -1.6393720070527073e-11 0.0 0.0 0.0;;; 5.595126561276888e-13 0.0 0.0 0.0; 1.1787154792703014e-13 0.0 0.0 0.0;;; 1.4751810727352213e-9 0.0 0.0 0.0; 0.0 0.0 0.0 0.0;;; -6.708468877007811e-8 0.0 0.0 0.0; 3.026408896811465e-8 0.0 0.0 0.0;;; 3.300201956443528e-8 0.0 0.0 0.0; -3.723787837566326e-8 0.0 0.0 0.0;;; -4.805896384168386e-9 0.0 0.0 0.0; 1.5289010248988546e-8 0.0 0.0 0.0;;; -3.4754378686830237e-10 0.0 0.0 0.0; -3.0152402096102556e-9 0.0 0.0 0.0;;; 1.81479446950944e-10 0.0 0.0 0.0; 3.07835507484142e-10 0.0 0.0 0.0;;; -2.1022722355702963e-11 0.0 0.0 0.0; -1.506079745721181e-11 0.0 0.0 0.0;;; 1.035467726604547e-12 0.0 0.0 0.0; 2.1721150343128253e-13 0.0 0.0 0.0;;; -1.8395652399672212e-14 0.0 0.0 0.0; 3.852107985386056e-15 0.0 0.0 0.0;;; -4.0609549751469545e-9 0.0 0.0 0.0; 0.0 0.0 0.0 0.0;;; -1.2482117506512014e-8 0.0 0.0 0.0; 5.638433267215509e-9 0.0 0.0 0.0;;; 7.905670401846284e-9 0.0 0.0 0.0; -8.929464291174913e-9 0.0 0.0 0.0;;; -1.3981037159203388e-9 0.0 0.0 0.0; 4.456832611637532e-9 0.0 0.0 0.0;;; -1.2492898592786975e-10 0.0 0.0 0.0; -1.0775898987149912e-9 0.0 0.0 0.0;;; 8.228418174866392e-11 0.0 0.0 0.0; 1.3933781542187157e-10 0.0 0.0 0.0;;; -1.2818959531322538e-11 0.0 0.0 0.0; -9.167593119553919e-12 0.0 0.0 0.0;;; 9.531511322909608e-13 0.0 0.0 0.0; 1.9904649104822634e-13 0.0 0.0 0.0;;; -3.402893923417914e-14 0.0 0.0 0.0; 7.160468249308245e-15 0.0 0.0 0.0;;; 4.465711619580417e-16 0.0 0.0 0.0; -3.1804698481402137e-16 0.0 0.0 0.0;;; -1.4824011156256148e-9 0.0 0.0 0.0; 0.0 0.0 0.0 0.0;;; -1.919833813101314e-9 0.0 0.0 0.0; 8.688650303478487e-10 0.0 0.0 0.0;;; 1.6200773449979007e-9 0.0 0.0 0.0; -1.8320295085795508e-9 0.0 0.0 0.0;;; -3.4396303878189783e-10 0.0 0.0 0.0; 1.0989242195318225e-9 0.0 0.0 0.0;;; -3.708192715671552e-11 0.0 0.0 0.0; -3.1787656903648365e-10 0.0 0.0 0.0;;; 2.962574958499155e-11 0.0 0.0 0.0; 5.007695694535446e-11 0.0 0.0 0.0;;; -5.8262901006951636e-12 0.0 0.0 0.0; -4.159070625703628e-12 0.0 0.0 0.0;;; 5.818200312063137e-13 0.0 0.0 0.0; 1.2092681969926203e-13 0.0 0.0 0.0;;; -3.1327310064348024e-14 0.0 0.0 0.0; 6.625396053471999e-15 0.0 0.0 0.0;;; 8.255548969437712e-16 0.0 0.0 0.0; -5.893192421106608e-16 0.0 0.0 0.0;;; -7.403943579266824e-18 0.0 0.0 0.0; 1.2317599122655843e-17 0.0 0.0 0.0;;; -3.515934968301032e-10 0.0 0.0 0.0; 0.0 0.0 0.0 0.0;;; -2.3594925434518903e-10 0.0 0.0 0.0; 1.07127957503847e-10 0.0 0.0 0.0;;; 2.8743136385190034e-10 0.0 0.0 0.0; -3.2548499126874616e-10 0.0 0.0 0.0;;; -7.297407969269781e-11 0.0 0.0 0.0; 2.3371958216102407e-10 0.0 0.0 0.0;;; -9.344496321234191e-12 0.0 0.0 0.0; -7.957394809536844e-11 0.0 0.0 0.0;;; 8.820872426341902e-12 0.0 0.0 0.0; 1.4881530759394605e-11 0.0 0.0 0.0;;; -2.10510183229474e-12 0.0 0.0 0.0; -1.4998051808823597e-12 0.0 0.0 0.0;;; 2.6492865541168975e-13 0.0 0.0 0.0; 5.4789022475871284e-14 0.0 0.0 0.0;;; -1.9135933669548904e-14 0.0 0.0 0.0; 4.0683844309145415e-15 0.0 0.0 0.0;;; 7.598586435383344e-16 0.0 0.0 0.0; -5.43729257580773e-16 0.0 0.0 0.0;;; -1.366948657228694e-17 0.0 0.0 0.0; 2.2803585253454954e-17 0.0 0.0 0.0;;; 4.513564991747822e-20 0.0 0.0 0.0; -3.389937360059917e-19 0.0 0.0 0.0;;; -6.573572046881621e-11 0.0 0.0 0.0; 0.0 0.0 0.0 0.0;;; -2.043863801299287e-11 0.0 0.0 0.0; 9.352606219154242e-12 0.0 0.0 0.0;;; 4.445562910300353e-11 0.0 0.0 0.0; -5.042557742977247e-11 0.0 0.0 0.0;;; -1.3543695917757064e-11 0.0 0.0 0.0; 4.3496602156693426e-11 0.0 0.0 0.0;;; -2.0397436768532417e-12 0.0 0.0 0.0; -1.7246583947075334e-11 0.0 0.0 0.0;;; 2.2333444964799568e-12 0.0 0.0 0.0; 3.76016740753104e-12 0.0 0.0 0.0;;; -6.297270638036358e-13 0.0 0.0 0.0; -4.477390964671143e-13 0.0 0.0 0.0;;; 9.597644376139124e-14 0.0 0.0 0.0; 1.974437557076672e-14 0.0 0.0 0.0;;; -8.724659175859616e-15 0.0 0.0 0.0; 1.8650669605364297e-15 0.0 0.0 0.0;;; 4.642657322793347e-16 0.0 0.0 0.0; -3.3304590676846863e-16 0.0 0.0 0.0;;; -1.2569421365745768e-17 0.0 0.0 0.0; 2.1028191987394226e-17 0.0 0.0 0.0;;; 8.266548146039412e-20 0.0 0.0 0.0; -6.271851481459042e-19 0.0 0.0 0.0;;; 2.019255448918255e-21 0.0 0.0 0.0; 7.17077959778311e-21 0.0 0.0 0.0;;; -1.0278419931152499e-11 0.0 0.0 0.0; 0.0 0.0 0.0 0.0;;; -4.203242166211679e-13 0.0 0.0 0.0; 2.1066117376944255e-13 0.0 0.0 0.0;;; 6.002693212564717e-12 0.0 0.0 0.0; -6.823256754648399e-12 0.0 0.0 0.0;;; -2.222360876797338e-12 0.0 0.0 0.0; 7.159343863048163e-12 0.0 0.0 0.0;;; -3.915530871566449e-13 0.0 0.0 0.0; -3.2854874061920314e-12 0.0 0.0 0.0;;; 4.906975535591016e-13 0.0 0.0 0.0; 8.243739866537346e-13 0.0 0.0 0.0;;; -1.603914412140437e-13 0.0 0.0 0.0; -1.1379290070808765e-13 0.0 0.0 0.0;;; 2.881187837425885e-14 0.0 0.0 0.0; 5.8944206681181154e-15 0.0 0.0 0.0;;; -3.166737553935407e-15 0.0 0.0 0.0; 6.808073260383793e-16 0.0 0.0 0.0;;; 2.11823342218566e-16 0.0 0.0 0.0; -1.5235001948604837e-16 0.0 0.0 0.0;;; -7.674869280589808e-18 0.0 0.0 0.0; 1.2877827358552214e-17 0.0 0.0 0.0;;; 7.540003064291612e-20 0.0 0.0 0.0; -5.781607461451349e-19 0.0 0.0 0.0;;; 3.753862216640856e-21 0.0 0.0 0.0; 1.3258129373615239e-20 0.0 0.0 0.0;;; -9.398251768197244e-23 0.0 0.0 0.0; -1.164083006095856e-22 0.0 0.0 0.0;;; -1.3762401781072233e-12 0.0 0.0 0.0; 0.0 0.0 0.0 0.0;;; 2.8087332169094867e-13 0.0 0.0 0.0; -1.2390042397442847e-13 0.0 0.0 0.0;;; 7.032457492620986e-13 0.0 0.0 0.0; -8.016843386759864e-13 0.0 0.0 0.0;;; -3.248637046058368e-13 0.0 0.0 0.0; 1.0502423099069615e-12 0.0 0.0 0.0;;; -6.687546341485134e-14 0.0 0.0 0.0; -5.565382595152538e-13 0.0 0.0 0.0;;; 9.500740818317414e-14 0.0 0.0 0.0; 1.5924469307571103e-13 0.0 0.0 0.0;;; -3.5499221407458265e-14 0.0 0.0 0.0; -2.5128180166288814e-14 0.0 0.0 0.0;;; 7.370959401207499e-15 0.0 0.0 0.0; 1.4991782388637415e-15 0.0 0.0 0.0;;; -9.53078284611239e-16 0.0 0.0 0.0; 2.0611156438331795e-16 0.0 0.0 0.0;;; 7.697575427474568e-17 0.0 0.0 0.0; -5.551354761137937e-17 0.0 0.0 0.0;;; -3.500656875141304e-18 0.0 0.0 0.0; 5.891922242973979e-18 0.0 0.0 0.0;;; 4.5664032899842155e-20 0.0 0.0 0.0; -3.540569052057126e-19 0.0 0.0 0.0;;; 3.4784857373427885e-21 0.0 0.0 0.0; 1.2216677736179948e-20 0.0 0.0 0.0;;; -1.7414905033630244e-22 0.0 0.0 0.0; -2.1502849742163935e-22 0.0 0.0 0.0;;; 2.460117981125405e-24 0.0 0.0 0.0; 1.3094830233058996e-24 0.0 0.0 0.0]
expansion_order = Val{14}()
fmm._B2M!_panel(multipole_expansion, qnm_prev, jnm_prev, inm_prev, R0, Ru, Rv, strength, normal, expansion_order, fmm.UniformSourcePanel())
for i in eachindex(multipole_expansion)
    @test isapprox(multipole_expansion[i], multipole_expansion_test[i]; atol=1e-11)
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
strength = 1.0
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
panel_type = fmm.UniformNormalDipolePanel()
fmm._B2M!_panel(branch.multipole_expansion, qnm_prev, jnm_prev, inm_prev, R0, Ru, Rv, strength, normal, expansion_order, panel_type)

# evaluate
target_potential_dipole = zeros(4)
irregular_harmonics_dipole = FastMultipole.initialize_harmonics(P,Float64)
multipole_expansion_dipole = branch.multipole_expansion
FastMultipole.M2B!(target_potential_dipole, target, displ, irregular_harmonics_dipole, multipole_expansion_dipole, P)

this_potential = 0.0015577438029241901

@test isapprox(target_potential_dipole[1], this_potential; atol=1e-12)

# reset and repeat for source panel
branch.multipole_expansion .= 0.0

panel_type = fmm.UniformSourcePanel()
fmm._B2M!_panel(branch.multipole_expansion, qnm_prev, jnm_prev, inm_prev, R0, Ru, Rv, strength, normal, expansion_order, panel_type)

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
probes = fmm.ProbeSystem(bodies[1:3,2:end]; scalar_potential=true, vector_potential=true, velocity=true, velocity_gradient=true)
mass = Gravitational(bodies)
fmm.fmm!(probes, mass; expansion_order=5, n_per_branch_source=50, n_per_branch_target=50, multipole_acceptance_criterion=0.4)
fmm.fmm!(mass; expansion_order=5, n_per_branch=50, multipole_acceptance_criterion=0.4)

for i_mass in 2:n_bodies
    @test isapprox(mass.potential[1,i_mass], probes.scalar_potential[i_mass-1]; atol=1e-12)
end

end