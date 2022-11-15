using Test
import Statistics
S = Statistics
import Symbolics
sym = Symbolics

import FLOWFMM
fmm = FLOWFMM

using LegendrePolynomials
import PyPlot
plt = PyPlot
using LaTeXStrings

test_dir = @__DIR__

#####
##### define gravitational kernel and mass elements
#####
include(joinpath(test_dir, "gravitational.jl"))

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

    bodies = vcat(x,m',zeros(3,length(m)))
    mass = Gravitational(bodies)

    mass.direct!(mass)
    V_tots_direct = mass.potential[1,:]

    for i in 1:length(V_tots)
        @test isapprox(V_tots[i], V_tots_direct[i]; atol=1e-4)
    end
end

# #####
# ##### define plotting functions
# #####
# function plot_leaf(elements, leaf, fig_name, save_name;
#     stl = "*", clr = "",
#     initialize_fig = false,
#     save_fig = false
# )
#     fig = plt.figure(fig_name)
#     if initialize_fig
#         fig.clear()
#         fig.add_subplot(111)
#     end
#     ax = fig.get_axes()[1]
#     n = length(leaf.children)
#     Xs = zeros(dims,n)
#     for (i,i_index) in enumerate(leaf.children)
#         Xs[:,i] .= fmm.get_X(mass,i_index)
#     end
#     if clr == ""
#         ax.scatter(Xs[1,:], Xs[2,:], marker=stl)
#     else
#         ax.scatter(Xs[1,:], Xs[2,:], marker=stl, color=clr)
#     end


#     if save_fig
#         fig.savefig(save_name)
#     end
# end

# function plot_leaves(elements, leaves, fig_name, save_name;
#     stls = ["v", "^", ">", "<", "+", "x", "*"],
#     clrs = ["b", "r", "g", "y", "c", "m"],
#     initialize_fig = false,
#     save_fig = false
# )

#     for (i,leaf) in enumerate(leaves)
#         initialize_fig = initialize_fig && i==1 ? true : false
#         stl = stls[i % length(stls) + 1]
#         clr = clrs[i % length(clrs) + 1]
#         this_save_fig = i == length(leaves) ? save_fig : false
#         plot_leaf(elements, leaf, fig_name, save_name;
#             stl, clr, initialize_fig, save_fig=this_save_fig
#         )
#     end
# end

# function plot_branch(elements, root, level, branch_i, fig_name, save_name;
#     stl = "+", clr = "b",
#     initialize_fig = false,
#     save_fig = false
# )
#     branch = root.branches[level][branch_i]
#     if initialize_fig
#         fig = plt.figure(fig_name)
#         fig.clear()
#         fig.add_subplot(111)
#     end

#     if level == 2
#         leaves = root.branches[1][branch.children]
#         plot_leaves(elements, leaves, fig_name, save_name;
#             stls = [stl], clrs = [clr],
#             initialize_fig = false,
#             save_fig
#         )
#     elseif level > 2
#         for branch_i in branch.children
#             plot_branch(elements, root, level-1, branch_i, fig_name, save_name;
#                 stl, clr,
#                 initialize_fig = false,
#                 save_fig
#             )
#         end
#     else
#         @error "requested plot_branch on level $level; must be >= 2"
#     end
# end

# function plot_branches(elements, root, level, fig_name, save_name;
#     stls = ["v", "^", ">", "<", "+", "x", "*"],
#     clrs = ["b", "r", "g", "y", "c", "m"],
#     initialize_fig = false,
#     save_fig = false
# )
#     if level > 1
#         branches = root.branches[level]
#         for branch_i in 1:length(branches)
#             stl = stls[branch_i % length(stls) + 1]
#             clr = clrs[branch_i % length(clrs) + 1]
#             this_initialize_fig = branch_i == 1 ? initialize_fig : false
#             plot_branch(elements, root, level, branch_i, fig_name, save_name;
#                 stl, clr,
#                 initialize_fig = this_initialize_fig,
#                 save_fig = save_fig
#             )
#         end
#     else
#         plot_leaves(elements, root.branches[level], fig_name, save_name;
#             stls,
#             clrs,
#             initialize_fig,
#             save_fig
#         )
#     end
# end

@testset "tree" begin

    # build list of elements to sort
    xs = [
        1.2 1.1 0.8;
        0.8 0.9 0.2;
        0.1 0.2 0.9;
        0.1 0.3 0.2;
        0.2 0.25 0.4
    ]
    ms = rand(size(xs)[1])
    bodies = vcat(xs',ms',zeros(3,length(ms)))
    elements = Gravitational(bodies)

    # test center_radius function
    center, radius = fmm.center_radius(elements; scale_radius = 1.00001)
    test_center = [0.65, 0.65, 0.55]
    test_radius = 0.5500055

    for i in 1:3
        @test isapprox(center[i], test_center[i]; atol=1e-4)
    end
    @test isapprox(radius, test_radius; atol=1e-4)

    # test branch! function
    tree = fmm.Tree(elements, 2, 1)

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
        @test isapprox(tree.branches[i_branch].n_bodies, test_branches[i_branch,1]; atol=1e-8)
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

bodies = zeros(7,length(ms))
for i in 1:length(ms)
    bodies[1:3,i] .= xs[i,1:3]
    bodies[4,i] = ms[i]
end
elements = Gravitational(bodies)

expansion_order = 2
tree = fmm.Tree(elements, expansion_order, 1)

i_mass = 1
i_branch = 5 # use the first mass

harmonics = zeros(Complex{Float64},(expansion_order+1)^2)
elements.B2M!(tree, tree.branches[i_branch], bodies[:,new_order_index[i_mass]], 1, harmonics)

# mp_expansion = zeros(Complex{Float64}, (expansion_order+1)^2)
# this_dx = xs[i_mass,:] - tree.branches[i_branch].center
# fmm.cartesian_2_spherical!(this_dx)
# fmm.regular_harmonic!(mp_expansion, this_dx..., expansion_order)

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

bodies = zeros(7,length(ms))
for i in 1:length(ms)
    bodies[1:3,i] .= xs[i,1:3]
    bodies[4,i] = ms[i]
end
elements = Gravitational(bodies)

expansion_order = 3
tree = fmm.Tree(elements, expansion_order, 1)

i_branch = 2 # contains 4th and 5th elements
i_branch_4 = 6 # use the fourth mass
# i_branch_5 = 7 # use the fifth mass
harmonics = zeros(Complex{Float64},(expansion_order+1)^2)
# using only the 4th mass: (note it has been reordered)
elements.B2M!(tree, tree.branches[i_branch_4], bodies[:,new_order_index[4]], 1, harmonics) # evaluate multipole coefficients
# fmm.B2M!(i_branch_5, tree, elements) # evaluate multipole coefficients
fmm.M2M!(tree, i_branch) # translate coefficients to the center of branch 2

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

bodies = zeros(7,length(ms))
for i in 1:length(ms)
    bodies[1:3,i] .= xs[i,1:3]
    bodies[4,i] = ms[i]
end
elements = Gravitational(bodies)

expansion_order = 20
tree = fmm.Tree(elements, expansion_order, 1)

branch_i = 2 # contains two elements; 4 and 5
target_i = new_order_index[4]
source_i = new_order_index[1] # just needs to be farther away than the target to ensure convergence

dx_source = fmm.cartesian_2_spherical(elements.bodies[i_POSITION,source_i] - tree.branches[branch_i].center)
dx_target = fmm.cartesian_2_spherical(elements.bodies[i_POSITION,target_i] - tree.branches[branch_i].center)

local_coefficients_theta = zeros(Complex{Float64}, ((expansion_order+1)*(expansion_order+2))>>1)
local_coefficients_expanded = zeros(Complex{Float64}, (expansion_order+1)^2)
local_coefficients_expanded_theta = zeros(Complex{Float64}, (expansion_order+1)^2)
fmm.irregular_harmonic!(local_coefficients_expanded, dx_source..., expansion_order)
local_coefficients_expanded .*= ms[1]
regular_harmonics_expanded = zeros(Complex{Float64}, (expansion_order+1)^2)
fmm.regular_harmonic!(regular_harmonics_expanded, dx_target..., expansion_order)

fmm.B2L!(tree, branch_i, elements.bodies[:,source_i])

harmonics = zeros(Complex{Float64},(expansion_order+1)^2)
harmonics_theta = zeros(Complex{Float64},(expansion_order+1)^2)
harmonics_theta_2 = zeros(Complex{Float64},(expansion_order+1)^2)
workspace = zeros(3,4)
spherical_potential = zeros(i_POTENTIAL_HESSIAN[end])
fmm.L2B!(view(elements.potential,:,target_i), harmonics, harmonics_theta, harmonics_theta_2, workspace, spherical_potential, elements.bodies[:,target_i], tree, tree.branches[branch_i])

u_fmm = elements.potential[1,target_i]

dx_direct = xs[4,:] - xs[1,:]
u_check = 1 / sqrt(dx_direct' * dx_direct)
u_check *= ms[1]

u_man = real(sum(regular_harmonics_expanded' * local_coefficients_expanded)) # appears to work

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

bodies = zeros(7,length(ms))
for i in 1:length(ms)
    bodies[1:3,i] .= xs[i,1:3]
    bodies[4,i] = ms[i]
end
elements = Gravitational(bodies)

expansion_order = 20
tree = fmm.Tree(elements, expansion_order, 1)

# local coefficient at branch 2 due to mass 1
fmm.B2L!(tree, 2, elements.bodies[:,new_order_index[1]])
# local_2 = deepcopy(tree.branches[2].local_expansion)

# check L2P now:
harmonics = zeros(Complex{Float64},(expansion_order+1)^2)
harmonics_theta = zeros(Complex{Float64},(expansion_order+1)^2)
harmonics_theta_2 = zeros(Complex{Float64},(expansion_order+1)^2)
workspace = zeros(3,4)
spherical_potential = zeros(i_POTENTIAL_HESSIAN[end])
fmm.L2B!(view(elements.potential,:,new_order_index[5]), harmonics, harmonics_theta, harmonics_theta_2, workspace, spherical_potential, elements.bodies[:,new_order_index[5]], tree, tree.branches[2])
# fmm.L2B!(elements[5], tree, tree.branches[2], harmonics, harmonics_theta, harmonics_theta_2, d_potential_cartesian, d_potential, dH_potential)
u_fmm_no_x = elements.potential[1,new_order_index[5]]
elements.potential[1,new_order_index[5]] *= 0

# translate local expansion to branch 7 (mass 5)
fmm.L2L!(tree, tree.branches[2], tree.branches[7], harmonics)

local_coefficients_check = zeros(Complex{Float64}, (expansion_order+1)^2)
dx_check = fmm.cartesian_2_spherical(elements.bodies[i_POSITION,new_order_index[1]] - tree.branches[7].center)
fmm.irregular_harmonic!(local_coefficients_check, dx_check..., expansion_order)
local_coefficients_check .*= ms[1]

# evaluate local expansion at mass 5
fmm.L2B!(tree, elements, 7)
u_fmm = elements.potential[1,new_order_index[5]]

dx_direct = elements.bodies[i_POSITION,new_order_index[5]] - elements.bodies[i_POSITION,new_order_index[1]]
u_check = ms[1] / sqrt(dx_direct' * dx_direct)

regular_harmonics = zeros(Complex{Float64}, (expansion_order+1)^2)
dx_target = fmm.cartesian_2_spherical(elements.bodies[i_POSITION,new_order_index[5]] - tree.branches[7].center)
fmm.regular_harmonic!(regular_harmonics, dx_target..., expansion_order)
u_man = real(sum(regular_harmonics' * local_coefficients_check))

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

bodies = zeros(7,length(ms))
for i in 1:length(ms)
    bodies[1:3,i] .= xs[i,1:3]
    bodies[4,i] = ms[i]
end
elements = Gravitational(bodies)

expansion_order = 30
tree = fmm.Tree(elements, expansion_order, 1)

i_branch_multipole = 7 # mass 5
i_branch_local = 5 # mass 1
harmonics = zeros(Complex{Float64}, (expansion_order+1)^2)
harmonics_theta = zeros(Complex{Float64}, (expansion_order+1)^2)
harmonics_theta_2 = zeros(Complex{Float64}, (expansion_order+1)^2)
workspace = zeros(3,4)

elements.B2M!(tree, tree.branches[i_branch_multipole], elements.bodies[:,new_order_index[5]], 1, harmonics)

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

fmm.M2L!(tree, i_branch_local, i_branch_multipole)

spherical_potential = zeros(i_POTENTIAL_HESSIAN[end])
fmm.L2B!(view(elements.potential,:,new_order_index[1]), harmonics, harmonics_theta, harmonics_theta_2, workspace, spherical_potential, elements.bodies[:,new_order_index[1]], tree, tree.branches[i_branch_local])
u_fmm = elements.potential[1,new_order_index[1]]

local_exp = tree.branches[i_branch_local].local_expansion[1]

# test local
# dx_l = xs[1,:] - tree.branches[i_branch_local].center
# fmm.regular_harmonic!(harmonics, dx_l..., expansion_order)
# u_check_local = real(sum(harmonics' * local_exp))

dx_direct = elements.bodies[i_POSITION,new_order_index[1]] - elements.bodies[i_POSITION, new_order_index[5]]
u_direct = elements.bodies[i_STRENGTH[1],new_order_index[5]] / sqrt(dx_direct' * dx_direct)

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

bodies = zeros(7,length(ms))
for i in 1:length(ms)
    bodies[1:3,i] .= xs[i,1:3]
    bodies[4,i] = ms[i]
end
elements = Gravitational(bodies)

expansion_order = 24
theta = 4
tree = fmm.Tree(elements, expansion_order, 1)

# perform upward pass
fmm.upward_pass!(tree, elements)

# m6 = tree.branches[6].multipole_expansion
target = [4.1,2.2,3.4]
dx_direct_6 = target - elements.bodies[i_POSITION,new_order_index[4]]
u_direct_6 = ms[4] / sqrt(dx_direct_6' * dx_direct_6)

mass_target_potential = zeros(4)
mass_target = target
fmm.M2B!(mass_target_potential, mass_target, 6, tree)
u_fmm_6 = mass_target_potential[1]

# add branches 6 and 7
dx_direct_7 = target - elements.bodies[i_POSITION,new_order_index[5]]
u_direct_67 = u_direct_6 + ms[5] / sqrt(dx_direct_7' * dx_direct_7)

# reset target potential
mass_target_potential *= 0

# use summed multipole expansion from branches 6 and 7 (summed at 2)
fmm.M2B!(mass_target_potential, mass_target, 2, tree)
u_fmm_67 = mass_target_potential[1]

# perform horizontal pass
fmm.horizontal_pass!(tree, elements, theta)

# consider the effect on branch 3 (mass 2)
elements.potential[i_POTENTIAL,new_order_index[2]] .*= 0 # reset potential at mass 2
# P2P is performed from branches 3 (mass 2), 4 (mass 3), and 5 (mass 1) to branch 3
# elements.direct!(elements.potential[i_POTENTIAL,new_order_index[2]], elements.bodies[i_POSITION,new_order_index[2]], elements.bodies[:,new_order_index[1]])
# elements.direct!(elements.potential[i_POTENTIAL,new_order_index[2]], elements.bodies[i_POSITION,new_order_index[2]], elements.bodies[:,new_order_index[2]])
# elements.direct!(elements.potential[i_POTENTIAL,new_order_index[2]], elements.bodies[i_POSITION,new_order_index[2]], elements.bodies[:,new_order_index[3]])
fmm.P2P!(tree, elements, 3, 3)
fmm.P2P!(tree, elements, 3, 4)
fmm.P2P!(tree, elements, 3, 5)
u_fmm_123 = elements.potential[i_POTENTIAL[1],new_order_index[2]]

dx_12 = elements.bodies[i_POSITION,new_order_index[2]] - elements.bodies[i_POSITION,new_order_index[1]]
u_direct_12 = elements.bodies[i_STRENGTH[1],new_order_index[1]] / sqrt(dx_12' * dx_12)
u_direct_22 = 0.0
dx_32 = elements.bodies[i_POSITION,new_order_index[2]] - elements.bodies[i_POSITION,new_order_index[3]]
u_direct_32 = elements.bodies[i_STRENGTH[1],new_order_index[3]] / sqrt(dx_32' * dx_32)

u_direct_123 = u_direct_12 + u_direct_22 + u_direct_32

# M2L is performed from branches 6, 7 to branch 3 (containing mass 2)
# fmm.L2B!(element, tree, branch, harmonics, harmonics_theta)
fmm.L2B!(tree, elements, 3)
u_fmm_12345 = elements.potential[i_POTENTIAL[1],new_order_index[2]]

dx_42 = elements.bodies[i_POSITION,new_order_index[4]] - elements.bodies[i_POSITION,new_order_index[2]]
u_direct_42 = elements.bodies[i_STRENGTH[1],new_order_index[4]] / sqrt(dx_42' * dx_42)
dx_52 = elements.bodies[i_POSITION,new_order_index[5]] - elements.bodies[i_POSITION,new_order_index[2]]
u_direct_52 = elements.bodies[i_STRENGTH[1],new_order_index[5]] / sqrt(dx_52' * dx_52)

u_direct_12345 = u_direct_123 + u_direct_42 + u_direct_52

@test isapprox(u_direct_123, u_fmm_123; atol=1e-12)
@test isapprox(u_direct_12345, u_fmm_12345; atol=1e-12)

# reset potentials
elements.potential .*= 0

# run fmm (reset potentials with reset_tree flag)
fmm.fmm!(tree, elements, theta; reset_tree=true)
u_fmm = deepcopy(elements.potential[1,:])

elements.potential .*= 0

fmm.direct!(elements)

u_direct = deepcopy(elements.potential[1,:])

for i in 1:size(elements.bodies)[2]
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
    @test isapprox(cartesian_hessian[i], cartesian_hessian_fd[i]; rtol=1e-5)
end

# now test FMM function
potential_hessian = zeros(3,3,4)
potential_jacobian = zeros(3,4)
for i in 1:4
    potential_hessian[:,:,i] .= spherical_hessian
    potential_jacobian[:,i] .= spherical_grad
end
workspace = zeros(3,4)

fmm.spherical_2_cartesian!(potential_jacobian, potential_hessian, workspace, r0, theta0, phi0)

for i in 1:3
    for j in 1:3
        @test isapprox(cartesian_hessian[i,j], potential_hessian[i,j,1]; rtol=1e-15)
    end
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
fmm.direct!(vortexparticles)
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

# manually build tree for testing
# Branch(n_branches, n_bodies, i_child, i_start, center, radius, multipole_expansion, local_expansion)
expansion_order = 9
x_branch_1 = [0.0,0,0]
branch_1 = fmm.Branch(2, 2, 2, 1, x_branch_1, 1/8, fmm.initialize_expansion(expansion_order), fmm.initialize_expansion(expansion_order))
x_branch_2 = xs[:,1] .+ [0.01, 0.02, -0.03]
branch_2 = fmm.Branch(-1, 1, -1, 1, x_branch_2, 1/8, fmm.initialize_expansion(expansion_order), fmm.initialize_expansion(expansion_order))
x_branch_3 = xs[:,2] .+ [0.02, -0.04, 0.01]
branch_3 = fmm.Branch(-1, 1, -1, 2, x_branch_3, 1/8, fmm.initialize_expansion(expansion_order), fmm.initialize_expansion(expansion_order))

# using FMM
# tree = fmm.Tree(branches, [expansion_order], n_per_branch, B2M!, P2P!)
tree = fmm.Tree([branch_1, branch_2, branch_3], [expansion_order], 1)
harmonics = zeros(Complex{Float64},(expansion_order+1)^2)
fmm.B2M!(tree, vortexparticles, 2)
fmm.B2M!(tree, vortexparticles, 3)

fmm.M2L!(tree, 2, 3)
fmm.M2L!(tree, 3, 2)
fmm.L2B!(tree, vortexparticles, 2)
fmm.L2B!(tree, vortexparticles, 3)
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

@testset "two 3D vortex particles" begin

bodies = [
    0.4 0.1
    0.1 -0.5
    -0.3 0.2
    0.3 -0.4
    -0.1 -0.2
    0.08 0.5
]

vortex_particles = VortexParticles(bodies)

#####
##### obtain psi, u, and stretching analytically
#####
function psi(target_x, source_x, source_gamma)
    dx = target_x - source_x
    dx_norm = sqrt(dx' * dx)
    return source_gamma ./ dx_norm
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
    hessian[:,:,1] = d2dr2 * source_gamma[1]
    hessian[:,:,2] = d2dr2 * source_gamma[2]
    hessian[:,:,3] = d2dr2 * source_gamma[3]
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

psis = zeros(3,2)
psis[:,1] = psi(bodies[1:3,1], bodies[1:3,2], bodies[4:6,2])
psis[:,2] = psi(bodies[1:3,2], bodies[1:3,1], bodies[4:6,1])
hessians = zeros(3,3,3,2)
hessians[:,:,:,1] = d2psidx2(bodies[1:3,1], bodies[1:3,2], bodies[4:6,2])
hessians[:,:,:,2] = d2psidx2(bodies[1:3,2], bodies[1:3,1], bodies[4:6,1])
us = zeros(3,2)
us[:,1] = u(bodies[1:3,1], bodies[1:3,2], bodies[4:6,2])
us[:,2] = u(bodies[1:3,2], bodies[1:3,1], bodies[4:6,1])
ss = zeros(3,2)
ss[:,1] = stretching(bodies[1:3,1], bodies[1:3,2], bodies[4:6,1], bodies[4:6,2])
ss[:,2] = stretching(bodies[1:3,2], bodies[1:3,1], bodies[4:6,2], bodies[4:6,1])

#####
##### use direct method
#####
fmm.direct!(vortex_particles)
update_velocity_stretching!(vortex_particles)

psis_direct = deepcopy(vortex_particles.potential[2:4,:])
for i in 1:length(psis_direct)
    @test isapprox(psis_direct[i], psis[i]; atol=1e-10)
end
hessians_direct = deepcopy(reshape(vortex_particles.potential[i_POTENTIAL_HESSIAN[10:end],:],3,3,3,2))
for i in 1:length(hessians)
    @test isapprox(hessians_direct[i], hessians[i]; atol=1e-10)
end
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
x_branch_1 = (bodies[1:3,1] + bodies[1:3,2])/2
branch_1 = fmm.Branch(2, 2, 2, 1, x_branch_1, 1/8, fmm.initialize_expansion(expansion_order), fmm.initialize_expansion(expansion_order))
x_branch_2 = bodies[1:3,1] .+ [0.01, 0.02, -0.03]
branch_2 = fmm.Branch(-1, 1, -1, 1, x_branch_2, 1/8, fmm.initialize_expansion(expansion_order), fmm.initialize_expansion(expansion_order))
x_branch_3 = bodies[1:3,2] .+ [0.02, -0.04, 0.01]
branch_3 = fmm.Branch(-1, 1, -1, 2, x_branch_3, 1/8, fmm.initialize_expansion(expansion_order), fmm.initialize_expansion(expansion_order))

# tree = fmm.Tree(branches, [expansion_order], n_per_branch)
tree = fmm.Tree([branch_1, branch_2, branch_3], [expansion_order], 1)
# fmm.B2M!(tree, vortex_particles, 2)
# fmm.B2M!(tree, vortex_particles, 3)

# fmm.M2L!(tree, 2, 3)
# fmm.M2L!(tree, 3, 2)
# fmm.L2B!(tree, vortex_particles, 2)
# fmm.L2B!(tree, vortex_particles, 3)
theta = 4
fmm.fmm!(tree, vortex_particles, theta)
update_velocity_stretching!(vortex_particles)

psis_fmm = deepcopy(vortex_particles.potential[2:4,:])
for i in 1:length(psis_fmm)
    @test isapprox(psis_fmm[i], psis[i]; atol=1e-10)
end
hessians_fmm = deepcopy(reshape(vortex_particles.potential[i_POTENTIAL_HESSIAN[10:end],:],3,3,3,2))
for i in 1:length(hessians)
    @test isapprox(hessians_fmm[i], hessians[i]; atol=1e-8)
end
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
    0.3 -0.4 0.2
    -0.1 -0.2 0.5
    0.08 0.5 1.1
]

vortex_particles = VortexParticles(bodies)

#####
##### obtain psi, u, and stretching analytically
#####
function psi(target_x, source_x, source_gamma)
    dx = target_x - source_x
    dx_norm = sqrt(dx' * dx)
    return source_gamma ./ dx_norm
end

function dpsidx(target_x, source_x, source_gamma)
    dx = target_x - source_x
    dx_norm = sqrt(dx' * dx)
    x, y, z = dx
    jacobian = [
        -x*source_gamma[1] -x*source_gamma[2] -x*source_gamma[3];
        -y*source_gamma[1] -y*source_gamma[2] -y*source_gamma[3];
        -z*source_gamma[1] -z*source_gamma[2] -z*source_gamma[3];
    ] ./ dx_norm^3
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
    hessian[:,:,1] = d2dr2 * source_gamma[1]
    hessian[:,:,2] = d2dr2 * source_gamma[2]
    hessian[:,:,3] = d2dr2 * source_gamma[3]
    return hessian
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
    hessian[:,:,1] = d2dr2 * source_gamma[1]
    hessian[:,:,2] = d2dr2 * source_gamma[2]
    hessian[:,:,3] = d2dr2 * source_gamma[3]
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

psis = zeros(3,3)
psis[:,1] = psi(bodies[1:3,1], bodies[1:3,2], bodies[4:6,2]) + psi(bodies[1:3,1], bodies[1:3,3], bodies[4:6,3])
psis[:,2] = psi(bodies[1:3,2], bodies[1:3,1], bodies[4:6,1]) + psi(bodies[1:3,2], bodies[1:3,3], bodies[4:6,3])
psis[:,3] = psi(bodies[1:3,3], bodies[1:3,1], bodies[4:6,1]) + psi(bodies[1:3,3], bodies[1:3,2], bodies[4:6,2])
jacobians = zeros(3,3,3)
jacobians[:,:,1] = dpsidx(bodies[1:3,1], bodies[1:3,2], bodies[4:6,2]) + dpsidx(bodies[1:3,1], bodies[1:3,3], bodies[4:6,3])
jacobians[:,:,2] = dpsidx(bodies[1:3,2], bodies[1:3,1], bodies[4:6,1]) + dpsidx(bodies[1:3,2], bodies[1:3,3], bodies[4:6,3])
jacobians[:,:,3] = dpsidx(bodies[1:3,3], bodies[1:3,1], bodies[4:6,1]) + dpsidx(bodies[1:3,3], bodies[1:3,2], bodies[4:6,2])
hessians = zeros(3,3,3,3)
hessians[:,:,:,1] = d2psidx2(bodies[1:3,1], bodies[1:3,2], bodies[4:6,2]) + d2psidx2(bodies[1:3,1], bodies[1:3,3], bodies[4:6,3])
hessians[:,:,:,2] = d2psidx2(bodies[1:3,2], bodies[1:3,1], bodies[4:6,1]) + d2psidx2(bodies[1:3,2], bodies[1:3,3], bodies[4:6,3])
hessians[:,:,:,3] = d2psidx2(bodies[1:3,3], bodies[1:3,1], bodies[4:6,1]) + d2psidx2(bodies[1:3,3], bodies[1:3,2], bodies[4:6,2])
us = zeros(3,3)
us[:,1] = u(bodies[1:3,1], bodies[1:3,2], bodies[4:6,2]) + u(bodies[1:3,1], bodies[1:3,3], bodies[4:6,3])
us[:,2] = u(bodies[1:3,2], bodies[1:3,1], bodies[4:6,1]) + u(bodies[1:3,2], bodies[1:3,3], bodies[4:6,3])
us[:,3] = u(bodies[1:3,3], bodies[1:3,1], bodies[4:6,1]) + u(bodies[1:3,3], bodies[1:3,2], bodies[4:6,2])
ss = zeros(3,3)
ss[:,1] = stretching(bodies[1:3,1], bodies[1:3,2], bodies[4:6,1], bodies[4:6,2]) + stretching(bodies[1:3,1], bodies[1:3,3], bodies[4:6,1], bodies[4:6,3])
ss[:,2] = stretching(bodies[1:3,2], bodies[1:3,1], bodies[4:6,2], bodies[4:6,1]) + stretching(bodies[1:3,2], bodies[1:3,3], bodies[4:6,2], bodies[4:6,3])
ss[:,3] = stretching(bodies[1:3,3], bodies[1:3,1], bodies[4:6,3], bodies[4:6,1]) + stretching(bodies[1:3,3], bodies[1:3,2], bodies[4:6,3], bodies[4:6,2])

#####
##### use direct method
#####
fmm.direct!(vortex_particles)
update_velocity_stretching!(vortex_particles)

psis_direct = deepcopy(vortex_particles.potential[2:4,:])
for i in 1:length(psis_direct)
    @test isapprox(psis_direct[i], psis[i]; atol=1e-10)
end
jacobians_direct = deepcopy(reshape(vortex_particles.potential[i_POTENTIAL_JACOBIAN[4:end],:],3,3,3))
for i in 1:length(jacobians_direct)
    @test isapprox(jacobians_direct[i], jacobians[i]; atol=1e-10)
end
hessians_direct = deepcopy(reshape(vortex_particles.potential[i_POTENTIAL_HESSIAN[10:end],:],3,3,3,3))
for i in 1:length(hessians)
    @test isapprox(hessians_direct[i], hessians[i]; atol=1e-10)
end
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

# branch = Branch(n_branches, n_bodies, i_child, i_start, center, radius, multipole_expansion, local_expansion)
expansion_order = 32
n_per_branch = 1
# x_branch_1 = (bodies[1:3,1] + bodies[1:3,2] + bodies[1:3,3])/3
# branch_1 = fmm.Branch(3, 3, 2, 1, x_branch_1, 1/8, fmm.initialize_expansion(expansion_order), fmm.initialize_expansion(expansion_order))
# x_branch_2 = bodies[1:3,1] .+ [0.01, 0.02, -0.03]
# branch_2 = fmm.Branch(-1, 1, -1, 1, x_branch_2, 1/8, fmm.initialize_expansion(expansion_order), fmm.initialize_expansion(expansion_order))
# x_branch_3 = bodies[1:3,2] .+ [0.02, -0.04, 0.01]
# branch_3 = fmm.Branch(-1, 1, -1, 2, x_branch_3, 1/8, fmm.initialize_expansion(expansion_order), fmm.initialize_expansion(expansion_order))
# x_branch_4 = bodies[1:3,3] .+ [-0.01, 0.01, -0.04]
# branch_4 = fmm.Branch(-1, 1, -1, 3, x_branch_4, 1/8, fmm.initialize_expansion(expansion_order), fmm.initialize_expansion(expansion_order))

# tree = fmm.Tree(branches, [expansion_order], n_per_branch)
# tree = fmm.Tree([branch_1, branch_2, branch_3, branch_4], [expansion_order], 1)
tree = fmm.Tree(vortex_particles, expansion_order, n_per_branch)
# fmm.B2M!(tree, vortex_particles, 2)
# fmm.B2M!(tree, vortex_particles, 3)

# fmm.M2L!(tree, 2, 3)
# fmm.M2L!(tree, 3, 2)
# fmm.L2B!(tree, vortex_particles, 2)
# fmm.L2B!(tree, vortex_particles, 3)
theta =2.0
fmm.fmm!(tree, vortex_particles, theta)

update_velocity_stretching!(vortex_particles)

psis_fmm = deepcopy(vortex_particles.potential[2:4,:])
for i in 1:length(psis_fmm)
    @test isapprox(psis_fmm[i], psis[i]; rtol=1e-12)
end
hessians_fmm = deepcopy(reshape(vortex_particles.potential[i_POTENTIAL_HESSIAN[10:end],:],3,3,3,3))
for i in 1:length(hessians)
    @test isapprox(hessians_fmm[i], hessians[i]; rtol=1e-12)
end
us_fmm = deepcopy(vortex_particles.velocity_stretching[1:3,:])
for i in 1:length(us)
    @test isapprox(us_fmm, us;rtol=1e-12)
end
ss_fmm = deepcopy(vortex_particles.velocity_stretching[4:6,:])
for i in 1:length(ss)
    @test isapprox(ss_fmm[i], ss[i];rtol=1e-12)
end

end
