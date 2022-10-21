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

    mass = [Mass(x[:,i], [m[i],0,0,0], zeros(4), zeros(3)) for i in 1:length(m)]

    fmm.direct!(mass)
    V_tots_direct = [m.potential[1] for m in mass]

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
    masses = Vector{Mass}(undef,length(ms))
    for i in 1:length(ms)
        x = xs[i,:]
        mass = rand(1)
        potential = zeros(1)
        force = zeros(3)
        masses[i] = Mass(x,mass,potential,force)
    end

    # test center_radius function
    center, radius = fmm.center_radius(masses; scale_radius = 1.00001)
    test_center = [0.65, 0.65, 0.55]
    test_radius = 0.5500055

    for i in 1:3
        @test isapprox(center[i], test_center[i]; atol=1e-4)
    end
    @test isapprox(radius, test_radius; atol=1e-4)

    # test branch! function
    tree = fmm.Tree(masses, derivatives)

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
        @test isapprox(tree.branches[i_branch].n_elements, test_branches[i_branch,1]; atol=1e-8)
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

masses = Vector{Mass}(undef,length(ms))
for i in 1:length(ms)
    x = xs[i,:]
    mass = [ms[i], 0, 0, 0]
    potential = zeros(4)
    force = zeros(3)
    masses[i] = Mass(x,mass,potential,force)
end

expansion_order = 2
tree = fmm.Tree(masses, P2M!; expansion_order)

i_mass = 1
i_branch = 5 # use the first mass
fmm.P2M!(tree, masses, i_branch)
center = tree.branches[i_branch].center

x_target = [10.1,-7.3,8.6]
target = Mass(x_target, [1.0, 0,0,0], [0.0, 0,0,0], zeros(3))
fmm.M2P!(target, i_branch, tree)

u_fmm = target.potential[1]

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

masses = Vector{Mass}(undef,length(ms))
for i in 1:length(ms)
    x = xs[i,:]
    mass = [ms[i], 0, 0, 0]
    potential = zeros(4)
    force = zeros(3)
    masses[i] = Mass(x,mass,potential,force)
end

expansion_order = 3
tree = fmm.Tree(masses, P2M!; expansion_order)

i_branch = 2 # contains 4th and 5th masses
i_branch_4 = 6 # use the fourth mass
# i_branch_5 = 7 # use the fifth mass
fmm.P2M!(tree, masses, i_branch_4) # evaluate multipole coefficients
# fmm.P2M!(i_branch_5, tree, masses) # evaluate multipole coefficients
fmm.M2M!(tree, i_branch) # translate coefficients to the center of branch 2

x_target = [8.3,1.4,-4.2]
target = Mass(x_target, [0.0,0,0,0], [0.0,0,0,0], zeros(3))
fmm.M2P!(target, i_branch, tree)
u_fmm = target.potential[1]

target.potential .*= 0
fmm.M2P!(target, i_branch_4, tree)
u_fmm_no_x = target.potential[1]

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

masses = Vector{Mass}(undef,length(ms))
for i in 1:length(ms)
    x = xs[i,:]
    mass = [ms[i], 0, 0, 0]
    potential = zeros(4)
    force = zeros(3)
    masses[i] = Mass(x,mass,potential,force)
end

expansion_order = 20
tree = fmm.Tree(masses, P2M!; expansion_order)

branch_i = 2 # contains two masses; 4 and 5
target_i = 4
source_i = 1 # just needs to be farther away than the target to ensure convergence

dx_source = fmm.cartesian_2_spherical(masses[source_i].position - tree.branches[branch_i].center)
dx_target = fmm.cartesian_2_spherical(masses[target_i].position - tree.branches[branch_i].center)

local_coefficients_theta = zeros(Complex{Float64}, ((expansion_order+1)*(expansion_order+2))>>1)
local_coefficients_expanded = zeros(Complex{Float64}, (expansion_order+1)^2)
local_coefficients_expanded_theta = zeros(Complex{Float64}, (expansion_order+1)^2)
fmm.irregular_harmonic!(local_coefficients_expanded, dx_source..., expansion_order)
local_coefficients_expanded .*= ms[source_i]
regular_harmonics_expanded = zeros(Complex{Float64}, (expansion_order+1)^2)
regular_harmonics_theta_expanded = zeros(Complex{Float64}, (expansion_order+1)^2)
fmm.regular_harmonic!(regular_harmonics_expanded, regular_harmonics_theta_expanded, dx_target..., expansion_order)

fmm.P2L!(tree, branch_i, masses[source_i])

harmonics = zeros(Complex{Float64},(expansion_order+1)^2)
harmonics_theta = zeros(Complex{Float64},(expansion_order+1)^2)
d_potential = zeros(4)
fmm.L2P!(masses[target_i], tree, tree.branches[branch_i], harmonics, harmonics_theta, d_potential)

u_fmm = masses[target_i].potential[1]

dx_direct = xs[target_i,:] - xs[source_i,:]
u_check = 1 / sqrt(dx_direct' * dx_direct)
u_check *= ms[source_i]

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

masses = Vector{Mass}(undef,length(ms))
for i in 1:length(ms)
    x = xs[i,:]
    mass = [ms[i], 0, 0, 0]
    potential = zeros(4)
    force = zeros(3)
    masses[i] = Mass(x,mass,potential,force)
end

expansion_order = 20
tree = fmm.Tree(masses, P2M!; expansion_order)

# local coefficient at branch 2 due to mass 1
fmm.P2L!(tree, 2, masses[1])
# local_2 = deepcopy(tree.branches[2].local_expansion)

# check L2P now:
harmonics = zeros(Complex{Float64},(expansion_order+1)^2)
harmonics_theta = zeros(Complex{Float64},(expansion_order+1)^2)
d_potential = zeros(4)
fmm.L2P!(masses[5], tree, tree.branches[2], harmonics, harmonics_theta, d_potential)
u_fmm_no_x = masses[5].potential[1]
masses[5].potential[1] *= 0

# translate local expansion to branch 7 (mass 5)
fmm.L2L!(tree, tree.branches[2], tree.branches[7], harmonics, harmonics_theta)

local_coefficients_check = zeros(Complex{Float64}, (expansion_order+1)^2)
dx_check = fmm.cartesian_2_spherical(masses[1].position - tree.branches[7].center)
fmm.irregular_harmonic!(local_coefficients_check, dx_check..., expansion_order)
local_coefficients_check .*= ms[1]

# evaluate local expansion at mass 5
fmm.L2P!(tree, masses, 7)
u_fmm = masses[5].potential[1]

dx_direct = masses[5].position - masses[1].position
u_check = ms[1] / sqrt(dx_direct' * dx_direct)

regular_harmonics = zeros(Complex{Float64}, (expansion_order+1)^2)
regular_harmonics_theta = zeros(Complex{Float64}, (expansion_order+1)^2)
dx_target = fmm.cartesian_2_spherical(masses[5].position - tree.branches[7].center)
fmm.regular_harmonic!(regular_harmonics, regular_harmonics_theta, dx_target..., expansion_order)
u_man = real(sum(regular_harmonics' * local_coefficients_check))

@test isapprox(u_check, u_fmm; atol=1e-12)
@test isapprox(u_check, u_fmm_no_x; atol=1e-12)
@test isapprox(u_check, u_man; atol=1e-12)

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

masses = Vector{Mass}(undef,length(ms))
for i in 1:length(ms)
    x = xs[i,:]
    mass = [ms[i], 0, 0, 0]
    potential = zeros(4)
    force = zeros(3)
    masses[i] = Mass(x,mass,potential,force)
end

expansion_order = 30
tree = fmm.Tree(masses, P2M!; expansion_order)

i_branch_multipole = 7 # mass 5
i_branch_local = 5 # mass 1
harmonics = zeros(Complex{Float64}, (expansion_order+1)^2)
harmonics_theta = zeros(Complex{Float64}, (expansion_order+1)^2)

tree.P2M!(tree, tree.branches[i_branch_multipole], masses[5], harmonics)
fmm.M2L!(tree, masses, i_branch_local, i_branch_multipole)
d_potential = zeros(4)
fmm.L2P!(masses[1], tree, tree.branches[i_branch_local], harmonics, harmonics_theta, d_potential)

u_fmm = masses[1].potential[1]

dx_direct = masses[1].position - masses[5].position
u_direct = masses[5].strength[1] / sqrt(dx_direct' * dx_direct)

@test isapprox(u_fmm, u_direct; atol=1e-12)

# now check coefficients
# tree_2 = fmm.Tree(masses, P2M!; expansion_order)
# fmm.P2L!(tree_2, i_branch_local, masses[5])
# tree_2.branches[i_branch_local].local_expansion tree.branches[i_branch_local].local_expansion

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

masses = Vector{Mass}(undef,length(ms))
for i in 1:length(ms)
    x = xs[i,:]
    mass = [ms[i], 0, 0, 0]
    potential = zeros(4)
    force = zeros(3)
    masses[i] = Mass(x,mass,potential,force)
end

expansion_order = 24
theta = 4
tree = fmm.Tree(masses, P2M!; expansion_order)

# perform upward pass
fmm.upward_pass!(tree, masses)

# m6 = tree.branches[6].multipole_expansion
target = [4.1,2.2,3.4]
dx_direct_6 = target - masses[4].position
u_direct_6 = ms[4] / sqrt(dx_direct_6' * dx_direct_6)

mass_target = Mass(target, [0.0,0,0,0], [0.0,0,0,0], [0.0])

fmm.M2P!(mass_target, 6, tree)
u_fmm_6 = mass_target.potential[1]

# add branches 6 and 7
dx_direct_7 = target - masses[5].position
u_direct_67 = u_direct_6 + ms[5] / sqrt(dx_direct_7' * dx_direct_7)

# reset target potential
mass_target.potential[1] *= 0

# use summed multipole expansion from branches 6 and 7 (summed at 2)
fmm.M2P!(mass_target, 2, tree)
u_fmm_67 = mass_target.potential[1]

# perform horizontal pass
fmm.horizontal_pass!(tree, masses, theta)

# consider the effect on branch 3 (mass 2)
masses[2].potential .*= 0 # reset potential at mass 2
# P2P is performed from branches 3 (mass 2), 4 (mass 3), and 5 (mass 1) to branch 3
fmm.P2P!(tree, masses, 3, 3)
fmm.P2P!(tree, masses, 3, 4)
fmm.P2P!(tree, masses, 3, 5)
u_fmm_123 = masses[2].potential[1]

dx_12 = masses[2].position - masses[1].position
u_direct_12 = masses[1].strength[1] / sqrt(dx_12' * dx_12)
u_direct_22 = 0.0
dx_32 = masses[2].position - masses[3].position
u_direct_32 = masses[3].strength[1] / sqrt(dx_32' * dx_32)

u_direct_123 = u_direct_12 + u_direct_22 + u_direct_32

# M2L is performed from branches 6, 7 to branch 3 (containing mass 2)
# fmm.L2P!(element, tree, branch, harmonics, harmonics_theta)
fmm.L2P!(tree, masses, 3)
u_fmm_12345 = masses[2].potential[1]

dx_42 = masses[4].position - masses[2].position
u_direct_42 = masses[4].strength[1] / sqrt(dx_42' * dx_42)
dx_52 = masses[5].position - masses[2].position
u_direct_52 = masses[5].strength[1] / sqrt(dx_52' * dx_52)

u_direct_12345 = u_direct_123 + u_direct_42 + u_direct_52

@test isapprox(u_direct_123, u_fmm_123; atol=1e-12)
@test isapprox(u_direct_12345, u_fmm_12345; atol=1e-12)

# reset potentials
for i in 1:length(masses)
    masses[i].potential .*= 0
end

# run fmm (reset potentials with reset_tree flag)
fmm.fmm!(tree, masses, theta; reset_tree=true)
u_fmm = [mass.potential[1] for mass in masses]

for mass in masses
    mass.potential .*= 0
end

fmm.direct!(masses)

u_direct = [mass.potential[1] for mass in masses]

for i in 1:length(masses)
    @test isapprox(u_fmm[i], u_direct[i]; atol=1e-12)
end
end

#####
##### vector potential
#####