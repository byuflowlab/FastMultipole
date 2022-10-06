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

#=
@testset "direct" begin

    function V(xi, xj, mj; G=1)
        Rho_ij = xi - xj
        rho_ij = sqrt(Rho_ij' * Rho_ij)
        vij = -G * mj / rho_ij
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

    mass = [Mass(x[:,i], [m[i]], zeros(1), zeros(3)) for i in 1:length(m)]

    fmm.direct!(mass)
    V_tots_direct = [m.potential[1] for m in mass]

    for i in 1:length(V_tots)
        @test isapprox(V_tots[i], V_tots_direct[i]; atol=1e-4)
    end
end

#####
##### define plotting functions
#####
function plot_leaf(elements, leaf, fig_name, save_name;
    stl = "*", clr = "",
    initialize_fig = false,
    save_fig = false
)
    fig = plt.figure(fig_name)
    if initialize_fig
        fig.clear()
        fig.add_subplot(111)
    end
    ax = fig.get_axes()[1]
    n = length(leaf.children)
    Xs = zeros(dims,n)
    for (i,i_index) in enumerate(leaf.children)
        Xs[:,i] .= fmm.get_X(mass,i_index)
    end
    if clr == ""
        ax.scatter(Xs[1,:], Xs[2,:], marker=stl)
    else
        ax.scatter(Xs[1,:], Xs[2,:], marker=stl, color=clr)
    end


    if save_fig
        fig.savefig(save_name)
    end
end

function plot_leaves(elements, leaves, fig_name, save_name;
    stls = ["v", "^", ">", "<", "+", "x", "*"],
    clrs = ["b", "r", "g", "y", "c", "m"],
    initialize_fig = false,
    save_fig = false
)

    for (i,leaf) in enumerate(leaves)
        initialize_fig = initialize_fig && i==1 ? true : false
        stl = stls[i % length(stls) + 1]
        clr = clrs[i % length(clrs) + 1]
        this_save_fig = i == length(leaves) ? save_fig : false
        plot_leaf(elements, leaf, fig_name, save_name;
            stl, clr, initialize_fig, save_fig=this_save_fig
        )
    end
end

function plot_branch(elements, root, level, branch_i, fig_name, save_name;
    stl = "+", clr = "b",
    initialize_fig = false,
    save_fig = false
)
    branch = root.branches[level][branch_i]
    if initialize_fig
        fig = plt.figure(fig_name)
        fig.clear()
        fig.add_subplot(111)
    end

    if level == 2
        leaves = root.branches[1][branch.children]
        plot_leaves(elements, leaves, fig_name, save_name;
            stls = [stl], clrs = [clr],
            initialize_fig = false,
            save_fig
        )
    elseif level > 2
        for branch_i in branch.children
            plot_branch(elements, root, level-1, branch_i, fig_name, save_name;
                stl, clr,
                initialize_fig = false,
                save_fig
            )
        end
    else
        @error "requested plot_branch on level $level; must be >= 2"
    end
end

function plot_branches(elements, root, level, fig_name, save_name;
    stls = ["v", "^", ">", "<", "+", "x", "*"],
    clrs = ["b", "r", "g", "y", "c", "m"],
    initialize_fig = false,
    save_fig = false
)
    if level > 1
        branches = root.branches[level]
        for branch_i in 1:length(branches)
            stl = stls[branch_i % length(stls) + 1]
            clr = clrs[branch_i % length(clrs) + 1]
            this_initialize_fig = branch_i == 1 ? initialize_fig : false
            plot_branch(elements, root, level, branch_i, fig_name, save_name;
                stl, clr,
                initialize_fig = this_initialize_fig,
                save_fig = save_fig
            )
        end
    else
        plot_leaves(elements, root.branches[level], fig_name, save_name;
            stls,
            clrs,
            initialize_fig,
            save_fig
        )
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
    basis = fmm.Cartesian()
    tree = fmm.Tree(masses, basis)

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

@testset "P2M!" begin
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
        mass = [ms[i]]
        potential = zeros(1)
        force = zeros(3)
        masses[i] = Mass(x,mass,potential,force)
    end

    basis = fmm.Cartesian()
    tree = fmm.Tree(masses, basis)
    fmm.upward_pass!(tree, masses, basis)

    tree.branches[1].multipole_expansion .*= 0
    fmm.P2M!(1, tree, masses, basis)

    check_multipole = [
        6.5
        1.735
        1.29
       -0.125
        0.734125
        1.222
       -0.14675
        0.52075
       -0.1775
        0.279125
    ]

    for i in 1:length(check_multipole)
        @test isapprox(check_multipole[i], tree.branches[1].multipole_expansion[i]; atol=1e-6)
    end

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
        mass = [ms[i]]
        potential = zeros(1)
        force = zeros(3)
        masses[i] = Mass(x,mass,potential,force)
    end

    basis = fmm.Cartesian()
    tree = fmm.Tree(masses, basis)
    fmm.upward_pass!(tree, masses, basis)

    leaf_coefficients = [
        [
            1.1
            0.137503025
            0.027503025
            0.082496975
            0.008594128129159
            0.003437953758
            0.01031234874
            0.000343825629159
            0.00206265124
            0.003093523129
        ],
        [
            2.2
            0.604993950
            0.38499395
            -0.16499395
            0.083185836258
            0.1058722775
            -0.0453728825
            0.033686441258
            -0.0288734875
            0.006187046258
        ],
        [
            0.8
        -0.2199978
        -0.1399978
            0.0200022
            0.030249395
            0.03849901
        -0.00550055
            0.012249615
        -0.00350033
            0.000250055
        ],
        [
            0.5
            0.0687479375
        -0.0312520625
        -0.0312520625
            0.0047262789
        -0.0042970297
        -0.00429702967899
            0.0009766914105
            0.001953382821
            0.0009766914105
        ],
        [
            1.9
            0.0712421625
           -0.0237578375
            0.0237473875
            0.0013356436099
           -0.00089082090517
            0.000890429073
            0.0001485354849
           -0.0002969403
            0.0001484048455
        ]
    ]

    # test leaf coefficients
    for (i,i_branch) in enumerate(3:7)
        for i_coeff in 1:10
            @test isapprox(leaf_coefficients[i][i_coeff], tree.branches[i_branch].multipole_expansion[i_coeff])
        end
    end
end

@testset "M2M!" begin
    #####
    ##### M2M!
    #####
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
        mass = [ms[i]]
        potential = zeros(1)
        force = zeros(3)
        masses[i] = Mass(x,mass,potential,force)
    end

    basis = fmm.Cartesian()
    tree = fmm.Tree(masses, basis)
    fmm.upward_pass!(tree, masses, basis)

    # test M2M (branch 2)
    coefficients_2_check = [
        2.4
        0.4699934
        0.2749934
        -0.20000659999999998
        0.04799870749026
        0.0518729510765
        -0.03125074252
        0.016249243750600002
        -0.02687520625
        0.01625055000059
    ]

    for i in 1:length(coefficients_2_check)
        @test isapprox(coefficients_2_check[i], tree.branches[2].multipole_expansion[i]; atol=1e-6)
    end
end

@testset "M2L!" begin
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
        mass = [ms[i]]
        potential = zeros(1)
        force = zeros(3)
        masses[i] = Mass(x,mass,potential,force)
    end

    basis = fmm.Cartesian()
    tree = fmm.Tree(masses, basis)
    fmm.upward_pass!(tree, masses, basis)

    i_local = 7
    j_multipole = 1

    fmm.M2L!(i_local, j_multipole, tree, masses, derivatives, basis)

    local_coeff_check = [
        -25.2752529
        -24.0857454
        -26.1525931539
        -11.2952850056
        -12.71154145
        -42.901452396
        -14.30048412
        -12.71154144
        -14.30048412
         25.4230829
    ]

    for i in 1:length(local_coeff_check)
        @test isapprox(local_coeff_check[i], tree.branches[7].local_expansion[i]; atol=1e-6)
    end
end

@testset "L2L!" begin
    # L2L
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
        mass = [ms[i]]
        potential = zeros(1)
        force = zeros(3)
        masses[i] = Mass(x,mass,potential,force)
    end

    basis = fmm.Cartesian()
    tree = fmm.Tree(masses, basis)
    fmm.upward_pass!(tree, masses, basis)
    fmm.M2L!(2, 1, tree, masses, derivatives, basis)

    multipole_1 = [
        6.5
        1.7350000000000003
        1.29
        -0.12499999999999994
        0.734125
        1.2220000000000002
        -0.14675000000000005
        0.5207500000000002
        -0.1775
        0.27912500000000007
    ]

    local_2 = [
        -29.33342945576157
        -27.321185704286492
        -31.439002480724298
        -44.53273447771192
        0.0
        -60.1478854985297
        -60.1478854985297
        0.0
        -60.1478854985297
        0.0
    ]

    for i in 1:length(local_2)
        @test isapprox(multipole_1[i], tree.branches[1].multipole_expansion[i]; atol=1e-5)
        @test isapprox(local_2[i], tree.branches[2].local_expansion[i]; atol=1e-5)
    end

    check_local_7_addition = [
        -26.2399413
        -27.3211857
        -31.43900248
        -27.991900852
        0.0
        -60.14788549
        -60.14788549
        0.0
        -60.14788549
        0.0
    ]

    local_7_before = deepcopy(tree.branches[7].local_expansion)
    fmm.L2L!(2,tree, basis)
    local_7_addition = tree.branches[7].local_expansion - local_7_before

    for i in 1:length(local_2)
        @test isapprox(local_7_addition[i], check_local_7_addition[i]; atol=1e-5)
    end
end

@testset "L2P!" begin
    # L2P
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
        mass = [ms[i]]
        potential = zeros(1)
        force = zeros(3)
        masses[i] = Mass(x,mass,potential,force)
    end

    basis = fmm.Cartesian()
    tree = fmm.Tree(masses, basis)
    fmm.upward_pass!(tree, masses, basis)
    fmm.M2L!(2, 1, tree, masses, derivatives, basis)
    fmm.L2L!(2, tree, basis)
    fmm.L2P!(7, tree, masses, basis)

    Phi_d = masses[tree.indices[tree.branches[7].first_element]].potential[1]
    check_Phi_d = -25.24935390275

    @test isapprox(Phi_d, check_Phi_d; atol=1e-6)

    # also test using test case branch 7 on branch 3
    tree.branches[3].local_expansion .= [
        -1.90290525167
        1.388471647877
        1.288137857762
        -0.2375652268339
        -0.94431802545980
        -2.9509938277734
        0.590198761262374
        -0.9443180254598
        0.59019876126237
        1.88863605091961
    ] # local expansion due to 7 centered about 3
    masses[2].potential .*= 0
    fmm.L2P!(3, tree, masses, basis)
    phi_b_due2_e = masses[2].potential[1]
    phi_b_due2_e_check = -2.09580318715645

    @test isapprox(phi_b_due2_e, phi_b_due2_e_check; atol=1e-8)
end

@testset "P2P!" begin

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
        mass = [ms[i]]
        potential = zeros(1)
        force = zeros(3)
        masses[i] = Mass(x,mass,potential,force)
    end

    basis = fmm.Cartesian()
    tree = fmm.Tree(masses, basis)

    fmm.P2P!(1,1,tree, masses)
    potential_p2p = [mass.potential[1] for mass in masses]
    for mass in masses; mass.potential .*= 0; end
    fmm.direct!(masses)
    potential_direct = [mass.potential[1] for mass in masses]

    for i in 1:length(potential_p2p)
        @test isapprox(potential_p2p[i], potential_direct[i]; atol=1e-8)
    end
end

@testset "upward pass" begin
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
        mass = [ms[i]]
        potential = zeros(1)
        force = zeros(3)
        masses[i] = Mass(x,mass,potential,force)
    end

    basis = fmm.Cartesian()
    tree = fmm.Tree(masses, basis)

    fmm.P2M!(3, tree, masses, basis)
    fmm.P2M!(4, tree, masses, basis)
    fmm.P2M!(5, tree, masses, basis)
    fmm.P2M!(6, tree, masses, basis)
    fmm.P2M!(7, tree, masses, basis)
    fmm.M2M!(2, tree, basis)
    fmm.M2M!(1, tree, basis)

    masses_2 = Vector{Mass}(undef,length(ms))
    for i in 1:length(ms)
        x = xs[i,:]
        mass = [ms[i]]
        potential = zeros(1)
        force = zeros(3)
        masses_2[i] = Mass(x,mass,potential,force)
    end
    tree_2 = fmm.Tree(masses_2, basis)
    fmm.upward_pass!(tree_2, masses_2, basis)

    for i_branch in 1:length(tree.branches)
        for i_coeff in 1:length(tree.branches[1].multipole_expansion)
            @test isapprox(tree.branches[i_branch].multipole_expansion[i_coeff], tree_2.branches[i_branch].multipole_expansion[i_coeff]; atol=1e-8)
        end
    end
end

@testset "horizontal pass" begin
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
        mass = [ms[i]]
        potential = zeros(1)
        force = zeros(3)
        masses[i] = Mass(x,mass,potential,force)
    end

    basis = fmm.Cartesian()
    tree = fmm.Tree(masses, basis)


    fmm.upward_pass!(tree, masses, basis)
    multipole_7 = deepcopy(tree.branches[7].multipole_expansion)

    multipole_7_check = [
        1.9
        0.0712421625
        -0.0237578375
        0.023747387
        0.0013356436099
        -0.00089082090517
        0.0008904290732765
        0.0001485354849
        -0.000296940301723
        0.000148404845546
    ]

    for i in 1:length(multipole_7)
        @test isapprox(multipole_7[i], multipole_7_check[i]; atol=1e-8)
    end

    # P2P
    fmm.P2P!(3, 3, tree, masses)
    fmm.P2P!(3, 4, tree, masses)
    fmm.P2P!(3, 5, tree, masses)
    fmm.P2P!(4, 3, tree, masses)
    fmm.P2P!(5, 3, tree, masses)
    fmm.P2P!(4, 4, tree, masses)
    fmm.P2P!(5, 4, tree, masses)
    fmm.P2P!(6, 4, tree, masses)
    fmm.P2P!(7, 4, tree, masses)
    fmm.P2P!(4, 5, tree, masses)
    fmm.P2P!(4, 6, tree, masses)
    fmm.P2P!(4, 7, tree, masses)
    fmm.P2P!(5, 5, tree, masses)
    fmm.P2P!(6, 6, tree, masses)
    fmm.P2P!(6, 7, tree, masses)
    fmm.P2P!(7, 6, tree, masses)
    fmm.P2P!(7, 7, tree, masses)

    # M2L

    fmm.M2L!(3, 6, tree, masses, derivatives, basis)

    # checking branch 3's local expansion due to branch 7
    local_3_due2_7_check = [
        -1.90290525167
        1.388471647877
        1.288137857762
        -0.2375652268339
        -0.94431802545980
        -2.9509938277734
        0.590198761262374
        -0.9443180254598
        0.59019876126237
        1.88863605091961
    ]
    local_3_before = deepcopy(tree.branches[3].local_expansion)
    fmm.M2L!(3, 7, tree, masses, derivatives, basis)
    local_3_after = deepcopy(tree.branches[3].local_expansion)
    local_3_due2_7 = local_3_after - local_3_before
    for i in 1:length(local_3_due2_7)
        @test isapprox(local_3_due2_7[i], local_3_due2_7_check[i]; atol=1e-8)
    end

    # resume...
    fmm.M2L!(6, 3, tree, masses, derivatives, basis)
    fmm.M2L!(7, 3, tree, masses, derivatives, basis)
    fmm.M2L!(5, 6, tree, masses, derivatives, basis)
    fmm.M2L!(6, 5, tree, masses, derivatives, basis)
    fmm.M2L!(5, 7, tree, masses, derivatives, basis)
    fmm.M2L!(7, 5, tree, masses, derivatives, basis)

    masses_2 = Vector{Mass}(undef,length(ms))
    for i in 1:length(ms)
        x = xs[i,:]
        mass = [ms[i]]
        potential = zeros(1)
        force = zeros(3)
        masses_2[i] = Mass(x,mass,potential,force)
    end

    tree_2 = fmm.Tree(masses_2, basis)
    theta = 4

    fmm.upward_pass!(tree_2, masses_2, basis)
    fmm.horizontal_pass!(tree_2, masses_2, derivatives, theta, basis)

    for i_branch in 1:length(tree.branches)
        for i_multipole in 1:length(tree.branches[1].multipole_expansion)
            @test isapprox(tree.branches[i_branch].multipole_expansion[i_multipole],
                tree_2.branches[i_branch].multipole_expansion[i_multipole]; atol=1e-8)
        end
        for i_local in 1:length(tree.branches[1].local_expansion)
            @test isapprox(tree.branches[i_branch].local_expansion[i_local],
                tree_2.branches[i_branch].local_expansion[i_local]; atol=1e-8)
        end
    end
end

@testset "downward pass" begin
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
        mass = [ms[i]]
        potential = zeros(1)
        force = zeros(3)
        masses[i] = Mass(x,mass,potential,force)
    end

    basis = fmm.Cartesian()
    tree = fmm.Tree(masses, basis; expansion_order=4)

    theta = 4
    fmm.upward_pass!(tree, masses, basis)
    fmm.horizontal_pass!(tree, masses, derivatives, theta, basis)

    fmm.L2P!(3, tree, masses, basis)
    fmm.L2P!(5, tree, masses, basis)
    fmm.L2P!(6, tree, masses, basis)
    fmm.L2P!(7, tree, masses, basis)

    masses_2 = Vector{Mass}(undef,length(ms))
    for i in 1:length(ms)
        x = xs[i,:]
        mass = [ms[i]]
        potential = zeros(1)
        force = zeros(3)
        masses_2[i] = Mass(x,mass,potential,force)
    end
    tree_2 = fmm.Tree(masses_2, basis; expansion_order=4)

    fmm.upward_pass!(tree_2, masses_2, basis)
    fmm.horizontal_pass!(tree_2, masses_2, derivatives, theta, basis)
    fmm.downward_pass!(tree_2, masses_2, basis)

    for i_branch in 1:length(tree.branches)
        for i_local in 1:length(tree.branches[1].local_expansion)
            @test isapprox(tree.branches[i_branch].local_expansion[i_local],
                tree_2.branches[i_branch].local_expansion[i_local]; atol=1e-8)
        end
        for i_multipole in 1:length(tree.branches[1].local_expansion)
            @test isapprox(tree.branches[i_branch].multipole_expansion[i_multipole],
                tree_2.branches[i_branch].multipole_expansion[i_multipole]; atol=1e-8)
        end
    end

    for i_mass in 1:length(masses)
        @test isapprox(masses[i_mass].potential[1], masses_2[i_mass].potential[1]; atol=1e-8)
    end
end

@testset "derivatives" begin
    # build symbolic kernel
    sym.@variables x y z
    phi = -1/sqrt(x^2 + y^2 + z^2)

    # get derivatives
    symbolic_derivatives = Array{typeof(phi),3}(undef,5,5,5)
    symbolic_derivatives[1,1,1] = phi
    ddx = sym.derivative(phi, x)
    symbolic_derivatives[2,1,1] = ddx
    ddy = sym.derivative(phi, y)
    symbolic_derivatives[1,2,1] = ddy
    ddz = sym.derivative(phi, z)
    symbolic_derivatives[1,1,2] = ddz
    d2dx2 = sym.derivative(ddx, x)
    symbolic_derivatives[3,1,1] = d2dx2
    d2dxdy = sym.derivative(ddx, y)
    symbolic_derivatives[2,2,1] = d2dxdy
    d2dxdz = sym.derivative(ddx, z)
    symbolic_derivatives[2,1,2] = d2dxdz
    d2dy2 = sym.derivative(ddy, y)
    symbolic_derivatives[1,3,1] = d2dy2
    d2dydz = sym.derivative(ddy, z)
    symbolic_derivatives[1,2,2] = d2dydz
    d2dz2 = sym.derivative(ddz, z)
    symbolic_derivatives[1,1,3] = d2dz2
    d3dx3 = sym.derivative(d2dx2, x)
    symbolic_derivatives[4,1,1] = d3dx3
    d3dx2dy = sym.derivative(d2dx2, y)
    symbolic_derivatives[3,2,1] = d3dx2dy
    d3dx2dz = sym.derivative(d2dx2, z)
    symbolic_derivatives[3,1,2] = d3dx2dz
    d3dxdy2 = sym.derivative(d2dxdy, y)
    symbolic_derivatives[2,3,1] = d3dxdy2
    d3dxdydz = sym.derivative(d2dxdy, z)
    symbolic_derivatives[2,2,2] = d3dxdydz
    d3dxdz2 = sym.derivative(d2dxdz, z)
    symbolic_derivatives[2,1,3] = d3dxdz2
    d3dy3 = sym.derivative(d2dy2, y)
    symbolic_derivatives[1,4,1] = d3dy3
    d3dy2dz = sym.derivative(d2dy2, z)
    symbolic_derivatives[1,3,2] = d3dy2dz
    d3dydz2 = sym.derivative(d2dydz, z)
    symbolic_derivatives[1,2,3] = d3dydz2
    d3dz3 = sym.derivative(d2dz2, z)
    symbolic_derivatives[1,1,4] = d3dz3
    d4dx4 = sym.derivative(d3dx3, x)
    symbolic_derivatives[5,1,1] = d4dx4
    d4dx3dy = sym.derivative(d3dx3, y)
    symbolic_derivatives[4,2,1] = d4dx3dy
    d4dx3dz = sym.derivative(d3dx3, z)
    symbolic_derivatives[4,1,2] = d4dx3dz
    d4dx2dy2 = sym.derivative(d3dx2dy, y)
    symbolic_derivatives[3,3,1] = d4dx2dy2
    d4dx2dydz = sym.derivative(d3dx2dy, z)
    symbolic_derivatives[3,2,2] = d4dx2dydz
    d4dx2dz2 = sym.derivative(d3dx2dz, z)
    symbolic_derivatives[3,1,3] = d4dx2dz2
    d4dxdy3 = sym.derivative(d3dxdy2, y)
    symbolic_derivatives[2,4,1] = d4dxdy3
    d4dxdy2dz = sym.derivative(d3dxdy2, z)
    symbolic_derivatives[2,3,2] = d4dxdy2dz
    d4dxdydz2 = sym.derivative(d3dxdz2, y)
    symbolic_derivatives[2,2,3] = d4dxdydz2
    d4dxdz3 = sym.derivative(d3dxdz2, z)
    symbolic_derivatives[2,1,4] = d4dxdz3
    d4dy4 = sym.derivative(d3dy3, y)
    symbolic_derivatives[1,5,1] = d4dy4
    d4dy3dz = sym.derivative(d3dy3, z)
    symbolic_derivatives[1,4,2] = d4dy3dz
    d4dy2dz2 = sym.derivative(d3dy2dz, z)
    symbolic_derivatives[1,3,3] = d4dy2dz2
    d4dydz3 = sym.derivative(d3dydz2, z)
    symbolic_derivatives[1,2,4] = d4dydz3
    d4dz4 = sym.derivative(d3dz3, z)
    symbolic_derivatives[1,1,5] = d4dz4

    # compare
    xyz = (rand(3) .- 0.5) * 5
    s = size(symbolic_derivatives)[1]
    kernel_values = Array{Float64, 3}(undef,s,s,s)
    symbolic_values = Array{Float64, 3}(undef,s,s,s)
    for i = 1:s
        for j = 1:s
            for k = 1:s
                if i+j+k <= s+1
                    this_index = fmm.ijk_2_index(i-1,j-1,k-1)
                    kernel_values[i,j,k] = derivatives[this_index](xyz)
                    r = sym.substitute(symbolic_derivatives[i,j,k],Dict(x=>xyz[1],y=>xyz[2],z=>xyz[3])).val
                    symbolic_values[i,j,k] = r
                    # r2 = mykernel(xyz, i-1, j-1, k-1)
                    # function_values[i,j,k] = r2
                    @test isapprox(kernel_values[i,j,k], symbolic_values[i,j,k]; atol=1e-8)
                    # @test isapprox(kernel_values[i,j,k], function_values[i,j,k]; atol=1e-8)
                end
            end
        end
    end
end
=#

#####
##### spherical based
#####
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

#=
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
fmm.regular_harmonic!(rh_fmm, rho, alpha, beta, P)

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

@testset "P2M" begin
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
    mass = [ms[i]]
    potential = zeros(1)
    force = zeros(3)
    masses[i] = Mass(x,mass,potential,force)
end

basis = fmm.Spherical()
expansion_order = 2
tree = fmm.Tree(masses, basis; expansion_order)

i_mass = 1
i_branch = 5 # use the first mass
fmm.P2M!(i_branch, tree, masses, fmm.Spherical())
center = tree.branches[i_branch].center

x_target = [10.1,-7.3,8.6]
target = Mass(x_target, [1.0], [0.0], zeros(3))
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

@testset "m2m" begin
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
    mass = [ms[i]]
    potential = zeros(1)
    force = zeros(3)
    masses[i] = Mass(x,mass,potential,force)
end

basis = fmm.Spherical()
expansion_order = 3
tree = fmm.Tree(masses, basis; expansion_order)

i_branch = 2 # contains 4th and 5th masses
i_branch_4 = 6 # use the fourth mass
# i_branch_5 = 7 # use the fifth mass
fmm.P2M!(i_branch_4, tree, masses, basis) # evaluate multipole coefficients
# fmm.P2M!(i_branch_5, tree, masses, basis) # evaluate multipole coefficients
fmm.M2M!(i_branch, tree, basis) # translate coefficients to the center of branch 2

x_target = [8.3,1.4,-4.2]
target = Mass(x_target, [0.0], [0.0], zeros(3))
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
