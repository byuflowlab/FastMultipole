using Test
import Statistics
S = Statistics
import Symbolics
sym = Symbolics

import FLOWFMM
fmm = FLOWFMM

import PyPlot
plt = PyPlot
using LaTeXStrings

struct Mass{TF}
    x::Array{TF,1}
    mass::Array{TF,1}
    potential::Array{TF,1}
    force::Array{TF,1}
end

function fmm.kernel!(target::Mass, source::Mass)
    dx = fmm.get_x(target) - fmm.get_x(source)
    r = sqrt(dx' * dx)
    if r > 0
        dV = -source.mass[1] / r
        target.potential[1] += dV
    end
end

function fmm.get_x(mass::Mass)
    mass.x
end

function fmm.get_q(mass::Mass)
    mass.mass[1]
end


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

    kernel = fmm.Gravitational(3)
    V_tots_direct = fmm.direct(mass, kernel)

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
    tree = fmm.Tree(masses)

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

    tree = fmm.Tree(masses)

    fmm.upward_pass!(tree, masses)

    tree.branches[1].multipole_expansion .*= 0
    fmm.P2M!(1, tree, masses)

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

    tree = fmm.Tree(masses)

    fmm.upward_pass!(tree, masses)

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

    tree = fmm.Tree(masses)

    fmm.upward_pass!(tree, masses)

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
        if !isapprox(coefficients_2_check[i], tree.branches[2].multipole_expansion[i]; atol=1e-6); @show i; end
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

    tree = fmm.Tree(masses)

    fmm.upward_pass!(tree, masses)

    i_local = 7
    j_multipole = 1
    kernel = fmm.Gravitational()
    # @show tree.branches[7].local_expansion
    fmm.M2L!(i_local, j_multipole, tree, masses, kernel)
    # @show tree.branches[7].local_expansion

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

    tree = fmm.Tree(masses)
    kernel = fmm.Gravitational()

    fmm.upward_pass!(tree, masses)
    fmm.M2L!(2, 1, tree, masses, kernel)

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
    fmm.L2L!(2,tree)
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

    tree = fmm.Tree(masses)
    kernel = fmm.Gravitational()

    fmm.upward_pass!(tree, masses)
    fmm.M2L!(2, 1, tree, masses, kernel)
    fmm.L2L!(2, tree)
    fmm.L2P!(7, tree, masses)

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
    fmm.L2P!(3, tree, masses)
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

    tree = fmm.Tree(masses)
    kernel = fmm.Gravitational()

    fmm.P2P!(1,1,tree, masses, kernel)
    potential_p2p = [mass.potential[1] for mass in masses]
    potential_direct = fmm.direct(masses, kernel)

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

    tree = fmm.Tree(masses)
    kernel = fmm.Gravitational()

    fmm.P2M!(3, tree, masses)
    fmm.P2M!(4, tree, masses)
    fmm.P2M!(5, tree, masses)
    fmm.P2M!(6, tree, masses)
    fmm.P2M!(7, tree, masses)
    fmm.M2M!(2, tree)
    fmm.M2M!(1, tree)

    masses_2 = Vector{Mass}(undef,length(ms))
    for i in 1:length(ms)
        x = xs[i,:]
        mass = [ms[i]]
        potential = zeros(1)
        force = zeros(3)
        masses_2[i] = Mass(x,mass,potential,force)
    end
    tree_2 = fmm.Tree(masses_2)
    fmm.upward_pass!(tree_2, masses_2)

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

    tree = fmm.Tree(masses)
    kernel = fmm.Gravitational()

    fmm.upward_pass!(tree, masses)
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
    fmm.P2P!(3, 3, tree, masses, kernel)
    fmm.P2P!(3, 4, tree, masses, kernel)
    fmm.P2P!(3, 5, tree, masses, kernel)
    fmm.P2P!(4, 3, tree, masses, kernel)
    fmm.P2P!(5, 3, tree, masses, kernel)
    fmm.P2P!(4, 4, tree, masses, kernel)
    fmm.P2P!(5, 4, tree, masses, kernel)
    fmm.P2P!(6, 4, tree, masses, kernel)
    fmm.P2P!(7, 4, tree, masses, kernel)
    fmm.P2P!(4, 5, tree, masses, kernel)
    fmm.P2P!(4, 6, tree, masses, kernel)
    fmm.P2P!(4, 7, tree, masses, kernel)
    fmm.P2P!(5, 5, tree, masses, kernel)
    fmm.P2P!(6, 6, tree, masses, kernel)
    fmm.P2P!(6, 7, tree, masses, kernel)
    fmm.P2P!(7, 6, tree, masses, kernel)
    fmm.P2P!(7, 7, tree, masses, kernel)

    # M2L

    fmm.M2L!(3, 6, tree, masses, kernel)

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
    fmm.M2L!(3, 7, tree, masses, kernel)
    local_3_after = deepcopy(tree.branches[3].local_expansion)
    local_3_due2_7 = local_3_after - local_3_before
    for i in 1:length(local_3_due2_7)
        @test isapprox(local_3_due2_7[i], local_3_due2_7_check[i]; atol=1e-8)
    end

    # resume...
    fmm.M2L!(6, 3, tree, masses, kernel)
    fmm.M2L!(7, 3, tree, masses, kernel)
    fmm.M2L!(5, 6, tree, masses, kernel)
    fmm.M2L!(5, 7, tree, masses, kernel)
    fmm.M2L!(6, 5, tree, masses, kernel)
    fmm.M2L!(7, 5, tree, masses, kernel)

    masses_2 = Vector{Mass}(undef,length(ms))
    for i in 1:length(ms)
        x = xs[i,:]
        mass = [ms[i]]
        potential = zeros(1)
        force = zeros(3)
        masses_2[i] = Mass(x,mass,potential,force)
    end

    tree_2 = fmm.Tree(masses_2)

    fmm.upward_pass!(tree_2, masses_2)
    fmm.horizontal_pass!(tree_2, masses_2, kernel)

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

    tree = fmm.Tree(masses; expansion_order=4)
    kernel = fmm.Gravitational()

    fmm.upward_pass!(tree, masses)
    fmm.horizontal_pass!(tree, masses, kernel)

    fmm.L2P!(3, tree, masses)
    fmm.L2P!(5, tree, masses)
    fmm.L2P!(6, tree, masses)
    fmm.L2P!(7, tree, masses)

    masses_2 = Vector{Mass}(undef,length(ms))
    for i in 1:length(ms)
        x = xs[i,:]
        mass = [ms[i]]
        potential = zeros(1)
        force = zeros(3)
        masses_2[i] = Mass(x,mass,potential,force)
    end
    tree_2 = fmm.Tree(masses_2; expansion_order=4)

    fmm.upward_pass!(tree_2, masses_2)
    fmm.horizontal_pass!(tree_2, masses_2, kernel)
    fmm.downward_pass!(tree_2, masses_2)

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

    # potential_fmm = [mass.potential[1] for mass in masses]
    # potential_direct = fmm.direct(masses,kernel)
    # @show potential_fmm potential_direct
end

# fmm
#=
n_bodies = 30
xs = [
    0.388392    0.106547   0.845142
    0.240594    0.505334   0.983601
    0.745662    0.0839445  0.922518
    0.869459    0.98942    0.355172
    0.116437    0.392353   0.411578
    0.130627    0.89462    0.215045
    0.726242    0.799088   0.451056
    0.793813    0.488155   0.476661
    0.550357    0.571989   0.847709
    0.929304    0.322994   0.132763
    0.958071    0.700619   0.684907
    0.709772    0.749387   0.801939
    0.77977     0.756256   0.902697
    0.859913    0.111364   0.118537
    0.497385    0.553377   0.691205
    0.169268    0.498795   0.114511
    0.812004    0.279347   0.25
    0.274416    0.410975   0.0954442
    0.227162    0.587327   0.93572
    0.0935186   0.203038   0.161665
    0.91486     0.590517   0.463325
    0.0516916   0.18584    0.927009
    0.767046    0.905011   0.776572
    0.957557    0.989017   0.766284
    0.81422     0.460356   0.0103359
    0.0865304   0.970258   0.753353
    0.524689    0.140235   0.266709
    0.642358    0.453406   0.937758
    0.00710105  0.587438   0.656069
    0.858541    0.426729   0.0686081
]
ms = [
    0.4220728227420478
    0.0747553836759225
    0.4918442683705029
    0.2776964778405684
    0.728065982900099
    0.26036166706827624
    0.40947592030106006
    0.2561866650710254
    0.06438382582261348
    0.6679735476119923
    0.04889662377468129
    0.1550291511132713
    0.16037100929047354
    0.18190416799126385
    0.9021939971600423
    0.6684499952695557
    0.7967050012287675
    0.43663719038929916
    0.06761496973869541
    0.43105003179953894
    0.4550424416035894
    0.9582284464593078
    0.39984740708504507
    0.4488960100421002
    0.6698005050460909
    0.001771998875738312
    0.6259215534257303
    0.33877748744249203
    0.6120473836581699
    0.197476919815341
]
masses = Vector{Mass}(undef,length(ms))
for i in 1:length(ms)
    x = xs[i,:]
    mass = [ms[i]]
    potential = zeros(1)
    force = zeros(3)
    masses[i] = Mass(x,mass,potential,force)
end

kernel = fmm.Gravitational()
tree = fmm.fmm!(masses, kernel, 1, 1)

potential_fmm = [mass.potential[1] for mass in masses]
@show potential_fmm
for i in 1:length(masses)
    masses[i].potential .*= 0
end
potential_direct = fmm.direct(masses, kernel)

err = potential_fmm - potential_direct
@show err
# rel_err = err ./ potential_direct
=#


#=
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

function benchmark_fmm(tree::fmm.Tree, elements, expansion_order::Int, true_potential)
    fmm.change_expansion_order!(tree, expansion_order)
    for i in 1:length(elements)
        elements[i].potential .*= 0
    end
    time = @elapsed fmm.fmm!(tree, masses, kernel)
    potential = [masses[i].potential[1] for i in 1:length(masses)]
    err = abs.((potential - true_potential) ./ true_potential)
    max_err = maximum(err)
    min_err = minimum(err)
    mean_err = S.mean(err)
    return time, max_err, min_err, mean_err, potential
end


function benchmark_fmm(tree, orders, elements)
    check_potential = fmm.direct(elements, kernel)

    times = zeros(length(orders))
    max_errs = zeros(length(orders))
    min_errs = zeros(length(orders))
    mean_errs = zeros(length(orders))
    potentials = zeros(length(elements), length(orders))

    for (i,expansion_order) in enumerate(orders)
        time, max_err, min_err, mean_err, potential = benchmark_fmm(tree, elements, expansion_order, check_potential)
        times[i] = time
        max_errs[i] = max_err
        min_errs[i] = min_err
        mean_errs[i] = mean_err
        potentials[:,i] .= potential
    end

    return times, max_errs, min_errs, mean_errs, potentials
end

kernel = fmm.Gravitational(; order=4)
masses = Vector{Mass}(undef,length(ms))
for i in 1:length(ms)
    local x = xs[i,:]
    mass = [ms[i]]
    potential = zeros(1)
    force = zeros(3)
    masses[i] = Mass(x,mass,potential,force)
end

tree = fmm.Tree(masses; expansion_order=4)
orders = 1:4
times, max_errs, min_errs, mean_errs, potentials = benchmark_fmm(tree, orders, masses)

fig = plt.figure("benchmark_fmm")
fig.clear()
fig.add_subplot(111,xlabel="order",ylabel="rel. err.")
ax = fig.get_axes()[1]
ax.plot(orders, mean_errs)

fig = plt.figure("convergence")
fig.clear()
fig.add_subplot(111,xlabel="order", ylabel="potential")
ax = fig.get_axes()[1]
for i in 1:length(masses)
    ax.plot(orders,potentials[i,:], label="mass $i, fmm")
    ax.plot(orders,true_potentials[i]*ones(length(orders)), label="mass $i, direct")
end
ax.legend()

fmm.fmm!(tree, masses, kernel)
fmm_potentials = [masses[i].potential[1] for i in 1:length(masses)]
true_potentials = fmm.direct(masses, kernel)
=#

@testset "derivatives" begin
    kernel = fmm.Gravitational(3) # 2-D

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
    # d5dx5 = sym.derivative(d4dx4, x)
    # symbolic_derivatives[6,1,1] = d5dx5
    # d5dx4dy = sym.derivative(d4dx4, y)
    # symbolic_derivatives[5,2] = d5dx4dy
    # d5dx3dy2 = sym.derivative(d4dx3dy, y)
    # symbolic_derivatives[4,3] = d5dx3dy2
    # d5dx2dy3 = sym.derivative(d4dx2dy2, y)
    # symbolic_derivatives[3,4] = d5dx2dy3
    # d5dxdy4 = sym.derivative(d4dxdy3, y)
    # symbolic_derivatives[2,5] = d5dxdy4
    # d5dy5 = sym.derivative(d4dy4, y)
    # symbolic_derivatives[1,6] = d5dy5

    # compare
    xyz = (rand(3) .- 0.5) * 5
    s = size(kernel.potential_derivatives)[1]
    kernel_values = Array{Float64, 3}(undef,s,s,s)
    symbolic_values = Array{Float64, 3}(undef,s,s,s)
    for i = 1:s
        for j = 1:s
            for k = 1:s
                if i+j+k <= s+1
                    kernel_values[i,j,k] = kernel.potential_derivatives[i,j,k](xyz, 1.0, 1.0)
                    r = sym.substitute(symbolic_derivatives[i,j,k],Dict(x=>xyz[1],y=>xyz[2],z=>xyz[3])).val
                    symbolic_values[i,j,k] = r
                    @test isapprox(kernel_values[i,j,k], symbolic_values[i,j,k]; atol=1e-8)
                end
            end
        end
    end

end
