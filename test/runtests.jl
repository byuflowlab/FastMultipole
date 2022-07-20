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

grav = fmm.Gravitational()


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

    # println("i | j | Rho_ij | rho_ij | V_ij")
    # for i in 1:length(m)
    #     for j in 1:length(m)
    #         println("$i | $j | $(Rho_ijs[:,i,j]) | $(rho_ijs[i,j]) | $(V_ijs[i,j])")
    #     end
    # end

    V_tots = zeros(length(m))
    # println("i | netV_i")
    # println("--- | ---")
    for i = 1:length(m)
        V_tots[i] = sum(V_ijs[i,:])
        # println("$i | $(V_tots[i])")
    end

    mass = [fmm.Mass(x[:,i], m[i]) for i in 1:length(m)]

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

# @testset "tree" begin

    xs = [
            0.5832647965804127 0.7497478131912669 0.19859728506727925 0.6051066800624201 0.2578237363332738 0.6386521914280865 0.9474258695376094 0.7867488168628627 0.10126976512042862 0.44443574016766707 0.47418866654956693 0.12532503316242738 0.07993348203423878 0.9579369962788504 0.8884127880963202 0.5548496272645929 0.993112692800262 0.1110670125629094 0.8695311380557229 0.2626250154929006 0.10803109341561612 0.39287820224805303 0.9184451751721363 0.8649598210204277 0.3361347448400116 0.5759351701024353 0.1143830935594854 0.23688814330076702 0.20950034640560378 0.14488252317537365 0.9845579041732264 0.34738507454862955 0.7517099911123692 0.9468304995704726 0.6769347452716619 0.8654390901912041 0.6399898275602325 0.2876774212674882 0.08818030258475607 0.958941787673953 0.6729264741895464 0.0224541236530289 0.7706044769245111 0.21401241634602308 0.2571270961431722 0.7137570627228023 0.8612958357582063 0.3796785984756783 0.05415400054250363 0.7762418689479329 0.3421931509773748 0.6527952292238293 0.28928740230330696 0.14999587246136015 0.4756477306213436 0.3912464303222212 0.9497649194438649 0.9670172984149368 0.9814863413383381 0.32767784422388124;
            0.045652586829020514 0.26780811909299684 0.3277259473919307 0.2013166777338533 0.02811958218487831 0.08670807712515072 0.821149472916209 0.7534896315470452 0.7652433229686657 0.878576633803519 0.39741718447966834 0.8185504067547684 0.8772176175209456 0.43631461744150446 0.4260742735132981 0.8038823878723291 0.9664746555405714 0.4360562609113938 0.7376303978710215 0.13982290491431582 0.6672035071628828 0.658277553730098 0.09765290824207917 0.7225044290855902 0.7377862346784358 0.18153780845685752 0.46926799634217486 0.642753827507496 0.8029456959626464 0.5111896870492305 0.11889103183916405 0.20855761515425786 0.7475128100558941 0.9295458626517763 0.1394687618702517 0.33213284460615244 0.46653042330122885 0.1162665005571184 0.7121415813452903 0.9503882988315557 0.5891410028844561 0.7007200924605026 0.37420441551384753 0.8894182036349167 0.4139268777109171 0.3363760779498075 0.40932406231638363 0.10660561998439544 0.8383075523907313 0.6251521172261196 0.9360624271789464 0.563061032685646 0.08169190818003691 0.22796946022740427 0.981004563674444 0.5402455961049075 0.6121572652108134 0.44421679997831465 0.992071818024292 0.6676005622654677
        ]
    ms = [
        0.10773865335658317
        0.37609596577279225
        0.46033829195567444
        0.8212600808301382
        0.6440665252617324
        0.6441070414941334
        0.4549585375002274
        0.38823325161217537
        0.219853728520913
        0.6645239246041219
        0.1723442610723276
        0.3727230729059432
        0.15869742574998003
        0.932795040481947
        0.24835826831831653
        0.8190184527329007
        0.10941759749178681
        0.9569475714975166
        0.5935537135078579
        0.18582705463746185
        0.6308497640802049
        0.42584681604691244
        0.021696619502994396
        0.17251550691570827
        0.6115986281409824
        0.8740516771548767
        0.862714916596107
        0.49684259078769744
        0.643240003848714
        0.19071908754724132
        0.6316500185768179
        0.5527737147205465
        0.6948536910765557
        0.30115735121670917
        0.8809059168627882
        0.26645124762120154
        0.24028337817814904
        0.30960208260723365
        0.25042853861159275
        0.7894468676059447
        0.16620452772936112
        0.35935745492582005
        0.06453039827550322
        0.20996193491364168
        0.7450170960850135
        0.7135301510689132
        0.9528230709796315
        0.1696701978810904
        0.83446066994928
        0.7179771538474993
        0.8086417332975755
        0.3497532991319767
        0.6820726796262995
        0.558194355491491
        0.130118738780459
        0.36715413498723826
        0.04534098084022564
        0.22450149797711694
        0.7701710369032886
        0.8814790673893307
    ]

    mass = [fmm.Mass(xs[:,i], ms[i]) for i in 1:length(ms)]
    dims = size(xs)[1]

    # max_branches = 4
    branch_limit = 4
    elements_per_leaf = 3
    p_expansion = 3

    root = fmm.Root(mass, branch_limit, elements_per_leaf, p_expansion)
    leaves = root.branches[1]

    # ensure all elements are accounted for at the leaf level
    check_indices = Int32[]
    for leaf in leaves
        for i in leaf.children
            push!(check_indices, i)
        end
    end

    sort!(check_indices)
    @test check_indices[end] == length(mass)

    # ensure all elements are accounted for at level 1
    check_indices_2 = Int32[]
    for branch in root.branches[2]
        for leaf in root.branches[1][branch.children]
            for i in leaf.children
                push!(check_indices_2, i)
            end
        end
    end

    @test length(unique(check_indices_2)) == length(mass)

    # plot_leaves(mass, root.branches[1], "leaves", "leaves.png";
    #         stls = ["v", "^", ">", "<", "+", "x", "*"],
    #         initialize_fig = true,
    #         save_fig = true
    #     )

    # plot_branches(mass, root, 2, "branches", "branches.png";
    #         stls = ["v", "^", ">", "<", "+", "x", "*"],
    #         clrs = ["b", "r", "g", "y", "c", "m"],
    #         initialize_fig = true,
    #         save_fig = true
    #     )

    # plot_branches(mass, root, 3, "top level", "top_level.png";
    #         stls = ["v", "^", ">", "<", "+", "x", "*"],
    #         clrs = ["b", "r", "g", "y", "c", "m"],
    #         initialize_fig = true,
    #         save_fig = true
    #     )

# end

@testset "interaction list" begin

    xs = [
            0.5832647965804127 0.7497478131912669 0.19859728506727925 0.6051066800624201 0.2578237363332738 0.6386521914280865 0.9474258695376094 0.7867488168628627 0.10126976512042862 0.44443574016766707 0.47418866654956693 0.12532503316242738 0.07993348203423878 0.9579369962788504 0.8884127880963202 0.5548496272645929 0.993112692800262 0.1110670125629094 0.8695311380557229 0.2626250154929006 0.10803109341561612 0.39287820224805303 0.9184451751721363 0.8649598210204277 0.3361347448400116 0.5759351701024353 0.1143830935594854 0.23688814330076702 0.20950034640560378 0.14488252317537365 0.9845579041732264 0.34738507454862955 0.7517099911123692 0.9468304995704726 0.6769347452716619 0.8654390901912041 0.6399898275602325 0.2876774212674882 0.08818030258475607 0.958941787673953 0.6729264741895464 0.0224541236530289 0.7706044769245111 0.21401241634602308 0.2571270961431722 0.7137570627228023 0.8612958357582063 0.3796785984756783 0.05415400054250363 0.7762418689479329 0.3421931509773748 0.6527952292238293 0.28928740230330696 0.14999587246136015 0.4756477306213436 0.3912464303222212 0.9497649194438649 0.9670172984149368 0.9814863413383381 0.32767784422388124;
            0.045652586829020514 0.26780811909299684 0.3277259473919307 0.2013166777338533 0.02811958218487831 0.08670807712515072 0.821149472916209 0.7534896315470452 0.7652433229686657 0.878576633803519 0.39741718447966834 0.8185504067547684 0.8772176175209456 0.43631461744150446 0.4260742735132981 0.8038823878723291 0.9664746555405714 0.4360562609113938 0.7376303978710215 0.13982290491431582 0.6672035071628828 0.658277553730098 0.09765290824207917 0.7225044290855902 0.7377862346784358 0.18153780845685752 0.46926799634217486 0.642753827507496 0.8029456959626464 0.5111896870492305 0.11889103183916405 0.20855761515425786 0.7475128100558941 0.9295458626517763 0.1394687618702517 0.33213284460615244 0.46653042330122885 0.1162665005571184 0.7121415813452903 0.9503882988315557 0.5891410028844561 0.7007200924605026 0.37420441551384753 0.8894182036349167 0.4139268777109171 0.3363760779498075 0.40932406231638363 0.10660561998439544 0.8383075523907313 0.6251521172261196 0.9360624271789464 0.563061032685646 0.08169190818003691 0.22796946022740427 0.981004563674444 0.5402455961049075 0.6121572652108134 0.44421679997831465 0.992071818024292 0.6676005622654677
        ]
    ms = [
        0.10773865335658317
        0.37609596577279225
        0.46033829195567444
        0.8212600808301382
        0.6440665252617324
        0.6441070414941334
        0.4549585375002274
        0.38823325161217537
        0.219853728520913
        0.6645239246041219
        0.1723442610723276
        0.3727230729059432
        0.15869742574998003
        0.932795040481947
        0.24835826831831653
        0.8190184527329007
        0.10941759749178681
        0.9569475714975166
        0.5935537135078579
        0.18582705463746185
        0.6308497640802049
        0.42584681604691244
        0.021696619502994396
        0.17251550691570827
        0.6115986281409824
        0.8740516771548767
        0.862714916596107
        0.49684259078769744
        0.643240003848714
        0.19071908754724132
        0.6316500185768179
        0.5527737147205465
        0.6948536910765557
        0.30115735121670917
        0.8809059168627882
        0.26645124762120154
        0.24028337817814904
        0.30960208260723365
        0.25042853861159275
        0.7894468676059447
        0.16620452772936112
        0.35935745492582005
        0.06453039827550322
        0.20996193491364168
        0.7450170960850135
        0.7135301510689132
        0.9528230709796315
        0.1696701978810904
        0.83446066994928
        0.7179771538474993
        0.8086417332975755
        0.3497532991319767
        0.6820726796262995
        0.558194355491491
        0.130118738780459
        0.36715413498723826
        0.04534098084022564
        0.22450149797711694
        0.7701710369032886
        0.8814790673893307
    ]

    mass = [fmm.Mass(xs[:,i], ms[i]) for i in 1:length(ms)]
    dims = size(xs)[1]

    branch_limit = 4
    elements_per_leaf = 1
    p_expansion = 4

    root = fmm.Root(mass, branch_limit, elements_per_leaf, p_expansion)

    # plot
    # plot_branches(mass, root, 4, "four_levels_4", "four_levels_4.png";
    #     stls = ["v", "^", ">", "<", "+", "x", "*"],
    #     clrs = ["b", "r", "g", "y", "c", "m"],
    #     initialize_fig = true,
    #     save_fig = true
    # )
    # plot_branches(mass, root, 3, "four_levels_3", "four_levels_3.png";
    #     stls = ["v", "^", ">", "<", "+", "x", "*"],
    #     clrs = ["b", "r", "g", "y", "c", "m"],
    #     initialize_fig = true,
    #     save_fig = true
    # )
    # plot_branches(mass, root, 2, "four_levels_2", "four_levels_2.png";
    #     stls = ["v", "^", ">", "<", "+", "x", "*"],
    #     clrs = ["b", "r", "g", "y", "c", "m"],
    #     initialize_fig = true,
    #     save_fig = true
    # )
    # plot_branches(mass, root, 1, "four_levels_1", "four_levels_1.png";
    #     stls = ["v", "^", ">", "<", "+", "x", "*"],
    #     clrs = ["b", "r", "g", "y", "c", "m"],
    #     initialize_fig = true,
    #     save_fig = true
    # )

    level = 1
    ci = CartesianIndex((3,2))
    nn = fmm.i_nearest_neighbors(root, level, ci)
    il = fmm.interaction_list(root, level, ci)

    # 8x8 grid for (3,2) should contain:

    il_check = CartesianIndex.([
        (1,1), (1,2), (1,3), (1,4),
        (2,4), (3,4), (4,4),
        (5,1), (5,2), (5,3), (5,4),
        (6,1), (6,2), (6,3), (6,4)
    ])

    for ci in il
        @test in(ci, il_check)
    end
end

@testset "upward pass" begin # test multipole coefficients
    xs = [
            0.5832647965804127 0.7497478131912669 0.19859728506727925 0.6051066800624201 0.2578237363332738 0.6386521914280865 0.9474258695376094 0.7867488168628627 0.10126976512042862 0.44443574016766707 0.47418866654956693 0.12532503316242738 0.07993348203423878 0.9579369962788504 0.8884127880963202 0.5548496272645929 0.993112692800262 0.1110670125629094 0.8695311380557229 0.2626250154929006 0.10803109341561612 0.39287820224805303 0.9184451751721363 0.8649598210204277 0.3361347448400116 0.5759351701024353 0.1143830935594854 0.23688814330076702 0.20950034640560378 0.14488252317537365 0.9845579041732264 0.34738507454862955 0.7517099911123692 0.9468304995704726 0.6769347452716619 0.8654390901912041 0.6399898275602325 0.2876774212674882 0.08818030258475607 0.958941787673953 0.6729264741895464 0.0224541236530289 0.7706044769245111 0.21401241634602308 0.2571270961431722 0.7137570627228023 0.8612958357582063 0.3796785984756783 0.05415400054250363 0.7762418689479329 0.3421931509773748 0.6527952292238293 0.28928740230330696 0.14999587246136015 0.4756477306213436 0.3912464303222212 0.9497649194438649 0.9670172984149368 0.9814863413383381 0.32767784422388124;
            0.045652586829020514 0.26780811909299684 0.3277259473919307 0.2013166777338533 0.02811958218487831 0.08670807712515072 0.821149472916209 0.7534896315470452 0.7652433229686657 0.878576633803519 0.39741718447966834 0.8185504067547684 0.8772176175209456 0.43631461744150446 0.4260742735132981 0.8038823878723291 0.9664746555405714 0.4360562609113938 0.7376303978710215 0.13982290491431582 0.6672035071628828 0.658277553730098 0.09765290824207917 0.7225044290855902 0.7377862346784358 0.18153780845685752 0.46926799634217486 0.642753827507496 0.8029456959626464 0.5111896870492305 0.11889103183916405 0.20855761515425786 0.7475128100558941 0.9295458626517763 0.1394687618702517 0.33213284460615244 0.46653042330122885 0.1162665005571184 0.7121415813452903 0.9503882988315557 0.5891410028844561 0.7007200924605026 0.37420441551384753 0.8894182036349167 0.4139268777109171 0.3363760779498075 0.40932406231638363 0.10660561998439544 0.8383075523907313 0.6251521172261196 0.9360624271789464 0.563061032685646 0.08169190818003691 0.22796946022740427 0.981004563674444 0.5402455961049075 0.6121572652108134 0.44421679997831465 0.992071818024292 0.6676005622654677
        ]
    ms = [
        0.10773865335658317
        0.37609596577279225
        0.46033829195567444
        0.8212600808301382
        0.6440665252617324
        0.6441070414941334
        0.4549585375002274
        0.38823325161217537
        0.219853728520913
        0.6645239246041219
        0.1723442610723276
        0.3727230729059432
        0.15869742574998003
        0.932795040481947
        0.24835826831831653
        0.8190184527329007
        0.10941759749178681
        0.9569475714975166
        0.5935537135078579
        0.18582705463746185
        0.6308497640802049
        0.42584681604691244
        0.021696619502994396
        0.17251550691570827
        0.6115986281409824
        0.8740516771548767
        0.862714916596107
        0.49684259078769744
        0.643240003848714
        0.19071908754724132
        0.6316500185768179
        0.5527737147205465
        0.6948536910765557
        0.30115735121670917
        0.8809059168627882
        0.26645124762120154
        0.24028337817814904
        0.30960208260723365
        0.25042853861159275
        0.7894468676059447
        0.16620452772936112
        0.35935745492582005
        0.06453039827550322
        0.20996193491364168
        0.7450170960850135
        0.7135301510689132
        0.9528230709796315
        0.1696701978810904
        0.83446066994928
        0.7179771538474993
        0.8086417332975755
        0.3497532991319767
        0.6820726796262995
        0.558194355491491
        0.130118738780459
        0.36715413498723826
        0.04534098084022564
        0.22450149797711694
        0.7701710369032886
        0.8814790673893307
    ]

    mass = [fmm.Mass(xs[:,i], ms[i]) for i in 1:length(ms)]
    dims = size(xs)[1]

    branch_limit = 4
    elements_per_leaf = 1
    p_expansion = 4

    root = fmm.Root(mass, branch_limit, elements_per_leaf, p_expansion)


    fmm.upward_pass!(root, mass)

    # test total charge
    @test isapprox(fmm.get_q(root.branches[end][1].element), sum(ms); atol=1e-8)

    # test multipole coefficients
    xs = [
        0.25 0.75 1.4 2.4 3.7 3.3;
        3.25 3.75 2.4 1.3 0.8 0.4;
    ]

    ms = [
        1.42, 0.89, 0.76, 1.10, 0.44, 0.63
    ]

    mass = [fmm.Mass(xs[:,i], ms[i]) for i in 1:length(ms)]
    dims = size(xs)[1]
    branch_limit = 4
    elements_per_leaf = 1
    p_expansion = 4
    n_divides = 3
    rect = [0.0 0.0; 4.0 4.0]
    kernel = fmm.Gravitational(2)
    root = fmm.Root(mass, branch_limit, elements_per_leaf, p_expansion; n_divides, rect)

    fmm.upward_pass!(root, mass)

    #####
    ##### test leaf level
    #####
    cells_check = CartesianIndex.([(1, 7), (2, 8), (3, 5), (5, 3), (8, 2), (7, 1)])
    multipole_coefficients_check = [
        [
        -0.0 - 0.0im
        -0.0 - 0.0im
        -0.0 - 0.0im
        -0.0 - 0.0im
        ],
        [
        -0.0 - 0.0im
        -0.0 - 0.0im
        -0.0 - 0.0im
        -0.0 - 0.0im
        ],
        [
        -0.11399999999999993 - 0.11399999999999993im
                        -0.0 - 0.01709999999999998im
        0.0017099999999999971 - 0.0017099999999999971im
        0.00038474999999999913 - 0.0im
        ],
        [
        -0.16499999999999992 - 0.055000000000000056im
        -0.010999999999999984 - 0.008250000000000004im
        -0.0008249999999999973 - 0.0011916666666666662im
        -4.812499999999963e-5 - 0.0001649999999999998im
        ],
        [
        0.021999999999999922 - 0.02200000000000002im
        4.865899350114944e-18 + 0.0010999999999999972im
        -3.666666666666676e-5 - 3.6666666666666286e-5im
        2.749999999999986e-6 - 2.4329496750574655e-20im
        ],
        [
        -0.03149999999999989 - 0.09450000000000001im
        0.006300000000000008 - 0.004724999999999984im
        0.0006824999999999981 + 0.00047250000000000184im
        -2.7562500000000396e-5 + 9.44999999999998e-5im
        ]
    ]
    for i_cell in 1:length(multipole_coefficients_check)
        for i in 1:p_expansion
            @test isapprox(root.branches[1][cells_check[i_cell]].multipole_coefficients[i], multipole_coefficients_check[i_cell][i]; atol=1e-8)
        end
    end

    #####
    ##### 2nd level
    #####
    cells_check = CartesianIndex.([(1, 4), (2, 3), (3, 2), (4, 1)])
    multipole_coefficients_check = [
        [
           0.13249999999999998 + 0.13249999999999998im
                          -0.0 - 0.144375im
        -0.0055208333333333325 + 0.0055208333333333325im
                  0.0090234375 - 0.0im
        ],
        [
           0.07600000000000007 + 0.07600000000000007im
                           0.0 - 0.0076000000000000095im
        -0.0005066666666666683 + 0.0005066666666666683im
          7.600000000000055e-5 + 0.0im
        ],
        [
            0.1100000000000001 + 0.21999999999999997im
          0.016499999999999987 - 0.02200000000000002im
         -0.004033333333333335 - 0.0007333333333333272im
         0.0001925000000000017 + 0.00066im
        ],
        [
          0.038000000000000034 - 0.06900000000000002im
         0.0015499999999999663 - 0.03900000000000004im
          0.007166666666666677 + 0.0009899999999999948im
         0.0014192500000000038 + 0.000941999999999999im
        ]
    ]
    for i_cell in 1:length(multipole_coefficients_check)
        for i in 1:p_expansion
            @test isapprox(root.branches[2][cells_check[i_cell]].multipole_coefficients[i], multipole_coefficients_check[i_cell][i]; atol=1e-8)
        end
    end

    # now check the evaluation of the multipole expansions
    # consider the influence of the (4,1) cell, which consists of just 2 particles
    target = [1.9, 2.3]
    phi_41_direct = fmm.direct(mass[end-1:end], target, kernel)

    function multipole(branch, target)
        dx_complex = fmm.complexify(target - fmm.get_X(branch))
        V = fmm.get_q(branch) * log(dx_complex)
        for k in 1:length(branch.multipole_coefficients)
            V += branch.multipole_coefficients[k] / dx_complex^k
        end
        return real(V)
    end

    phi_41_multipole = multipole(root.branches[2][4,1], target)

    @test isapprox(phi_41_direct, phi_41_multipole; atol=1e-4)
end

@testset "downward sweep" begin

    xs = [
        0.25 0.75 1.4 2.4 3.7 3.3;
        3.25 3.75 2.4 1.3 0.8 0.4;
    ]

    ms = [
        1.42, 0.89, 0.76, 1.10, 0.44, 0.63
    ]

    mass = [fmm.Mass(xs[:,i], ms[i]) for i in 1:length(ms)]
    dims = size(xs)[1]
    branch_limit = 4
    elements_per_leaf = 1
    p_expansion = 4
    n_divides = 3
    rect = [0.0 0.0; 4.0 4.0]
    kernel = fmm.Gravitational(2)
    root = fmm.Root(mass, branch_limit, elements_per_leaf, p_expansion; n_divides, rect)

    fmm.upward_pass!(root, mass)
    fmm.downward_pass!(root, mass, kernel)

    #####
    ##### second level (interaction lists are empty at 3rd and 4th)
    #####
    cells_check = CartesianIndex.([(1, 4), (2, 3), (3, 2), (4, 1)])
    local_coefficients_check = [
        [
              2.704544249875981 + 5.037790035229733im
            -0.4287457799318416 - 0.4660806557034466im
           0.009193899498456792 - 0.09472583902659894im
           0.015934137881069246 - 0.011435465184452875im
           0.004533099657060751 + 0.0010950858566779556im
        ],
        [
           1.0908589448507957 + 2.5288641149224755im
            -0.27345589453125 - 0.26970611328125im
       -0.0008782421874999977 - 0.0682271728515625im
         0.011109712727864585 - 0.01150348876953125im
         0.004169183044433594 - 8.29418945312495e-5im
        ],
        [
             2.419660864429273 - 1.7483648095314388im
            0.5849212646484375 + 0.5523138427734374im
         -0.008022460937499998 - 0.1377836608886719im
         -0.019932276407877605 + 0.02385724385579427im
          0.007690629959106445 - 0.0009542846679687499im
        ],
        [
             4.137508961926001 - 2.329139599284487im
            0.5886573854166668 + 0.5550849008487655im
         -0.007145874485596711 - 0.10999056777263375im
         -0.013059289676926154 + 0.016213466345879057im
         0.0045568182794460165 - 0.0007179117726909012im
        ]
    ]
    for i_cell in 1:length(local_coefficients_check)
        for i in 1:p_expansion+1
            @test isapprox(root.branches[2][cells_check[i_cell]].local_coefficients[i], local_coefficients_check[i_cell][i]; atol=1e-8)
        end
    end

    #####
    ##### leaf (first) level
    #####
    cells_check = CartesianIndex.([(1, 7), (2,8), (3,5)])#, (5,3), (8,2), (7,1)])
    local_coefficients_check = [
        [
                2.9789707557638114 + 7.165626283309736im
               -0.9030737015609823 - 0.7336940109590224im
              -0.06658961719319323 - 0.2723559563062657im
               0.04139789491401981 - 0.09869042153152492im
              0.042187240282060764 - 0.024198664143322033im
        ],
        [
              3.032813634912179 + 6.350909092105569im
            -0.5932570323199999 - 0.95943551776im
                  0.13441210368 - 0.22027076975999996im
            0.09276878395733333 + 0.010910466432000034im
           0.012969991084799983 + 0.037797178566400004im
        ],
        [
              2.4399976877728635 + 3.13318739308609im
            0.017141812499999853 + 0.5428188125im
             0.07794314999999999 - 0.7693699437500001im
            -0.12016757624999996 - 0.03694487708333337im
             0.13458028156250001 + 0.013450870599999996im
        ]
    ]
    for i_cell in 1:length(local_coefficients_check)
        for i in 1:p_expansion+1
            @test isapprox(root.branches[1][cells_check[i_cell]].local_coefficients[i], local_coefficients_check[i_cell][i]; atol=1e-5)
        end
    end

end

function err_2_table(direct, fmm)
    println("Direct | FMM | Err. | % Err.")
    println("--- | --- | --- | ---")
    for i in 1:length(direct)
        println(direct[i], " | ", fmm[i], " | ", direct[i] - fmm[i], " | ", (direct[i] - fmm[i])/direct[i])
    end
end

function test_fmm(p_expansion=4)
    xs = [
        0.25 0.75 1.4 2.4 3.7 3.3;
        3.25 3.75 2.4 1.3 0.8 0.4;
    ]

    ms = [
        1.42, 0.89, 0.76, 1.10, 0.44, 0.63
    ]

    mass = [fmm.Mass(xs[:,i], ms[i]) for i in 1:length(ms)]
    dims = size(xs)[1]
    branch_limit = 4
    elements_per_leaf = 1
    n_divides = 3
    rect = [0.0 0.0; 4.0 4.0]
    kernel = fmm.Gravitational(2)
    root = fmm.Root(mass, branch_limit, elements_per_leaf, p_expansion; n_divides, rect)

    fmm.upward_pass!(root, mass)
    fmm.downward_pass!(root, mass, kernel)

    fmm.evaluate!(root, mass, kernel)

    V_tots_direct = fmm.direct(mass, kernel)
    V_tots_fmm = [fmm.get_V(mass,i) for i in 1:length(mass)]

    return V_tots_direct, V_tots_fmm
end

@testset "eval" begin
    V_tots_direct, V_tots_fmm = test_fmm()

    for i in 1:length(V_tots_direct)
        @test isapprox(V_tots_direct[i], V_tots_fmm[i]; atol=1e-5)
    end
end

function create_root(n_elements, p_expansion;
        kernel = fmm.Gravitational(2)
    )
    xs = rand(2,n_elements)
    ms = rand(n_elements)
    mass = [fmm.Mass(xs[:,i], ms[i]) for i in 1:length(ms)]
    branch_limit = 4
    elements_per_leaf = 1
    @time root = fmm.Root(mass, branch_limit, elements_per_leaf, p_expansion)
    return root, mass, kernel
end

function evaluate_time(root, mass, kernel; run_direct=true, run_fmm=true)
    fmm.reset_potential!(mass)
    if run_direct
        t_direct = @elapsed fmm.direct(mass, kernel)
    else
        t_direct = NaN
    end
    if run_fmm
        t_fmm = @elapsed fmm.fmm!(root, mass, kernel)
    else
        t_fmm = NaN
    end
    return t_direct, t_fmm
end

function benchmark_fmm!(add_elements; run_direct=true, run_fmm=true, p_expansion=4)
    n_element_list = Int64[]
    ts_direct = Float64[]
    ts_fmm = Float64[]
    benchmark_fmm!(n_element_list, ts_direct, ts_fmm, add_elements; run_direct, run_fmm, p_expansion)
    return n_element_list, ts_direct, ts_fmm
end

function benchmark_fmm!(n_element_list, ts_direct, ts_fmm, add_elements; run_direct=true, run_fmm=true, p_expansion=4)
    for n_elements in add_elements
        println("----- $n_elements ELEMENTS -----")
        root, mass, kernel = create_root(n_elements, p_expansion)
        @time t_direct, t_fmm = evaluate_time(root, mass, kernel; run_direct, run_fmm)
        push!(n_element_list, n_elements)
        push!(ts_direct, t_direct)
        push!(ts_fmm, t_fmm)
        @show t_direct t_fmm
    end
end

@testset "derivatives" begin
    kernel = fmm.Gravitational(2) # 2-D

    # build symbolic kernel
    sym.@variables x y
    phi = log(sqrt(x^2 + y^2))

    # get derivatives
    symbolic_derivatives = Array{typeof(phi),2}(undef,6,6)
    symbolic_derivatives[1,1] = phi
    ddx = sym.derivative(phi, x)
    symbolic_derivatives[2,1] = ddx
    ddy = sym.derivative(phi, y)
    symbolic_derivatives[1,2] = ddy
    d2dx2 = sym.derivative(ddx, x)
    symbolic_derivatives[3,1] = d2dx2
    d2dxdy = sym.derivative(ddx, y)
    symbolic_derivatives[2,2] = d2dxdy
    d2dy2 = sym.derivative(ddy, y)
    symbolic_derivatives[1,3] = d2dy2
    d3dx3 = sym.derivative(d2dx2, x)
    symbolic_derivatives[4,1] = d3dx3
    d3dx2dy = sym.derivative(d2dx2, y)
    symbolic_derivatives[3,2] = d3dx2dy
    d3dxdy2 = sym.derivative(d2dxdy, y)
    symbolic_derivatives[2,3] = d3dxdy2
    d3dy3 = sym.derivative(d2dy2, y)
    symbolic_derivatives[1,4] = d3dy3
    d4dx4 = sym.derivative(d3dx3, x)
    symbolic_derivatives[5,1] = d4dx4
    d4dx3dy = sym.derivative(d3dx3, y)
    symbolic_derivatives[4,2] = d4dx3dy
    d4dx2dy2 = sym.derivative(d3dx2dy, y)
    symbolic_derivatives[3,3] = d4dx2dy2
    d4dxdy3 = sym.derivative(d3dxdy2, y)
    symbolic_derivatives[2,4] = d4dxdy3
    d4dy4 = sym.derivative(d3dy3, y)
    symbolic_derivatives[1,5] = d4dy4
    d5dx5 = sym.derivative(d4dx4, x)
    symbolic_derivatives[6,1] = d5dx5
    d5dx4dy = sym.derivative(d4dx4, y)
    symbolic_derivatives[5,2] = d5dx4dy
    d5dx3dy2 = sym.derivative(d4dx3dy, y)
    symbolic_derivatives[4,3] = d5dx3dy2
    d5dx2dy3 = sym.derivative(d4dx2dy2, y)
    symbolic_derivatives[3,4] = d5dx2dy3
    d5dxdy4 = sym.derivative(d4dxdy3, y)
    symbolic_derivatives[2,5] = d5dxdy4
    d5dy5 = sym.derivative(d4dy4, y)
    symbolic_derivatives[1,6] = d5dy5

    # compare
    xy = (rand(2) .- 0.5) * 5
    s = size(kernel.potential_derivatives)[1]
    kernel_values = Array{Float64, 2}(undef,s,s)
    symbolic_values = Array{Float64, 2}(undef,s,s)
    for i = 1:s
        for j = 1:s
            if i+j <= s+1
                kernel_values[i,j] = kernel.potential_derivatives[i,j](xy, 1.0, 1.0)
                r = sym.substitute(symbolic_derivatives[i,j],Dict(x=>xy[1],y=>xy[2])).val
                symbolic_values[i,j] = r
                @test isapprox(kernel_values[i,j], symbolic_values[i,j]; atol=1e-8)
            end
        end
    end

end
