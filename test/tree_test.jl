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

@testset "unsort/resort" begin

n_bodies = 101
bodies = rand(8,n_bodies)
x_presorted = deepcopy(bodies[1:3,:])
bodies[5:8,2:n_bodies] .= 0.0
system = Gravitational(bodies)
tree, m2l_list, direct_list, switch = FastMultipole.fmm!(system; unsort_bodies=false)
x_sorted = deepcopy(bodies[1:3,:])
FastMultipole.unsort!(system, tree)
x_unsorted = deepcopy(bodies[1:3,:])
FastMultipole.resort!(system, tree)
x_resorted = deepcopy(bodies[1:3,:])

@test isapprox(x_presorted, x_unsorted)
@test isapprox(x_sorted, x_resorted)

end

@testset "leaf index" begin

n_bodies = 101
bodies = rand(8,n_bodies)
system = Gravitational(bodies)

tree = fmm.Tree(system; leaf_size=5)
leaf_index = tree.leaf_index

i_leaf = 1
for branch in tree.branches
    if branch.n_branches == 0
        @test i_leaf == branch.i_leaf
        i_leaf += 1
    end
end

end

