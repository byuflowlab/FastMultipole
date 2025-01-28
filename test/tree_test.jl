@testset "octree creation" begin

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
    center, box = FastMultipole.center_box((elements,))

    # manually
    test_center = [0.65, 0.65, 0.55]
    x_min, x_max = 0.1, 1.2
    y_min, y_max = 0.2, 1.1
    z_min, z_max = 0.2, 0.9
    test_box = [max(x_max-test_center[1], test_center[1]-x_min),
                max(y_max-test_center[2], test_center[2]-y_min),
                max(z_max-test_center[3], test_center[3]-z_min)]

    for i in 1:3
        @test isapprox(center[i], test_center[i]; atol=1e-12)
    end
    for i in 1:3
        @test isapprox(test_box[i], box[i]; atol=1e-12)
    end

    # test branch! function
    expansion_order, leaf_size, multipole_threshold = 2, SVector{1}(1), 0.5
    tree = FastMultipole.Tree((elements,); expansion_order, leaf_size, shrink_recenter=false)

    r1 = max(max(tree.branches[1].target_box[1],tree.branches[1].target_box[2]),tree.branches[1].target_box[3])
    r2 = max(max(tree.branches[2].target_box[1],tree.branches[2].target_box[2]),tree.branches[2].target_box[3])

    test_branches = [
        5 0.65 0.65 0.55 sqrt(3)*r1;
        2 0.65-r1/2 0.65-r1/2 0.55-r1/2 sqrt(3)*r2;
        1 0.65+r1/2 0.65+r1/2 0.55-r1/2 sqrt(3)*r2;
        1 0.65-r1/2 0.65-r1/2 0.55+r1/2 sqrt(3)*r2;
        1 0.65+r1/2 0.65+r1/2 0.55+r1/2 sqrt(3)*r2;
        1 tree.branches[2].target_center[1]-r2/2 tree.branches[2].target_center[2]-r2/2 tree.branches[2].target_center[3]-r2/2 sqrt(3)*r2/2;
        1 tree.branches[2].target_center[1]-r2/2 tree.branches[2].target_center[2]-r2/2 tree.branches[2].target_center[3]+r2/2 sqrt(3)*r2/2
    ]

    @test length(tree.branches) == size(test_branches)[1]

    for i_branch in 1:length(tree.branches)
        @test isapprox(length(tree.branches[i_branch].bodies_index[1]), test_branches[i_branch,1]; atol=1e-8)
        for i in 1:3
            @test isapprox(tree.branches[i_branch].target_center[i], test_branches[i_branch,1+i]; atol=1e-7)
        end
        @test isapprox(tree.branches[i_branch].target_radius, test_branches[i_branch,5]; atol=1e-7)
    end
end

@testset "octree creation: shrinking" begin

    # build list of elements to sort
    xs = [
        1.2 1.1 0.8;
        0.8 0.9 0.2;
        0.1 0.2 0.9;
        0.1 0.3 0.2;
        0.2 0.25 0.4
    ]
    Random.seed!(123)
    radii = rand(size(xs,1))
    ms = rand(size(xs)[1])
    bodies = vcat(xs',radii',ms',zeros(3,length(ms)))
    elements = Gravitational(bodies)

    # test target box
    center, box = FastMultipole.center_box((elements,))
    source_center, source_box = FastMultipole.source_center_box((elements,), FastMultipole.get_bodies_index((elements,)))

    # test source box
    target_radius, source_radius = FastMultipole.shrink_radius(center, source_center, (elements,), SVector{1}((1:5,)))

    # manually
    x_min, x_max, y_min, y_max, z_min, z_max = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    for i in 1:5
        x_min = min(x_min, elements.bodies[i].position[1] - elements.bodies[i].radius)
        y_min = min(y_min, elements.bodies[i].position[2] - elements.bodies[i].radius)
        z_min = min(z_min, elements.bodies[i].position[3] - elements.bodies[i].radius)
        x_max = max(x_max, elements.bodies[i].position[1] + elements.bodies[i].radius)
        y_max = max(y_max, elements.bodies[i].position[2] + elements.bodies[i].radius)
        z_max = max(z_max, elements.bodies[i].position[3] + elements.bodies[i].radius)
    end

    test_source_center = SVector{3}((x_min+x_max)*0.5, (y_min+y_max)*0.5, (z_min+z_max)*0.5)
    test_source_box = [(x_max-x_min)*0.5, (y_max-y_min)*0.5, (z_max-z_min)*0.5]
    test_source_radius = 0.0
    for i_body in 1:5
        test_source_radius = max(test_source_radius, norm(elements.bodies[i_body].position - test_source_center) + elements.bodies[i_body].radius)
    end

    @test isapprox(test_source_center, source_center; atol=1e-12)
    @test isapprox(test_source_radius, source_radius; atol=1e-12)
    @test isapprox(test_source_box, source_box; atol=1e-12)

    # test branch! function
    expansion_order, leaf_size, multipole_threshold = 2, SVector{1}(1), 0.5
    tree = FastMultipole.Tree((elements,); expansion_order, leaf_size, shrink_recenter=true)

    # test branch 3-5 (leaf branches)
    function test_leaf(branches, elements, i_leaf)
        branch = tree.branches[i_leaf]
        i_body = branch.bodies_index[1][1]
        test_target_center = elements.bodies[i_body].position
        test_source_center = elements.bodies[i_body].position
        test_target_radius = 0.0
        test_source_radius = elements.bodies[i_body].radius
        test_target_box = SVector{3}(0.0,0.0,0.0)
        test_source_box = SVector{3}(test_source_radius for _ in 1:3)

        @test isapprox(test_target_center, branch.target_center; atol=1e-12)
        @test isapprox(test_source_center, branch.source_center; atol=1e-12)
        @test isapprox(test_target_radius, branch.target_radius; atol=1e-12)
        @test isapprox(test_source_radius, branch.source_radius; atol=1e-12)
        @test isapprox(test_target_box, branch.target_box; atol=1e-12)
        @test isapprox(test_source_box, branch.source_box; atol=1e-12)
    end

    test_leaf(tree.branches, elements, 3)
    test_leaf(tree.branches, elements, 4)
    test_leaf(tree.branches, elements, 5)
    test_leaf(tree.branches, elements, 6)
    test_leaf(tree.branches, elements, 7)

    # test branch 2 (mid level)
    # NOTE: this test works by comparing to body locations
    #       it only works because the child branches of
    #       branch 2 have a target radius of zero
    #       (exactly 1 body each at their centers)
    branch = tree.branches[2]
    bodies_index = branch.bodies_index[1]
    centers = [elements.bodies[i].position for i in bodies_index]
    test_target_center = sum(centers) / length(centers)

    # source center
    min_x_source = minimum([centers[i][1] - elements.bodies[i].radius for i in bodies_index])
    min_y_source = minimum([centers[i][2] - elements.bodies[i].radius for i in bodies_index])
    min_z_source = minimum([centers[i][3] - elements.bodies[i].radius for i in bodies_index])
    max_x_source = maximum([centers[i][1] + elements.bodies[i].radius for i in bodies_index])
    max_y_source = maximum([centers[i][2] + elements.bodies[i].radius for i in bodies_index])
    max_z_source = maximum([centers[i][3] + elements.bodies[i].radius for i in bodies_index])

    test_source_center = SVector{3}((min_x_source+max_x_source)*0.5, (min_y_source+max_y_source)*0.5, (min_z_source+max_z_source)*0.5)
    test_source_box = SVector{3}((max_x_source-min_x_source)*0.5, (max_y_source-min_y_source)*0.5, (max_z_source-min_z_source)*0.5)

    test_target_radius = maximum([norm(elements.bodies[i].position - branch.target_center) for i in bodies_index])
    test_source_radius = maximum([norm(elements.bodies[i].position - branch.source_center) + elements.bodies[i].radius for i in bodies_index])
    bodies_x = [c[1] for c in centers]
    bodies_y = [c[2] for c in centers]
    bodies_z = [c[3] for c in centers]
    min_x = minimum(bodies_x)
    max_x = maximum(bodies_x)
    min_y = minimum(bodies_y)
    max_y = maximum(bodies_y)
    min_z = minimum(bodies_z)
    max_z = maximum(bodies_z)
    test_target_box = SVector{3}((max_x - min_x)/2, (max_y - min_y)/2, (max_z - min_z)/2)

    @test isapprox(test_target_center, branch.target_center; atol=1e-12)
    @test isapprox(test_source_center, branch.source_center; atol=1e-12)
    @test isapprox(test_target_radius, branch.target_radius; atol=1e-12)
    @test isapprox(test_source_radius, branch.source_radius; atol=1e-12)
    @test isapprox(test_target_box, branch.target_box; atol=1e-12)
    @test isapprox(test_source_box, branch.source_box; atol=1e-12)

end

@testset "octree creation: unsort/resort" begin

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

@testset "octree creation: leaf index" begin

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

