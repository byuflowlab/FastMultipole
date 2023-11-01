include("../test/gravitational.jl")
using BenchmarkTools
using Random
using WriteVTK

function generate_gravitational(seed, n_bodies; radius_factor=0.1)
    Random.seed!(123)
    bodies = rand(8,n_bodies)
    # bodies[1:3,3] .=  0.811770914672987, 0.15526131946379113, 0.30656077208169424
    # bodies[1:3,3] .=   0.7427186184997012, 0.2351893322824516, 0.3380666354208596
    bodies[4,:] ./= (n_bodies^(1/3)*2)
    bodies[4,:] .*= radius_factor
    system = Gravitational(bodies)
end

function bm_fmm()
    expansion_order, n_per_branch, theta = 5, 100, 0.4
    n_bodies = 5000
    system = generate_gravitational(123, n_bodies)
    # tree = fmm.Tree(system; expansion_order, n_per_branch, shrink_recenter=true)
    fmm.fmm!(system; expansion_order=expansion_order, n_per_branch=n_per_branch, theta=theta, nearfield=true, farfield=true, unsort_bodies=true, shrink_recenter=false)
    return nothing
end

function bm_direct()    
    n_bodies = 5000
    system2 = generate_gravitational(123, n_bodies)
    fmm.direct!(system2, 1:n_bodies, system2, 1:n_bodies)
    return sum(system2.potential)
end

function bm_fmm_accuracy(expansion_order, n_per_branch, theta, n_bodies, shrink_recenter)
    system = (generate_gravitational(123, n_bodies),)
    println("Create tree")
    @time tree = fmm.Tree(system; expansion_order, n_per_branch, shrink_recenter=shrink_recenter)
    println("Run fmm")
    @time fmm.fmm!(tree, system; theta=theta, reset_tree=true, nearfield=true, farfield=true, unsort_bodies=true)
    println("BEGIN DIRECT")
    system2 = generate_gravitational(123, n_bodies)
    @time fmm.direct!(system2, 1:n_bodies, system2, 1:n_bodies)
    phi = system[1].potential[1,:]
    phi2 = system2.potential[1,:]
    return maximum(abs.(phi2 - phi)), system, tree, system2
end

function bm_fmm_accuracy_dual_tree(expansion_order, n_per_branch, theta, n_bodies, shrink_recenter)
    system = (generate_gravitational(123, n_bodies),)
    println("Create trees")
    @time source_tree = fmm.Tree(system; expansion_order, n_per_branch, shrink_recenter=shrink_recenter)
    @time target_tree = fmm.Tree(system; expansion_order, n_per_branch, shrink_recenter=shrink_recenter)
    println("Run fmm")
    @time fmm.fmm!(target_tree, system, source_tree, system; theta=theta, nearfield=true, farfield=true)
    println("BEGIN DIRECT")
    system2 = generate_gravitational(123, n_bodies)
    @time fmm.direct!(system2, 1:n_bodies, system2, 1:n_bodies)
    phi = system[1].potential[1,:]
    phi2 = system2.potential[1,:]
    return maximum(abs.(phi2 - phi)), system, tree, system2
end

function visualize_tree(name, system, tree; probe_indices=[])
    #####
    ##### branches
    #####
    branches = tree.branches
    n_branches = length(branches)
    branch_locations = Array{Float64,4}(undef,3,n_branches,1,1)
    branch_radii = Array{Float64,3}(undef,n_branches,1,1)
    for (i,branch) in enumerate(branches)
        branch_locations[:,i,1,1] .= branch.center
        branch_radii[i,1,1] = branch.radius
    end
    vtk_grid(name*"_branches", branch_locations) do vtk
        vtk["radius"] = branch_radii
    end

    #####
    ##### bodies
    #####
    n_bodies = length(system)
    body_locations = Array{Float64,4}(undef,3,n_bodies,1,1)
    body_radii = Array{Float64,3}(undef,n_bodies,1,1)
    scalar_strength = Array{Float64,3}(undef,n_bodies,1,1)
    vector_strength = Array{Float64,4}(undef,3,n_bodies,1,1)
    for i in 1:n_bodies
        body_locations[:,i,1,1] .= system[i,fmm.POSITION]
        body_radii[i,1,1] = system[i,fmm.RADIUS]
        scalar_strength[i,1,1] = system[i,fmm.SCALAR_STRENGTH]
        vector_strength[:,i,1,1] .= system[i,fmm.VECTOR_STRENGTH]
    end
    vtk_grid(name*"_bodies", body_locations) do vtk
        vtk["radius"] = body_radii
        vtk["scalar strength"] = scalar_strength
        vtk["vector strength"] = vector_strength
    end
    
    #####
    ##### probes
    #####
    n_probes = length(probe_indices)
    if n_probes > 0
        body_locations = Array{Float64,4}(undef,3,n_probes,1,1)
        body_radii = Array{Float64,3}(undef,n_probes,1,1)
        scalar_strength = Array{Float64,3}(undef,n_probes,1,1)
        vector_strength = Array{Float64,4}(undef,3,n_probes,1,1)
        for (i_probe,i_body) in enumerate(probe_indices)
            body_locations[:,i_probe,1,1] .= system[i_body,fmm.POSITION]
            body_radii[i_probe,1,1] = system[i_body,fmm.RADIUS]
            scalar_strength[i_probe,1,1] = system[i_body,fmm.SCALAR_STRENGTH]
            vector_strength[:,i_probe,1,1] .= system[i_body,fmm.VECTOR_STRENGTH]
        end
        vtk_grid(name*"_probes", body_locations) do vtk
            vtk["radius"] = body_radii
            vtk["scalar strength"] = scalar_strength
            vtk["vector strength"] = vector_strength
            vtk["index"] = probe_indices
        end
    end

    return nothing
end
# println("generate bodies")
# @time generate_gravitational(123, 5000)
# @time generate_gravitational(123, 5000)

# note: old_bodies[tree.index_list[1]] = systems[1].bodies
println("Run FMM:")
@time bm_fmm()
@time bm_fmm()

println("Run Direct:")
@time bm_direct()
@time bm_direct()
# @btime fmm.fmm!($tree, $systems, $options; unsort_bodies=true)
# println("Calculating accuracy:")
expansion_order, n_per_branch, theta = 8, 100, 0.31
n_bodies = 5000
shrink_recenter, ndivisions = true, 5
# sys = generate_gravitational(123,n_bodies)
# function bmtree()# let sys=sys, expansion_order=expansion_order, n_per_branch=n_per_branch, ndivisions=ndivisions, shrink_recenter=shrink_recenter
#         return fmm.Tree(sys; expansion_order, n_per_branch, ndivisions=ndivisions, shrink_recenter=shrink_recenter)
#     # end
# end
# @time fmm.Tree(sys; expansion_order=expansion_order, n_per_branch=n_per_branch, ndivisions=ndivisions, shrink_recenter=shrink_recenter)
# fmm.unsort!(sys, tree)
# @time bmtree()
# fmm.unsort!(sys, tree)
# @time bmtree()
# @time bmtree()
# @time bmtree()

# println("done")

# sys_noshrinking = generate_gravitational(123,n_bodies)
# tree_noshrinking = fmm.Tree(sys_noshrinking; expansion_order, n_per_branch, ndivisions=5, shrink_recenter=false)

# println("done")
run_bm_accuracy() = bm_fmm_accuracy(expansion_order, n_per_branch, theta, n_bodies, shrink_recenter)
accuracy, system, tree, system2 = run_bm_accuracy()
accuracy, system, tree, system2 = run_bm_accuracy()
println("single tree accuracy: $accuracy")
run_bm_accuracy_dual_tree() = bm_fmm_accuracy_dual_tree(expansion_order, n_per_branch, theta, n_bodies, shrink_recenter)
accuracy, system, tree, system2 = run_bm_accuracy_dual_tree()
accuracy, system, tree, system2 = run_bm_accuracy_dual_tree()
println("dual tree accuracy: $accuracy")

# visualize tree
# visualize_tree("test_fmm", system, tree; probe_indices=[11])
# visualize_tree("test_direct", system2, tree)
# visualize_tree("test_shrinking", sys, tree)
# visualize_tree("test_noshrinking", sys_noshrinking, tree_noshrinking)

# test various single/dual tree, single/multi branch
n_bodies = 5000
seed = 123
validation_system = generate_gravitational(seed, n_bodies; radius_factor=0.1)
fmm.direct!(validation_system)
validation_potential = validation_system.potential[1,:]

system3 = generate_gravitational(seed, n_bodies; radius_factor=0.1)
fmm.fmm!(system3; n_per_branch=n_per_branch, theta=0.34, nearfield=true, farfield=true, unsort_bodies=true)
potential3 = system3.potential[1,:]
println("Case 3 err:")
@show maximum(potential3 - validation_potential)
system4 = generate_gravitational(seed, n_bodies; radius_factor=0.1)
fmm.fmm!((system4,); n_per_branch=n_per_branch, theta=0.34, nearfield=true, farfield=true, unsort_bodies=true)
potential4 = system4.potential[1,:]
println("Case 4 err:")
@show maximum(potential4 - validation_potential)
system5 = generate_gravitational(seed, n_bodies; radius_factor=0.1)
fmm.fmm!(system5, system5; n_per_branch_source=n_per_branch, n_per_branch_target=n_per_branch, theta=0.34, nearfield=true, farfield=true, unsort_source_bodies=true, unsort_target_bodies=true)
potential5 = system5.potential[1,:]
println("Case 5 err:")
@show maximum(potential5 - validation_potential)
system6 = generate_gravitational(seed, n_bodies; radius_factor=0.1)
fmm.fmm!((system6,), system6; n_per_branch_source=n_per_branch, n_per_branch_target=n_per_branch, theta=0.34, nearfield=true, farfield=true, unsort_source_bodies=true, unsort_target_bodies=true)
potential6 = system6.potential[1,:]
println("Case 6 err:")
@show maximum(potential6 - validation_potential)
system7 = generate_gravitational(seed, n_bodies; radius_factor=0.1)
fmm.fmm!((system7,), (system7,); n_per_branch_source=n_per_branch, n_per_branch_target=n_per_branch, theta=0.34, nearfield=true, farfield=true, unsort_source_bodies=true, unsort_target_bodies=true)
potential7 = system7.potential[1,:]
println("Case 7 err:")
@show maximum(potential7 - validation_potential)

# test SortWrapper
system8 = fmm.SortWrapper(generate_gravitational(seed, n_bodies; radius_factor=0.1))
fmm.fmm!(system8; n_per_branch=n_per_branch, theta=0.34, nearfield=true, farfield=true, unsort_bodies=true)
potential8 = system8.system.potential[1,:]
println("Case 8 err:")
@show maximum(potential8 - validation_potential)
system9 = fmm.SortWrapper(generate_gravitational(seed, n_bodies; radius_factor=0.1))
fmm.fmm!((system9,); n_per_branch=n_per_branch, theta=0.34, nearfield=true, farfield=true, unsort_bodies=true)
potential9 = system9.system.potential[1,:]
println("Case 9 err:")
@show maximum(potential9 - validation_potential)
system10 = fmm.SortWrapper(generate_gravitational(seed, n_bodies; radius_factor=0.1))
fmm.fmm!(system10, system10; n_per_branch_source=n_per_branch, n_per_branch_target=n_per_branch, theta=0.34, nearfield=true, farfield=true, unsort_source_bodies=true, unsort_target_bodies=true)
potential10 = system10.system.potential[1,:]
println("Case 10 err:")
@show maximum(potential10 - validation_potential)
system11 = fmm.SortWrapper(generate_gravitational(seed, n_bodies; radius_factor=0.1))
fmm.fmm!((system11,), system11; n_per_branch_source=n_per_branch, n_per_branch_target=n_per_branch, theta=0.34, nearfield=true, farfield=true, unsort_source_bodies=true, unsort_target_bodies=true)
potential11 = system11.system.potential[1,:]
println("Case 11 err:")
@show maximum(potential11 - validation_potential)
system12 = fmm.SortWrapper(generate_gravitational(seed, n_bodies; radius_factor=0.1))
fmm.fmm!((system12,), (system12,); n_per_branch_source=n_per_branch, n_per_branch_target=n_per_branch, theta=0.34, nearfield=true, farfield=true, unsort_source_bodies=true, unsort_target_bodies=true)
potential12 = system12.system.potential[1,:]
println("Case 12 err:")
@show maximum(potential12 - validation_potential)
println("done.")
