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
    tree = fmm.Tree(system, expansion_order, n_per_branch; shrinking=true)
    fmm.fmm!(tree, system; theta=theta, reset_tree=true, nearfield=true, farfield=true, unsort_bodies=true)
    return nothing
end

function bm_direct()    
    n_bodies = 5000
    system2 = generate_gravitational(123, n_bodies)
    fmm.direct!(system2, 1:n_bodies, system2, 1:n_bodies)
    return nothing
end

function bm_fmm_accuracy(expansion_order, n_per_branch, theta, n_bodies, shrinking)
    system = generate_gravitational(123, n_bodies)
    tree = fmm.Tree(system, expansion_order, n_per_branch; shrinking=shrinking)
    fmm.fmm!(tree, system; theta=theta, reset_tree=true, nearfield=true, farfield=true, unsort_bodies=true)
    println("BEGIN DIRECT")
    system2 = generate_gravitational(123, n_bodies)
    fmm.direct!(system2, 1:n_bodies, system2, 1:n_bodies)
    phi = system.potential[1,:]
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

# note: old_bodies[tree.index_list[1]] = systems[1].bodies
# println("Run FMM:")
# @time bm_fmm()
# @time bm_fmm()

# println("Run Direct:")
# @time bm_direct()
# @time bm_direct()
# @btime fmm.fmm!($tree, $systems, $options; unsort_bodies=true)
println("Calculating accuracy:")
expansion_order, n_per_branch, theta = 13, 1, 0.4
n_bodies = 5000
shrinking = false
accuracy, system, tree, system2 = bm_fmm_accuracy(expansion_order, n_per_branch, theta, n_bodies, shrinking)
println("accuracy: $accuracy")

# visualize tree
visualize_tree("test_fmm", system, tree; probe_indices=[11])
visualize_tree("test_direct", system2, tree)

# run some tests
# system3 = Gravitational(bodies)
# fmm.fmm!(system3; n_per_branch=n_per_branch, theta=0.34, nearfield=true, farfield=true, unsort_bodies=true)
# @time fmm.fmm!(system3; n_per_branch=n_per_branch, theta=0.34, nearfield=true, farfield=true, unsort_bodies=true)
# system4 = Gravitational(bodies)
# fmm.fmm!((system4,); n_per_branch=n_per_branch, theta=0.34, nearfield=true, farfield=true, unsort_bodies=true)
# @time fmm.fmm!((system4,); n_per_branch=n_per_branch, theta=0.34, nearfield=true, farfield=true, unsort_bodies=true)
# system5 = Gravitational(bodies)
# fmm.fmm!(system5, system5; n_per_branch_source=n_per_branch, n_per_branch_target=n_per_branch, theta=0.34, nearfield=true, farfield=true, unsort_source_bodies=true, unsort_target_bodies=true)
# @time fmm.fmm!(system5, system5; n_per_branch_source=n_per_branch, n_per_branch_target=n_per_branch, theta=0.34, nearfield=true, farfield=true, unsort_source_bodies=true, unsort_target_bodies=true)
# system6 = Gravitational(bodies)
# fmm.fmm!((system6,), system6; n_per_branch_source=n_per_branch, n_per_branch_target=n_per_branch, theta=0.34, nearfield=true, farfield=true, unsort_source_bodies=true, unsort_target_bodies=true)
# @time fmm.fmm!((system6,), system6; n_per_branch_source=n_per_branch, n_per_branch_target=n_per_branch, theta=0.34, nearfield=true, farfield=true, unsort_source_bodies=true, unsort_target_bodies=true)
# system7 = Gravitational(bodies)
# fmm.fmm!((system7,), (system7,); n_per_branch_source=n_per_branch, n_per_branch_target=n_per_branch, theta=0.34, nearfield=true, farfield=true, unsort_source_bodies=true, unsort_target_bodies=true)
# @time fmm.fmm!((system7,), (system7,); n_per_branch_source=n_per_branch, n_per_branch_target=n_per_branch, theta=0.34, nearfield=true, farfield=true, unsort_source_bodies=true, unsort_target_bodies=true)
# println("done.")

# n_bodies = 30000
# bodies = rand(8,n_bodies)
# systems = (Gravitational(bodies),)
# # systems = (fmm.SortWrapper(Gravitational(bodies)),)
# options = fmm.Options(13,700,10.0)
# old_bodies = deepcopy(systems[1].bodies)
# @time tree = fmm.Tree(systems, options)
# # note: old_bodies[tree.index_list[1]] = systems[1].bodies
# println("Run FMM:")
# # @btime fmm.fmm!($tree, $systems, $options; unsort_bodies=true)
# @time fmm.fmm!(tree, systems, options; unsort_bodies=true)
# println("done.")
# systems2 = (Gravitational(bodies),)
# println("Run direct:")
# @time fmm.direct!(systems2[1], 1:n_bodies, systems2[1], 1:n_bodies)
# println("done.")
# phi = systems[1].potential[1,:]
# phi2 = systems2[1].potential[1,:]
# @show maximum(abs.(phi2 - phi))
