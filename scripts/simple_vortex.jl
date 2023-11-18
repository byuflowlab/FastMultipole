include("../test/vortex.jl")
using Random
using WriteVTK

function generate_vortices(seed, n_bodies; radius_factor=0.1)
    Random.seed!(123)
    bodies = rand(7,n_bodies)
    bodies[4,:] ./= (n_bodies^(1/3)*2)
    bodies[4,:] .*= radius_factor
    system = VortexParticles(bodies)
    return system
end

function bm_fmm()
    expansion_order, n_per_branch, theta = 5, 5, 0.4
    n_bodies = 50
    system = generate_vortices(123, n_bodies)
    fmm.fmm!(system; expansion_order=expansion_order, n_per_branch=n_per_branch, theta=theta, nearfield=true, farfield=true, unsort_bodies=true, shrink_recenter=false)
    return system.potential[2:4,:], system.velocity_stretching[1:3,:]
end

function bm_fmm(expansion_order, n_per_branch, theta, n_bodies, shrink_recenter)
    system = generate_vortices(123, n_bodies)
    tree = fmm.fmm!(system; expansion_order=expansion_order, n_per_branch=n_per_branch, theta=theta, nearfield=true, farfield=true, unsort_bodies=true, shrink_recenter=shrink_recenter)
    return system.potential[2:4,:], system.velocity_stretching[1:3,:], system, tree
end

function bm_direct()
    n_bodies = 50
    system = generate_vortices(123, n_bodies)
    fmm.direct!(system, 1:n_bodies, system, 1:n_bodies)
    return system.potential[2:4,:], system.velocity_stretching[1:3,:]
end

function bm_direct(n_bodies)
    system = generate_vortices(123, n_bodies)
    fmm.direct!(system, 1:n_bodies, system, 1:n_bodies)
    return system.potential[2:4,:], system.velocity_stretching[1:3,:], system
end

function bm_fmm_accuracy(expansion_order, n_per_branch, theta, n_bodies, shrink_recenter)
    system = generate_vortices(123, n_bodies)
    tree = fmm.fmm!(system; expansion_order=expansion_order, n_per_branch=n_per_branch, shrink_recenter=shrink_recenter, theta=theta, nearfield=true, farfield=true, unsort_bodies=true)
    system2 = generate_vortices(123, n_bodies)
    fmm.direct!(system2, 1:n_bodies, system2, 1:n_bodies)
    phi = system.potential[2:4,:]
    phi2 = system2.potential[2:4,:]
    v = system.velocity_stretching[1:3,:]
    v2 = system2.velocity_stretching[1:3,:]
    return maximum(abs.(phi2 - phi))/maximum(abs.(phi2)), maximum(abs.(v2 - v))/maximum(abs.(v2)), system, tree, system2
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
    scalar_potential = Array{Float64,3}(undef,n_bodies,1,1)
    vector_potential = Array{Float64,4}(undef,3,n_bodies,1,1)
    for i in 1:n_bodies
        body_locations[:,i,1,1] .= system[i,fmm.POSITION]
        body_radii[i,1,1] = system[i,fmm.RADIUS]
        scalar_strength[i,1,1] = system[i,fmm.SCALAR_STRENGTH]
        vector_strength[:,i,1,1] .= system[i,fmm.VECTOR_STRENGTH]
        scalar_potential[i,1,1] = system[i,fmm.SCALAR_POTENTIAL]
        vector_potential[:,i,1,1] .= system[i,fmm.VECTOR_POTENTIAL]
    end
    vtk_grid(name*"_bodies", body_locations) do vtk
        vtk["radius"] = body_radii
        vtk["scalar strength"] = scalar_strength
        vtk["vector strength"] = vector_strength
        vtk["scalar potential"] = scalar_potential
        vtk["vector potential"] = vector_potential
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

function debug()
    expansion_order, n_per_branch, theta, n_bodies, shrink_recenter = 20, 2, 0.4, 14, false

    potential, v, system, tree = bm_fmm(expansion_order, n_per_branch, theta, n_bodies, shrink_recenter)
    potential_direct, v_direct, system_direct = bm_direct(n_bodies)

    err_potential = zeros(size(potential,2))
    for i in 1:length(err_potential)
        err_potential[i] = max(abs.(potential[:,i] - potential_direct[:,i])...)
    end

    @show err_potential
    visualize_tree("test_vpm_tree_20231103", system, tree; probe_indices=[13,14])
end

# err_potential, err_velocity, sys, tree, sys2 = bm_fmm_accuracy(expansion_order, n_per_branch, theta, n_bodies, shrink_recenter)
# println("relative potential error: $err_potential")
# println("relative velocity error: $err_velocity")

system = generate_vortices(123, 50000)
fmm.get_nearfield_parameters(system)