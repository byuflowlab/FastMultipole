function visualize_bodies(name, system, probe_indices=())
    n_bodies = get_n_bodies(system)
    body_locations = zeros(3,n_bodies,1,1)
    body_radii = zeros(n_bodies,1,1)
    scalar_strength = zeros(n_bodies,1,1)
    vector_strength = zeros(3,n_bodies,1,1)
    scalar_potential = zeros(n_bodies,1,1)
    vector_potential = zeros(3,n_bodies,1,1)
    velocity = zeros(3,n_bodies,1,1)
    for i in 1:n_bodies
        body_locations[:,i,1,1] .= system[i,POSITION]
        body_radii[i,1,1] = system[i,RADIUS]
        scalar_strength[i,1,1] = system[i,SCALAR_STRENGTH]
        vector_strength[:,i,1,1] .= system[i,VECTOR_STRENGTH]
        scalar_potential[i,1,1] = system[i,SCALAR_POTENTIAL]
        vector_potential[:,i,1,1] .= system[i,VECTOR_POTENTIAL]
        velocity[:,i,1,1] .= system[i,VELOCITY]
    end
    vtk_grid(name*"_bodies", body_locations) do vtk
        vtk["radius"] = body_radii
        vtk["scalar strength"] = scalar_strength
        vtk["vector strength"] = vector_strength
        vtk["scalar potential"] = scalar_potential
        vtk["vector potential"] = vector_potential
        vtk["velocity"] = velocity
    end

    #####
    ##### probes
    #####
    n_probes = length(probe_indices)
    if n_probes > 0
        body_locations = zeros(3,n_probes,1,1)
        body_radii = zeros(n_probes,1,1)
        scalar_strength = zeros(n_probes,1,1)
        vector_strength = zeros(3,n_probes,1,1)
        for (i_probe,i_body) in enumerate(probe_indices)
            body_locations[:,i_probe,1,1] .= system[i_body,POSITION]
            body_radii[i_probe,1,1] = system[i_body,RADIUS]
            scalar_strength[i_probe,1,1] = system[i_body,SCALAR_STRENGTH]
            vector_strength[:,i_probe,1,1] .= system[i_body,VECTOR_STRENGTH]
        end
        vtk_grid(name*"_probes", body_locations) do vtk
            vtk["radius"] = body_radii
            vtk["scalar strength"] = scalar_strength
            vtk["vector strength"] = vector_strength
            vtk["index"] = probe_indices
        end
    end
end

function visualize_bodies(name, systems::Tuple, probe_indices_list)
    for (system, probe_indices) in zip(systems, probe_indices_list)
        visualize_bodies(name, system, probe_indices)
    end
end

function visualize_bodies(name, systems::Tuple, probe_indices_list::Nothing)
    for system in systems
        visualize_bodies(name, system)
    end
end

function visualize(name, system, tree; probe_indices=nothing, toggle_branches=true, toggle_bodies=false)
    #####
    ##### branches
    #####
    if toggle_branches
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
    end

    #####
    ##### bodies
    #####
    if toggle_bodies
        visualize_bodies(name, system, probe_indices)
    end

    return nothing
end
