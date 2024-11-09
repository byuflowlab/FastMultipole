#------- direct interactions -------#

function nearfield_singlethread!(target_system, target_branches, derivatives_switch, source_systems, source_branches::Vector{<:MultiBranch}, direct_list)
    for (i_source_system, source_system) in enumerate(source_systems)
        nearfield_singlethread!(target_system, target_branches, derivatives_switch, source_system, i_source_system, source_branches, direct_list)
    end
end

function nearfield_singlethread!(target_system, target_branches, derivatives_switch, source_system, i_source_system, source_branches, direct_list)
    for (i_target, i_source) in direct_list
        target_index = target_branches[i_target].bodies_index
        source_index = source_branches[i_source].bodies_index[i_source_system]
        nearfield_singlethread!(target_system, target_index, derivatives_switch, source_system, source_index)
    end
end

function nearfield_singlethread!(target_system, target_branches, derivatives_switch, source_system, source_branches::Vector{<:SingleBranch}, direct_list)
    for (i_target, i_source) in direct_list
        target_index = target_branches[i_target].bodies_index
        source_index = source_branches[i_source].bodies_index
        nearfield_singlethread!(target_system, target_index, derivatives_switch, source_system, source_index)
    end
end

function nearfield_singlethread!(target_system, target_index, derivatives_switch, source_system, source_index)
    _direct!(target_system, target_index, derivatives_switch, source_system, source_index)
end

function nearfield_singlethread!(target_systems, target_indices::SVector, derivatives_switches::Tuple, source_system, source_index)
    for (target_system, target_index, derivatives_switch) in zip(target_systems, target_indices, derivatives_switches)
        _direct!(target_system, target_index, derivatives_switch, source_system, source_index)
    end
end

function get_n_interactions(target_systems, target_branches, source_system, i_source_system, source_branches::Vector{<:MultiBranch}, direct_list)
    n_interactions = 0
    for (i_target, i_source) in direct_list
        n_interactions += get_n_bodies(target_branches[i_target]) * length(source_branches[i_source].bodies_index[i_source_system])
    end

    return n_interactions
end

function get_n_interactions(target_systems, target_branches, source_system, source_branches::Vector{<:SingleBranch}, direct_list)
    n_interactions = 0
    for (i_target, i_source) in direct_list
        n_interactions += get_n_bodies(target_branches[i_target]) * length(source_branches[i_source].bodies_index)
    end

    return n_interactions
end

function make_assignments!(assignments, target_branches, i_source_system, source_branches, direct_list, n_threads, n_per_thread)
    i_start = 1
    i_end = 1
    i_thread = 1
    n_interactions = 0

    # loop over interaction list
    for (i_target, i_source) in direct_list
        # update number of interactions in the current assignment
        n_interactions += get_n_bodies(target_branches[i_target]) * length(source_branches[i_source].bodies_index[i_source_system])

        # if we exceed n_per_thread, finish this assignment and reset counters
        if n_interactions >= n_per_thread
            assignments[i_thread] = i_start:i_end
            i_start = i_end + 1
            i_thread += 1
            n_interactions = 0
        end

        i_end += 1
    end

    i_thread <= n_threads && (assignments[i_thread] = i_start:length(direct_list))
end

function make_assignments!(assignments, target_branches, source_branches, direct_list, n_threads, n_per_thread)
    i_start = 1
    i_end = 1
    i_thread = 1
    n_interactions = 0

    # loop over interaction list
    for (i_target, i_source) in direct_list

        # update number of interactions in the current assignment
        n_interactions += get_n_bodies(target_branches[i_target]) * length(source_branches[i_source].bodies_index)

        # if we exceed n_per_thread, finish this assignment and reset counters
        if n_interactions >= n_per_thread
            assignments[i_thread] = i_start:i_end
            i_start = i_end + 1
            i_thread += 1
            n_interactions = 0
        end

        i_end += 1
    end

    i_thread <= n_threads && (assignments[i_thread] = i_start:length(direct_list))
end

function execute_assignment!(target_systems, target_branches, derivatives_switches, source_system, i_source_system, source_branches, direct_list, assignment)
    for i_interaction in assignment
        i_target, i_source = direct_list[i_interaction]
        nearfield_singlethread!(target_systems, target_branches[i_target].bodies_index, derivatives_switches, source_system, source_branches[i_source].bodies_index[i_source_system])
    end
end

function execute_assignment!(target_systems, target_branches, derivatives_switches, source_system, source_branches, direct_list, assignment)
    for i_interaction in assignment
        i_target, i_source = direct_list[i_interaction]
        nearfield_singlethread!(target_systems, target_branches[i_target].bodies_index, derivatives_switches, source_system, source_branches[i_source].bodies_index)
    end
end

function nearfield_multithread!(target_system, target_branches, derivatives_switch, source_system, source_branches, direct_list, n_threads)
    if n_threads == 1 # single thread

        nearfield_singlethread!(target_system, target_branches, derivatives_switch, source_system, source_branches, direct_list)

    else # multithread

        _nearfield_multithread!(target_system, target_branches, derivatives_switch, source_system, source_branches, direct_list, n_threads)

    end

    return nothing
end

function _nearfield_multithread!(target_system, target_branches, derivatives_switch, source_systems::Tuple, source_branches, direct_list, n_threads)
    for (i_source_system, source_system) in enumerate(source_systems)
        _nearfield_multithread!(target_system, target_branches, derivatives_switch, source_system, i_source_system, source_branches, direct_list, n_threads)
    end
end

function _nearfield_multithread!(target_systems, target_branches, derivatives_switches, source_system, i_source_system, source_branches, direct_list, n_threads)

    #--- load balance ---#

    # total number of interactions
    n_interactions = get_n_interactions(target_systems, target_branches, source_system, i_source_system, source_branches, direct_list)

    # interactions per thread
    n_per_thread, rem = divrem(n_interactions, n_threads)
    rem > 0 && (n_per_thread += 1)

    # if there are too many threads, we'll actually hurt performance
    n_per_thread < MIN_NPT_NF && (n_per_thread = MIN_NPT_NF)

    # create assignments
    assignments = Vector{UnitRange{Int64}}(undef,n_threads)
    for i in eachindex(assignments)
        assignments[i] = 1:0
    end
    make_assignments!(assignments, target_branches, i_source_system, source_branches, direct_list, n_threads, n_per_thread)

    # execute tasks
    Threads.@threads for assignment in assignments
        execute_assignment!(target_systems, target_branches, derivatives_switches, source_system, i_source_system, source_branches, direct_list, assignment)
    end

end

function _nearfield_multithread!(target_systems, target_branches, derivatives_switches, source_system, source_branches, direct_list, n_threads)

    #--- load balance ---#

    # total number of interactions
    n_interactions = get_n_interactions(target_systems, target_branches, source_system, source_branches, direct_list)

    # interactions per thread
    n_per_thread, rem = divrem(n_interactions, n_threads)
    rem > 0 && (n_per_thread += 1)

    # if there are too many threads, we'll actually hurt performance
    n_per_thread < MIN_NPT_NF && (n_per_thread = MIN_NPT_NF)

    # create assignments
    assignments = Vector{UnitRange{Int64}}(undef,n_threads)
    for i in eachindex(assignments)
        assignments[i] = 1:0
    end
    make_assignments!(assignments, target_branches, source_branches, direct_list, n_threads, n_per_thread)

    # execute tasks
    Threads.@threads for assignment in assignments
        execute_assignment!(target_systems, target_branches, derivatives_switches, source_system, source_branches, direct_list, assignment)
    end

end

#------- UPWARD PASS -------#

function upward_pass_singlethread!(branches::AbstractVector{<:Branch{TF}}, systems, expansion_order::Val{P}, lamb_helmholtz, leaf_index) where {TF,P}

    # try preallocating one container to be reused
    Ts = zeros(TF, length_Ts(P))
    eimϕs = zeros(TF, 2, P+1)
    weights_tmp_1 = initialize_expansion(P, TF)
    weights_tmp_2 = initialize_expansion(P, TF)

    # multipole expansions
    for i_branch in leaf_index
        branch = branches[i_branch]
        body_to_multipole!(branch, systems, branch.harmonics, expansion_order)
    end

    # if predicting error, accumulate charge TODO


    # loop over branches
    for i_branch in length(branches):-1:1 # no need to create a multipole expansion at the very top level
        branch = branches[i_branch]

        if branch.n_branches !== 0 # branch is not a leaf
            # iterate over children
            for i_child in branch.branch_index
                child_branch = branches[i_child]
                multipole_to_multipole!(branch, child_branch, weights_tmp_1, weights_tmp_2, Ts, eimϕs, ζs_mag, Hs_π2, expansion_order, lamb_helmholtz)
            end
        end
    end
end

function body_to_multipole_multithread!(branches::Vector{<:MultiBranch}, systems::Tuple, expansion_order::Val{P}, leaf_index, n_threads) where P
    ## load balance
    leaf_assignments = fill(1:0, length(systems), n_threads)
    for (i_system,system) in enumerate(systems)

        # total number of bodies
        n_bodies = 0
        for i_leaf in leaf_index
            n_bodies += length(branches[i_leaf].bodies_index[i_system])
        end

        # number of bodies per thread
        n_per_thread, rem = divrem(n_bodies, n_threads)
        rem > 0 && (n_per_thread += 1)

        # if there are too many threads, we'll actually hurt performance
        n_per_thread < MIN_NPT_B2M && (n_per_thread = MIN_NPT_B2M)

        # create chunks
        i_start = 1
        i_thread = 1
        n_bodies = 0
        for (i_end,i_leaf) in enumerate(leaf_index)
            n_bodies += length(branches[i_leaf].bodies_index[i_system])
            if n_bodies >= n_per_thread
                leaf_assignments[i_system,i_thread] = i_start:i_end
                i_start = i_end+1
                i_thread += 1
                n_bodies = 0
            end
        end
        i_thread <= n_threads && (leaf_assignments[i_system,i_thread] = i_start:length(leaf_index))
    end

    ## compute multipole expansion coefficients
    Threads.@threads for i_thread in 1:n_threads
        for (i_system,system) in enumerate(systems)
            leaf_assignment = leaf_assignments[i_system,i_thread]
            for i_task in leaf_assignment
                branch = branches[leaf_index[i_task]]
                Threads.lock(branch.lock) do
                    body_to_multipole!(system, branch, branch.bodies_index[i_system], branch.harmonics, expansion_order)
                end
            end
        end
    end
end

function body_to_multipole_multithread!(branches::Vector{<:SingleBranch}, system, expansion_order::Val{P}, leaf_index, n_threads) where P
    ## load balance
    leaf_assignments = fill(1:0, n_threads)

    # total number of bodies
    n_bodies = get_n_bodies(system)

    # number of bodies per thread
    n_per_thread, rem = divrem(n_bodies, n_threads)
    rem > 0 && (n_per_thread += 1)

    # if there are too many threads, we'll actually hurt performance
    n_per_thread < MIN_NPT_B2M && (n_per_thread = MIN_NPT_B2M)

    # create chunks
    i_start = 1
    i_thread = 1
    n_bodies = 0
    for (i_end,i_leaf) in enumerate(leaf_index)
        n_bodies += get_n_bodies(branches[i_leaf])
        if n_bodies >= n_per_thread
            leaf_assignments[i_thread] = i_start:i_end
            i_start = i_end+1
            i_thread += 1
            n_bodies = 0
        end
    end
    i_thread <= n_threads && (leaf_assignments[i_thread] = i_start:length(leaf_index))

    ## compute multipole expansion coefficients
    Threads.@threads for i_thread in 1:n_threads
        for i_assignment in leaf_assignments[i_thread]
            branch = branches[leaf_index[i_assignment]]
            body_to_multipole!(system, branch, branch.bodies_index, branch.harmonics, expansion_order)
        end
    end
end

function translate_multipoles_multithread!(branches::Vector{<:Branch{TF}}, expansion_order::Val{P}, lamb_helmholtz, levels_index, n_threads) where {TF,P}

    # try preallocating one set of containers to be reused
    Ts = [zeros(TF, length_Ts(P)) for _ in 1:n_threads]
    eimϕs = [zeros(TF, 2, P+1) for _ in 1:n_threads]
    weights_tmp_1 = [initialize_expansion(P, TF) for _ in 1:n_threads]
    weights_tmp_2 = [initialize_expansion(P, TF) for _ in 1:n_threads]

    # iterate over levels
    for i_level in length(levels_index):-1:2
        level_index = levels_index[i_level]

        # load balance
        n_branches = length(level_index)
        n_per_thread, rem = divrem(n_branches, n_threads)
        rem > 0 && (n_per_thread += 1)

        # if there are too many threads, we'll actually hurt performance
        n_per_thread < MIN_NPT_M2M && (n_per_thread = MIN_NPT_M2M)

        # assign thread start branches
        i_starts = 1:n_per_thread:n_branches
        Threads.@threads for i_task in 1:length(i_starts)
            # get first branch
            i_start = i_starts[i_task]

            # get final branch
            i_end = min(n_branches, i_start+n_per_thread-1)

            # loop over branches
            for i in i_start:i_end
                i_branch = level_index[i]
                child_branch = branches[i_branch]
                parent_branch = branches[child_branch.i_parent]
                Threads.lock(parent_branch.lock) do
                    multipole_to_multipole!(parent_branch, child_branch, weights_tmp_1[i_task], weights_tmp_2[i_task], Ts[i_task], eimϕs[i_task], ζs_mag, Hs_π2, expansion_order, lamb_helmholtz)
                end
            end
        end
    end
end

function upward_pass_multithread!(branches, systems, expansion_order, lamb_helmholtz, levels_index, leaf_index, n_threads)

    if n_threads == 1

        # single threaded version
        upward_pass_singlethread!(branches, systems, expansion_order, lamb_helmholtz, leaf_index)

    else
        # create multipole expansions
        body_to_multipole_multithread!(branches, systems, expansion_order, leaf_index, n_threads)

        # m2m translation
        translate_multipoles_multithread!(branches, expansion_order, lamb_helmholtz, levels_index, n_threads)
    end
end

#------- direct interaction matrix -------#

# TODO: add influence matrix approach to direct interactions

#------- horizontal pass -------#

function horizontal_pass_singlethread!(target_branches::Vector{<:Branch{TF}}, source_branches, m2l_list, lamb_helmholtz, Pmax::Val{P}) where {TF,P}

    # preallocate containers to be reused
    weights_tmp_1 = initialize_expansion(P, TF)
    weights_tmp_2 = initialize_expansion(P, TF)
    Ts = zeros(TF, length_Ts(P))
    eimϕs = zeros(TF, 2, P+1)

    for (i_target, j_source, expansion_order) in m2l_list
        target_branch = target_branches[i_target]
        source_branch = source_branches[j_source]
        #weights_tmp_1 = target_branch.expansion_storage
        #weights_tmp_2 = source_branch.expansion_storage
        multipole_to_local!(target_branch, source_branch, weights_tmp_1, weights_tmp_2, Ts, eimϕs, ζs_mag, ηs_mag, Hs_π2, expansion_order, lamb_helmholtz)
    end

end

function horizontal_pass_multithread!(target_branches, source_branches::Vector{<:Branch{TF}}, m2l_list, lamb_helmholtz, Pmax::Val{P}, n_threads) where {TF,P}

    if n_threads == 1 # single thread

        horizontal_pass_singlethread!(target_branches, source_branches, m2l_list, lamb_helmholtz, Pmax)

    else # multithread

        # number of translations per thread
        n_per_thread, rem = divrem(length(m2l_list),n_threads)
        rem > 0 && (n_per_thread += 1)
        assignments = 1:n_per_thread:length(m2l_list)

        # try preallocating one set of containers to be reused
        Ts = [zeros(TF, length_Ts(P)) for _ in 1:n_threads]
        eimϕs = [zeros(TF, 2, P+1) for _ in 1:n_threads]
        weights_tmp_1 = [initialize_expansion(P, TF) for _ in 1:n_threads]
        weights_tmp_2 = [initialize_expansion(P, TF) for _ in 1:n_threads]

        # execute tasks
        Threads.@threads for i_thread in 1:length(assignments)
            i_start = assignments[i_thread]
            i_stop = min(i_start+n_per_thread-1, length(m2l_list))
            for (i_target, j_source, expansion_order) in m2l_list[i_start:i_stop]
                Threads.lock(target_branches[i_target].lock) do
                    multipole_to_local!(target_branches[i_target], source_branches[j_source], weights_tmp_1[i_thread], weights_tmp_2[i_thread], Ts[i_thread], eimϕs[i_thread], ζs_mag, ηs_mag, Hs_π2, expansion_order, lamb_helmholtz)
                end
            end
        end
    end

    return nothing
end

#------- DOWNWARD PASS -------#

function downward_pass_singlethread!(branches::AbstractVector{<:Branch{TF}}, systems, expansion_order::Val{P}, lamb_helmholtz, derivatives_switches) where {TF,P}
    # try preallocating one container to be reused
    Ts = zeros(TF, length_Ts(P))
    eimϕs = zeros(TF, 2, P+1)
    weights_tmp_1 = initialize_expansion(P, TF)
    weights_tmp_2 = initialize_expansion(P, TF)
    velocity_n_m = zeros(eltype(branches[1]), 2, 3, size(weights_tmp_1, 3))

    # loop over branches
    for branch in branches
        if branch.n_branches == 0 # leaf level
            evaluate_local!(systems, branch, branch.harmonics, velocity_n_m, expansion_order, lamb_helmholtz, derivatives_switches)
        else
            for i_child_branch in branch.branch_index
                local_to_local!(branches[i_child_branch], branch, weights_tmp_1, weights_tmp_2, Ts, eimϕs, ηs_mag, Hs_π2, expansion_order, lamb_helmholtz)
            end
        end
    end
end

function translate_locals_multithread!(branches::Vector{<:Branch{TF}}, expansion_order::Val{P}, lamb_helmholtz, levels_index, n_threads) where {TF,P}

    # try preallocating one set of containers to be reused
    Ts = [zeros(TF, length_Ts(P)) for _ in 1:n_threads]
    eimϕs = [zeros(TF, 2, P+1) for _ in 1:n_threads]
    weights_tmp_1 = [initialize_expansion(P, TF) for _ in 1:n_threads]
    weights_tmp_2 = [initialize_expansion(P, TF) for _ in 1:n_threads]

    # iterate over levels
    for i_level in 2:length(levels_index)
        level_index = levels_index[i_level]

        # divide chunks
        n_per_thread, rem = divrem(length(level_index),n_threads)
        rem > 0 && (n_per_thread += 1)

        # if there are too many threads, we'll actually hurt performance
        n_per_thread < MIN_NPT_L2L && (n_per_thread = MIN_NPT_L2L)

        # loop over branches
        i_starts = 1:n_per_thread:length(level_index)
        #Threads.@threads for (i_task,i_start) in enumerate(1:n_per_thread:length(level_index))
        Threads.@threads for i_task in 1:length(i_starts)
            i_start = i_starts[i_task]
            i_stop = min(i_start+n_per_thread-1,length(level_index))

            # loop over branches
            for i_child in level_index[i_start]:level_index[i_stop]
                child_branch = branches[i_child]
                local_to_local!(child_branch, branches[child_branch.i_parent], weights_tmp_1[i_task], weights_tmp_2[i_task], Ts[i_task], eimϕs[i_task], ηs_mag, Hs_π2, expansion_order, lamb_helmholtz)
            end
        end
    end
    return nothing
end

function local_to_body_multithread!(branches::Vector{<:Branch{TF}}, systems, derivatives_switches, expansion_order::Val{P}, lamb_helmholtz, leaf_index, n_threads) where {TF,P}
    # try preallocating one set of containers to be reused
    velocity_n_m = [initialize_velocity_n_m(P, TF) for _ in 1:n_threads]

    # create assignments
    n_bodies = 0
    for i_leaf in leaf_index
        n_bodies += get_n_bodies(branches[i_leaf])
    end

    n_per_thread, rem = divrem(n_bodies,n_threads)
    rem > 0 && (n_per_thread += 1)

    assignments = fill(1:0,n_threads)
    i_start = 1
    i_thread = 0
    n_bodies = 0
    for (i_end,i_leaf) in enumerate(leaf_index)
        n_bodies += get_n_bodies(branches[i_leaf])
        if n_bodies >= n_per_thread
            i_thread += 1
            assignments[i_thread] = i_start:i_end
            i_start = i_end+1
            n_bodies = 0
        end
    end
    if i_start <= length(leaf_index)
        i_thread += 1
        assignments[i_thread] = i_start:length(leaf_index)
    end
    resize!(assignments, i_thread)

    # spread remainder across rem chunks
    Threads.@threads for i_thread in eachindex(assignments)
        assignment = assignments[i_thread]
        for i_task in assignment
            leaf = branches[leaf_index[i_task]]
            evaluate_local!(systems, leaf, leaf.harmonics, velocity_n_m[i_thread], expansion_order, lamb_helmholtz, derivatives_switches)
        end
    end
end

function downward_pass_multithread!(branches, systems, derivatives_switch, expansion_order, lamb_helmholtz, levels_index, leaf_index, n_threads)
    if n_threads == 1

        # single thread
        downward_pass_singlethread!(branches, systems, expansion_order, lamb_helmholtz, derivatives_switch)

    else # multithread

        # m2m translation
	    translate_locals_multithread!(branches, expansion_order, lamb_helmholtz, levels_index, n_threads)

        # local to body interaction
        local_to_body_multithread!(branches, systems, derivatives_switch, expansion_order, lamb_helmholtz, leaf_index, n_threads)

    end
end

#--- running FMM ---#

"""
    fmm!(target_systems, source_systems; kwargs...)

Apply all interactions of `source_systems` acting on `target_systems` using the fast multipole method. Assumes compatibility functions have been overloaded for both source and target systems.

# Arguments

- `target_systems`: either

    - a system object for which compatibility functions have been overloaded, or - a tuple of system objects for which compatibility functions have been overloaded

- `source_systems`: either

    - a system object for which compatibility functions have been overloaded, or
    - a tuple of system objects for which compatibility functions have been overloaded

# Optional Arguments

- `expansion_order::Int`: the expansion order to be used
- `leaf_size_source::Int`: maximum number of bodies from `source_systems` allowed in a leaf-level branch
- `leaf_size_target::Int`: maximum number of bodies from `target_systems` allowed in a leaf-level branch
- `multipole_threshold::Float64`: number between 0 and 1 (often denoted theta in [0,1]) controls the accuracy by determining the non-dimensional distance after which multipoles are used; 0 means an infinite distance (no error, high cost), and 1 means barely convergent (high error, low cost)
- `lamb_helmholtz::Bool`: determines whether or not to calculate the induced velocity due to a vector potential using the Lamb-Helmholtz decomposition; erroroneous velocity and gradient will result if `lamb_helmholtz==false` and a vector potential is used.
- `scalar_potential::Bool`: either a `::Bool` or a `::AbstractVector{Bool}` of length `length(target_systems)` indicating whether each system should receive a scalar potential from `source_systems`
- `velocity::Bool`: either a `::Bool` or a `::AbstractVector{Bool}` of length `length(target_systems)` indicating whether each system should receive a velocity from `source_systems`
- `velocity_gradient::Bool`: either a `::Bool` or a `::AbstractVector{Bool}` of length `length(target_systems)` indicating whether each system should receive a velocity gradient from `source_systems`
- `upward_pass::Bool`: whether or not to form the multipole expansions from source bodies and translate them upward in the source tree
- `horizontal_pass::Bool`: whether or not to transform multipole expansions from the source tree into local expansions in the target tree
- `downward_pass::Bool`: whether or not to translate local expansions down to the leaf level of the target tree and evaluate them
- `nearfield::Bool`: indicates whether near-field (comuted without multipoles) interactions should be included
- `farfield::Bool`: indicates whether far-field (comuted with multipoles) interactions should be included
- `self_induced::Bool`: indicates whether to include the interactions of each leaf-level branch on itself
- `unsort_source_bodies::Bool`: indicates whether or not to undo the sort operation used to generate the octree for the `source_systems`
- `unsort_target_bodies::Bool`: indicates whether or not to undo the sort operation used to generate the octree for the `target_systems`
- `source_shink_recenter::Bool`: indicates whether or not to resize branches for the `source_systems` octree after it is created to increase computational efficiency
- `target_shink_recenter::Bool`: indicates whether or not to resize branches for the `target_systems` octree after it is created to increase computational efficiency
- `save_tree_source::Bool`: indicates whether or not to save a VTK file for visualizing the source octree
- `save_tree_target::Bool`: indicates whether or not to save a VTK file for visualizing the target octree
- `save_name_source::String`: name and path of the source octree visualization if `save_tree == true`
- `save_name_target::String`: name and path of the target octree visualization if `save_tree == true`

"""
function fmm!(target_systems, source_systems;
    error_method::ErrorMethod=UnequalBoxes(), predict_error::Val{PE}=Val(false),
    expansion_order::Union{Int,Dynamic{<:Any,<:Any}}=5, leaf_size_source=50, leaf_size_target=50, multipole_threshold=0.4,
    lamb_helmholtz::Bool=false,
    scalar_potential=true, velocity=true, velocity_gradient=true,
    upward_pass::Bool=true, horizontal_pass::Bool=true, downward_pass::Bool=true,
    nearfield::Bool=true, farfield::Bool=true, self_induced::Bool=true,
    unsort_source_bodies=true, unsort_target_bodies=true,
    source_shrink_recenter::Bool=false, target_shrink_recenter::Bool=false,
    save_tree_source=false, save_tree_target=false, save_name_source="source_tree", save_name_target="target_tree",
    nearfield_user::Bool=false
) where PE
    # check for duplicate systems
    target_systems = wrap_duplicates(target_systems, source_systems)

    # create trees
    source_tree = Tree(source_systems; expansion_order=get_Pmax(expansion_order), leaf_size=leaf_size_source, shrink_recenter=source_shrink_recenter)
    target_tree = Tree(target_systems; expansion_order=get_Pmax(expansion_order), leaf_size=leaf_size_target, shrink_recenter=target_shrink_recenter)

    # perform fmm
    m2l_list, direct_list, derivatives_switches = fmm!(target_tree, target_systems, source_tree, source_systems;
        expansion_order, error_method, predict_error,
        lamb_helmholtz,
        scalar_potential, velocity, velocity_gradient,
        multipole_threshold,
        reset_source_tree=false, reset_target_tree=false,
        upward_pass, horizontal_pass, downward_pass,
        nearfield, farfield, self_induced,
        unsort_source_bodies, unsort_target_bodies, nearfield_user
    )

    # visualize
    save_tree_source && (visualize(save_name_source, source_systems, source_tree))
    save_tree_target && (visualize(save_name_target, target_systems, target_tree))

    return source_tree, target_tree, m2l_list, direct_list, derivatives_switches
end

"""
    fmm!(systems; kwargs...)

Apply all interactions of `systems` acting on itself using the fast multipole method. Assumes compatibility functions have been overloaded for both source and target systems.

# Arguments

- `systems`: either

    - a system object for which compatibility functions have been overloaded, or
    - a tuple of system objects for which compatibility functions have been overloaded

# Optional Arguments

- `expansion_order::Int`: the expansion order to be used
- `leaf_size::Int`: maximum number of bodies from `systems` allowed in a leaf-level branch
- `multipole_threshold::Float64`: number between 0 and 1 (often denoted theta in [0,1]) controls the accuracy by determining the non-dimensional distance after which multipoles are used; 0 means an infinite distance (no error, high cost), and 1 means barely convergent (high error, low cost)
- `scalar_potential::Bool`: either a `::Bool` or a `::AbstractVector{Bool}` of length `length(systems)` indicating whether each system should receive a scalar potential from `source_systems`
- `velocity::Bool`: either a `::Bool` or a `::AbstractVector{Bool}` of length `length(systems)` indicating whether each system should receive a velocity from `source_systems`
- `velocity_gradient::Bool`: either a `::Bool` or a `::AbstractVector{Bool}` of length `length(systems)` indicating whether each system should receive a velocity gradient from `source_systems`
- `upward_pass::Bool`: whether or not to form the multipole expansions from source bodies and translate them upward in the source tree
- `horizontal_pass::Bool`: whether or not to transform multipole expansions from the source tree into local expansions in the target tree
- `downward_pass::Bool`: whether or not to translate local expansions down to the leaf level of the target tree and evaluate them
- `nearfield::Bool`: indicates whether near-field (comuted without multipoles) interactions should be included in the direct_list
- `farfield::Bool`: indicates whether far-field (comuted with multipoles) interactions should be included in the m2l_list
- `self_induced::Bool`: indicates whether to include the interactions of each leaf-level branch on itself in the direct_list
- `unsort_bodies::Bool`: indicates whether or not to undo the sort operation used to generate the octree for `systems`
- `shink_recenter::Bool`: indicates whether or not to resize branches for the octree after it is created to increase computational efficiency
- `save_tree::Bool`: indicates whether or not to save a VTK file for visualizing the octree
- `save_name::String`: name and path of the octree visualization if `save_tree == true`
- `nearfield_user::Bool`: indicates whether or not the `nearfield_user!` function should be used for nearfield interactions

"""
function fmm!(systems;
    error_method::ErrorMethod=UnequalBoxes(), predict_error::Val{PE}=Val(false),
    expansion_order::Union{Int,Dynamic{<:Any,<:Any}}=5, leaf_size=50, multipole_threshold=0.4,
    lamb_helmholtz::Bool=false,
    scalar_potential=true, velocity=true, velocity_gradient=true,
    upward_pass::Bool=true, horizontal_pass::Bool=true, downward_pass::Bool=true,
    nearfield::Bool=true, farfield::Bool=true, self_induced::Bool=true,
    unsort_bodies::Bool=true, shrink_recenter::Bool=false,
    save_tree::Bool=false, save_name="tree", nearfield_user::Bool=false
) where PE

    # create tree
    tree = Tree(systems; expansion_order=get_Pmax(expansion_order), leaf_size, shrink_recenter)

    # perform fmm
    m2l_list, direct_list, derivatives_switches = fmm!(tree, systems;
        error_method, predict_error,
        expansion_order,
        lamb_helmholtz,
        scalar_potential, velocity, velocity_gradient,
        multipole_threshold, reset_tree=false,
        upward_pass, horizontal_pass, downward_pass,
        nearfield, farfield, self_induced,
        unsort_bodies, nearfield_user
    )

    # visualize
    save_tree && (visualize(save_name, systems, tree))

    return tree, m2l_list, direct_list, derivatives_switches
end

"""
    fmm!(tree, systems; kwargs...)

Dispatches `fmm!` using an existing `::Tree`.

# Arguments

- `tree::Tree`: a `<:Tree` object (see [`Tree`](@ref))
- `systems`: either

    - a system object for which compatibility functions have been overloaded, or
    - a tuple of system objects for which compatibility functions have been overloaded

# Optional Arguments

- `multipole_threshold::Float64`: number between 0 and 1 (often denoted theta in [0,1]) controls the accuracy by determining the non-dimensional distance after which multipoles are used; 0 means an infinite distance (no error, high cost), and 1 means barely convergent (high error, low cost)
- `scalar_potential::Bool`: either a `::Bool` or a `::AbstractVector{Bool}` of length `length(systems)` indicating whether each system should receive a scalar potential from `source_systems`
- `velocity::Bool`: either a `::Bool` or a `::AbstractVector{Bool}` of length `length(systems)` indicating whether each system should receive a velocity from `source_systems`
- `velocity_gradient::Bool`: either a `::Bool` or a `::AbstractVector{Bool}` of length `length(systems)` indicating whether each system should receive a velocity gradient from `source_systems`
- `upward_pass::Bool`: whether or not to form the multipole expansions from source bodies and translate them upward in the source tree
- `horizontal_pass::Bool`: whether or not to transform multipole expansions from the source tree into local expansions in the target tree
- `downward_pass::Bool`: whether or not to translate local expansions down to the leaf level of the target tree and evaluate them
- `nearfield::Bool`: indicates whether near-field (computed without multipoles) interactions should be included in the direct_list
- `farfield::Bool`: indicates whether far-field (computed with multipoles) interactions should be included in the m2l_list
- `self_induced::Bool`: indicates whether to include the interactions of each leaf-level branch on itself in the direct_list
- `unsort_bodies::Bool`: indicates whether or not to undo the sort operation used to generate the octree for `systems`
- `nearfield_user::Bool`: indicates whether or not the `nearfield_user!` function should be used for nearfield interactions

"""
function fmm!(tree::Tree, systems;
    error_method::ErrorMethod=UnequalBoxes(), predict_error::Val{PE}=Val(false),
    expansion_order::Union{Int,Dynamic{<:Any,<:Any}}=5,
    multipole_threshold=0.4, reset_tree::Bool=true,
    lamb_helmholtz::Bool=false,
    scalar_potential=true, velocity=true, velocity_gradient=true,
    upward_pass::Bool=true, horizontal_pass::Bool=true, downward_pass::Bool=true,
    nearfield::Bool=true, farfield::Bool=true, self_induced::Bool=true,
    unsort_bodies::Bool=true, nearfield_user::Bool=false
) where PE

    # assemble derivatives switch
    derivatives_switches = DerivativesSwitch(scalar_potential, velocity, velocity_gradient, systems)

    # create interaction lists
    m2l_list, direct_list = build_interaction_lists(tree.branches, tree.branches, tree.leaf_index, multipole_threshold, farfield, nearfield, self_induced, error_method, expansion_order)

    # run fmm
    Pmax = get_Pmax(expansion_order)
    fmm!(tree, systems, m2l_list, direct_list, derivatives_switches, Val(Pmax);
        predict_error,
        lamb_helmholtz,
        reset_tree,
        nearfield, upward_pass, horizontal_pass, downward_pass,
        unsort_bodies, nearfield_user
    )

    return m2l_list, direct_list, derivatives_switches
end

"""
    fmm!(target_tree, target_systems, source_tree, source_systems; kwargs...)

Dispatches `fmm!` using existing `::Tree` objects.

# Arguments

- `target_tree::Tree`: a `<:Tree` object (see [`Tree`](@ref))
- `target_systems`: either

    - a system object for which compatibility functions have been overloaded, or
    - a tuple of system objects for which compatibility functions have been overloaded

- `source_tree::Tree`: a `<:Tree` object (see [`Tree`](@ref))
- `source_systems`: either

    - a system object for which compatibility functions have been overloaded, or
    - a tuple of system objects for which compatibility functions have been overloaded

# Optional Arguments

- `multipole_threshold::Float64`: number between 0 and 1 (often denoted theta in [0,1]) controls the accuracy by determining the non-dimensional distance after which multipoles are used; 0 means an infinite distance (no error, high cost), and 1 means barely convergent (high error, low cost)
- `scalar_potential::Bool`: either a `::Bool` or a `::AbstractVector{Bool}` of length `length(target_systems)` indicating whether each system should receive a scalar potential from `source_systems`
- `velocity::Bool`: either a `::Bool` or a `::AbstractVector{Bool}` of length `length(target_systems)` indicating whether each system should receive a velocity from `source_systems`
- `velocity_gradient::Bool`: either a `::Bool` or a `::AbstractVector{Bool}` of length `length(target_systems)` indicating whether each system should receive a velocity gradient from `source_systems`
- `reset_source_tree::Bool`: either a `::Bool` or a `::AbstractVector{Bool}` of length `length(source_systems)` indicating whether or not to reset the expansions of each source tree
- `reset_target_tree::Bool`: either a `::Bool` or a `::AbstractVector{Bool}` of length `length(source_systems)` indicating whether or not to reset the expansions of each target tree
- `upward_pass::Bool`: whether or not to form the multipole expansions from source bodies and translate them upward in the source tree
- `horizontal_pass::Bool`: whether or not to transform multipole expansions from the source tree into local expansions in the target tree
- `downward_pass::Bool`: whether or not to translate local expansions down to the leaf level of the target tree and evaluate them
- `nearfield::Bool`: indicates whether near-field (computed without multipoles) interactions should be included in the direct_list
- `farfield::Bool`: indicates whether far-field (computed with multipoles) interactions should be included in the m2l_list
- `self_induced::Bool`: indicates whether to include the interactions of each leaf-level branch on itself in the direct_list
- `unsort_source_bodies::Bool`: indicates whether or not to undo the sort operation used to generate the octree for `source_systems`
- `unsort_target_bodies::Bool`: indicates whether or not to undo the sort operation used to generate the octree for `target_systems`
- `nearfield_user::Bool`: indicates whether or not the `nearfield_user!` function should be used for nearfield interactions

"""
function fmm!(target_tree::Tree, target_systems, source_tree::Tree, source_systems;
    error_method::ErrorMethod=UnequalBoxes(), predict_error::Val{PE}=Val(false),
    expansion_order::Union{Int,Dynamic{<:Any,<:Any}}=5, multipole_threshold=0.4,
    scalar_potential=true, velocity=true, velocity_gradient=true,
    lamb_helmholtz::Bool=false,
    reset_source_tree::Bool=true, reset_target_tree::Bool=true,
    upward_pass::Bool=true, horizontal_pass::Bool=true, downward_pass::Bool=true,
    nearfield::Bool=true, farfield::Bool=true, self_induced::Bool=true,
    unsort_source_bodies::Bool=true, unsort_target_bodies::Bool=true,
    nearfield_user::Bool=false
) where PE

    # assemble derivatives switch
    derivatives_switches = DerivativesSwitch(scalar_potential, velocity, velocity_gradient, target_systems)

    # create interaction lists
    m2l_list, direct_list = build_interaction_lists(target_tree.branches, source_tree.branches, source_tree.leaf_index, multipole_threshold, farfield, nearfield, self_induced, error_method, expansion_order)

    # extract max expansion order
    Pmax = get_Pmax(expansion_order)

    # run fmm
    fmm!(target_tree, target_systems, source_tree, source_systems, m2l_list, direct_list, derivatives_switches, Val(Pmax);
        predict_error,
        lamb_helmholtz,
        reset_source_tree, reset_target_tree,
        nearfield, upward_pass, horizontal_pass, downward_pass,
        unsort_source_bodies, unsort_target_bodies, nearfield_user
    )

    return m2l_list, direct_list, derivatives_switches
end

"""
    fmm!(tree, systems; kwargs...)

Dispatches `fmm!` using an existing `::Tree`.

# Arguments

- `tree::Tree`: a `<:Tree` object (see [`Tree`](@ref))
- `systems`: either

    - a system object for which compatibility functions have been overloaded, or
    - a tuple of system objects for which compatibility functions have been overloaded

- `m2l_list::Vector{SVector{2,Int32}}`: list of branch index pairs `[i_target, i_source]` for which multipole expansions of the source branch are to be transformed to local expansions at the target branch
- `direct_list::Union{Vector{SVector{2,Int32}}, InteractionList}`: list of branch index pairs `[i_target, i_source]` for which interactions are to be evaluted without multipole expansion (i.e., directly); if `typeof(direct_list) <: InteractionList`, then prepared influence matrices are used rather than computing direct influences on the fly
- `derivatives_switches::Union{DerivativesSwitch, Tuple{<:DerivativesSwitch,...}}`: switch determining which of scalar potential, vector potential, velocity, and/or velocity gradient are to be computed for each target system

# Optional Arguments

- `upward_pass::Bool`: whether or not to form the multipole expansions from source bodies and translate them upward in the source tree
- `horizontal_pass::Bool`: whether or not to transform multipole expansions from the source tree into local expansions in the target tree
- `downward_pass::Bool`: whether or not to translate local expansions down to the leaf level of the target tree and evaluate them
- `unsort_bodies::Bool`: indicates whether or not to undo the sort operation used to generate the octree for `systems`

"""
function fmm!(tree::Tree, systems, m2l_list, direct_list, derivatives_switches, Pmax::Val;
    error_method::ErrorMethod=UnequalBoxes(), predict_error::Val{PE}=Val(false),
    lamb_helmholtz::Bool=false,
    reset_tree::Bool=true,
    nearfield::Bool=true, upward_pass::Bool=true, horizontal_pass::Bool=true, downward_pass::Bool=true,
    unsort_bodies::Bool=true, nearfield_user::Bool=false
) where PE

    fmm!(tree, systems, tree, systems, m2l_list, direct_list, derivatives_switches, Pmax;
        predict_error,
        lamb_helmholtz,
        reset_source_tree=reset_tree, reset_target_tree=false,
        nearfield, upward_pass, horizontal_pass, downward_pass,
        unsort_source_bodies=unsort_bodies, unsort_target_bodies=false,
        nearfield_user
    )

end

const WARNING_FLAG_EMPTY_SOURCE = Array{Bool,0}(undef)
WARNING_FLAG_EMPTY_SOURCE[] = false
const WARNING_FLAG_EMPTY_TARGET = Array{Bool,0}(undef)
WARNING_FLAG_EMPTY_TARGET[] = false

"""
    fmm!(target_tree, target_systems, source_tree, source_systems, m2l_list, direct_list, derivatives_switches; kwargs...)

Dispatches `fmm!` using existing `::Tree` objects.

# Arguments

- `target_tree::Tree`: a `<:Tree` object (see [`Tree`](@ref))
- `target_systems`: either

    - a system object for which compatibility functions have been overloaded, or
    - a tuple of system objects for which compatibility functions have been overloaded

- `source_tree::Tree`: a `<:Tree` object (see [`Tree`](@ref))
- `source_systems`: either

    - a system object for which compatibility functions have been overloaded, or
    - a tuple of system objects for which compatibility functions have been overloaded

- `m2l_list::Vector{SVector{2,Int32}}`: list of branch index pairs `[i_target, i_source]` for which multipole expansions of the source branch are to be transformed to local expansions at the target branch
- `direct_list::Union{Vector{SVector{2,Int32}}, InteractionList}`: list of branch index pairs `[i_target, i_source]` for which interactions are to be evaluted without multipole expansion (i.e., directly); if `typeof(direct_list) <: InteractionList`, then prepared influence matrices are used rather than computing direct influences on the fly
- `derivatives_switches::Union{DerivativesSwitch, Tuple{<:DerivativesSwitch,...}}`: switch determining which of scalar potential, vector potential, velocity, and/or velocity gradient are to be computed for each target system

# Optional Arguments

- `multipole_threshold::Float64`: number between 0 and 1 (often denoted theta in [0,1]) controls the accuracy by determining the non-dimensional distance after which multipoles are used; 0 means an infinite distance (no error, high cost), and 1 means barely convergent (high error, low cost)
- `scalar_potential::Bool`: either a `::Bool` or a `::AbstractVector{Bool}` of length `length(target_systems)` indicating whether each system should receive a scalar potential from `source_systems`
- `velocity::Bool`: either a `::Bool` or a `::AbstractVector{Bool}` of length `length(target_systems)` indicating whether each system should receive a velocity from `source_systems`
- `velocity_gradient::Bool`: either a `::Bool` or a `::AbstractVector{Bool}` of length `length(target_systems)` indicating whether each system should receive a velocity gradient from `source_systems`
- `upward_pass::Bool`: whether or not to form the multipole expansions from source bodies and translate them upward in the source tree
- `horizontal_pass::Bool`: whether or not to transform multipole expansions from the source tree into local expansions in the target tree
- `downward_pass::Bool`: whether or not to translate local expansions down to the leaf level of the target tree and evaluate them
- `nearfield::Bool`: indicates whether near-field (comuted without multipoles) interactions should be included
- `farfield::Bool`: indicates whether far-field (comuted with multipoles) interactions should be included
- `self_induced::Bool`: indicates whether to include the interactions of each leaf-level branch on itself
- `unsort_source_bodies::Bool`: indicates whether or not to undo the sort operation used to generate the octree for `source_systems`
- `unsort_target_bodies::Bool`: indicates whether or not to undo the sort operation used to generate the octree for `target_systems`

"""
function fmm!(target_tree::Tree, target_systems, source_tree::Tree, source_systems, m2l_list, direct_list, derivatives_switches, Pmax::Val;
    error_method::ErrorMethod=UnequalBoxes(), predict_error::Val{PE}=Val(false),
    lamb_helmholtz::Bool=false,
    reset_source_tree::Bool=true, reset_target_tree::Bool=true,
    nearfield::Bool=true, upward_pass::Bool=true, horizontal_pass::Bool=true, downward_pass::Bool=true,
    unsort_source_bodies::Bool=true, unsort_target_bodies::Bool=true,
    nearfield_user::Bool=false
) where PE

    # check if systems are empty
    n_sources = get_n_bodies(source_systems)
    n_targets = get_n_bodies(target_systems)

    # wrap lamb_helmholtz in Val
    lamb_helmholtz = Val(lamb_helmholtz)


    if n_sources > 0 && n_targets > 0

        # precompute y-axis rotation by π/2 matrices (if not already done)
        update_Hs_π2!(Hs_π2, Pmax)

        # precompute y-axis Wigner matrix normalization (if not already done)
        update_ζs_mag!(ζs_mag, Pmax)
        update_ηs_mag!(ηs_mag, Pmax)

        # reset multipole/local expansions
        reset_target_tree && (reset_expansions!(source_tree))
        reset_source_tree && (reset_expansions!(source_tree))

        # available threads
        n_threads = Threads.nthreads()


        # begin FMM
        if nearfield_user

            # allow nearfield_user! to be called concurrently with upward and horizontal passes
            t1 = Threads.@spawn nearfield && nearfield_user!(target_systems, target_tree.branches, derivatives_switches, source_systems, source_tree.branches, direct_list)
            n_threads_multipole = n_threads == 1 ? n_threads : n_threads - 1
            t2 = Threads.@spawn begin
                    upward_pass && upward_pass_multithread!(source_tree.branches, source_systems, Pmax, source_tree.levels_index, source_tree.leaf_index, n_threads_multipole)
                    horizontal_pass && length(m2l_list) > 0 && horizontal_pass_multithread!(target_tree.branches, source_tree.branches, m2l_list, lamb_helmholtz, Pmax, n_threads_multipole)
	                downward_pass && translate_locals_multithread!(target_tree.branches, Pmax, target_tree.levels_index, n_threads_multipole)
                end

            fetch(t1)
            fetch(t2)

            # local to body interaction
            downward_pass && local_to_body_multithread!(branches, systems, derivatives_switch, Pmax, lamb_helmholtz, leaf_index, n_threads)

        else # standard nearfield function

            nearfield && nearfield_multithread!(target_systems, target_tree.branches, derivatives_switches, source_systems, source_tree.branches, direct_list, n_threads)
            upward_pass && upward_pass_multithread!(source_tree.branches, source_systems, Pmax, lamb_helmholtz, source_tree.levels_index, source_tree.leaf_index, n_threads)
            horizontal_pass && length(m2l_list) > 0 && horizontal_pass_multithread!(target_tree.branches, source_tree.branches, m2l_list, lamb_helmholtz, Pmax, n_threads)
            downward_pass && downward_pass_multithread!(target_tree.branches, target_systems, derivatives_switches, Pmax, lamb_helmholtz, target_tree.levels_index, target_tree.leaf_index, n_threads)

        end

    else
        if n_sources == 0 && !WARNING_FLAG_EMPTY_SOURCE[]
            @warn "fmm! called but the source system is empty; foregoing calculation"
            WARNING_FLAG_EMPTY_SOURCE[] = true
        end
        if n_targets == 0 && !WARNING_FLAG_EMPTY_TARGET[]
            @warn "fmm! called but the target system is empty; foregoing calculation"
            WARNING_FLAG_EMPTY_TARGET[] = true
        end
    end

   # unsort bodies
    n_sources > 0 && unsort_source_bodies && unsort!(source_systems, source_tree)
    n_targets > 0 && unsort_target_bodies && unsort!(target_systems, target_tree)

end
