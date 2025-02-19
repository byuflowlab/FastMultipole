#------- direct interactions -------#

function nearfield_singlethread!(targets, target_branches, source_systems, source_branches, derivatives_switches, direct_list)
    # loop over sources
    t_nf = @MVector zeros(length(targets))
    for (i_source_system, source_system) in enumerate(source_systems)

        # perform direct interactions
        t_elapsed = @elapsed nearfield_loop!(targets, target_branches, source_system, i_source_system, source_branches, direct_list, derivatives_switches)
        t_nf[i_source_system] = t_elapsed
    end

    return t_nf
end

function nearfield_loop!(targets, target_branches, source_system, i_source_system, source_branches, direct_list, derivatives_switches)
    # loop over target systems
    for (i_target_system, target_system) in enumerate(targets)

        # extract derivatives switch
        derivatives_switch = derivatives_switches[i_target_system]

        # loop over direct list
        for (i_target, i_source) in direct_list

            # identify sources
            source_index = source_branches[i_source].bodies_index[i_source_system]

            # identify targets
            target_index = target_branches[i_target].bodies_index[i_target_system]

            # compute interaction
            _direct!(target_system, target_index, derivatives_switch, source_system, source_index)

        end

    end
end

function get_n_interactions(target_systems, target_branches, source_system, i_source_system, source_branches::Vector{<:Branch}, direct_list)
    n_interactions = 0
    for (i_target, i_source) in direct_list
        n_interactions += get_n_bodies(target_branches[i_target]) * length(source_branches[i_source].bodies_index[i_source_system])
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
        target_branch = target_branches[i_target]
        Threads.lock(target_branch.lock) do
            nearfield_singlethread!(target_systems, target_branch.bodies_index, derivatives_switches, source_system, source_branches[i_source].bodies_index[i_source_system])
        end
    end
end

function execute_assignment!(target_systems, target_branches, derivatives_switches, source_system, source_branches, direct_list, assignment)
    for i_interaction in assignment
        i_target, i_source = direct_list[i_interaction]
        target_branch = target_branches[i_target]
        Threads.lock(target_branch.lock) do
            nearfield_singlethread!(target_systems, target_branch.bodies_index, derivatives_switches, source_system, source_branches[i_source].bodies_index)
        end
    end
end

function nearfield_multithread!(systems, branches, derivatives_switch, direct_list, n_threads)
    if n_threads == 1 # single thread

        nearfield_singlethread!(systems, branches, derivatives_switch, direct_list)

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
    if DEBUG[]
        @show assignments
    end
    Threads.@threads for assignment in assignments
        execute_assignment!(target_systems, target_branches, derivatives_switches, source_system, source_branches, direct_list, assignment)
    end

end

"""
    nearfield_device!(target_systems, target_tree, derivatives_switches, source_systems, source_tree, direct_list)

User-defined function used to offload nearfield calculations to a device, such as GPU.

# Arguments

* `target_systems`: user-defined system on which `source_system` acts
* `target_tree::Tree`: octree object used to sort `target_systems`
* `derivatives_switches::Union{DerivativesSwitch, NTuple{N,DerivativesSwitch}}`: determines whether the scalar potential, velocity, and or velocity gradient should be calculated
* `source_systems`: user-defined system acting on `target_system`
* `source_tree::Tree`: octree object used to sort `target_systems`
* `direct_list::Vector{SVector{2,Int32}}`: each element `[i,j]` maps nearfield interaction from `source_tree.branches[j]` on `target_tree.branches[i]`

"""
function nearfield_device!(target_systems, target_tree::Tree, derivatives_switches, source_systems, source_tree::Tree, direct_list)
    @warn "nearfield_device! was called but hasn't been overloaded by the user"
end

"""
    nearfield_device!(target_systems, derivatives_switches, source_systems)

Dispatches `nearfield_device!` without having to build a `::Tree`. Performs all interactions.

# Arguments

* `target_systems`: user-defined system on which `source_system` acts
* `derivatives_switches::Union{DerivativesSwitch, NTuple{N,DerivativesSwitch}}`: determines whether the scalar potential, velocity, and or velocity gradient should be calculated
* `source_systems`: user-defined system acting on `target_system`

"""
function nearfield_device!(target_systems, derivatives_switches, source_systems)

    # get type
    x_source, _, _ = first_body_position(source_systems)
    TF_source = typeof(x_source)
    x_target, _, _ = first_body_position(target_systems)
    TF_target = typeof(x_target)
    TF = promote_type(TF_source, TF_target)

    # build target tree
    target_bodies_index = get_bodies_index(target_systems)
    n_branches, branch_index, i_parent, i_leaf_index = 0, 1, -1, 1
    center, source_radius, target_radius = SVector{3,TF}(0.0,0,0), zero(TF), zero(TF)
    source_box, target_box, expansion_order = SVector{6,TF}(0.0,0,0,0,0,0), SVector{3,TF}(0.0,0,0), 0
    target_branch = Branch(target_bodies_index, n_branches, branch_index, i_parent, i_leaf_index, center, source_radius, target_radius, source_box, target_box, expansion_order)
    levels_index, leaf_index, sort_index, inverse_sort_index, leaf_size = [1:1], [1], dummy_sort_index(target_systems), dummy_sort_index(target_systems), full_leaf_size(target_systems)
    target_tree = Tree([target_branch], levels_index, leaf_index, sort_index, inverse_sort_index, buffer, Val(expansion_order), leaf_size)

    if target_systems === source_systems
        source_tree = target_tree
    else
        # build source tree
        source_bodies_index = get_bodies_index(source_systems)
        source_branch = Branch(source_bodies_index, n_branches, branch_index, i_parent, i_leaf_index, center, source_radius, target_radius, source_box, target_box, expansion_order)
        sort_index, inverse_sort_index, leaf_size = dummy_sort_index(source_systems), dummy_sort_index(source_systems), full_leaf_size(source_systems)
        source_tree = Tree([source_branch], levels_index, leaf_index, sort_index, inverse_sort_index, buffer, expansion_order, leaf_size)
    end

    # build direct_list
    direct_list = [SVector{2,Int32}(1,1)]

    # call user-defined function
    nearfield_device!(target_systems, target_tree, derivatives_switches, source_systems, source_tree, direct_list)

end

#------- UPWARD PASS -------#

function upward_pass_singlethread_1!(tree::Tree{TF, <:Any, <:Any}, systems, expansion_order) where TF

    harmonics = initialize_harmonics(expansion_order, TF)

    # body_to_multipole
    for (i_system, system) in enumerate(systems)
        body_to_multipole!(tree, harmonics, system, i_system, expansion_order)
    end
end

function upward_pass_singlethread_2!(tree::Tree{TF,<:Any,<:Any}, expansion_order, lamb_helmholtz) where TF

    # try preallocating one container to be reused
    Ts = zeros(TF, length_Ts(expansion_order))
    eimϕs = zeros(TF, 2, expansion_order+1)
    weights_tmp_1 = initialize_expansion(expansion_order, TF)
    weights_tmp_2 = initialize_expansion(expansion_order, TF)

    # loop over branches
    for i_branch in length(tree.branches):-1:1 # no need to create a multipole expansion at the very top level
        branch = tree.branches[i_branch]
        expansion = view(tree.expansions, :, :, :, i_branch)

        if branch.n_branches !== 0 # branch is not a leaf
            # iterate over children
            for i_child in branch.branch_index
                child_branch = tree.branches[i_child]
                child_expansion = view(tree.expansions, :, :, :, i_child)
                multipole_to_multipole!(expansion, branch, child_expansion, child_branch, weights_tmp_1, weights_tmp_2, Ts, eimϕs, ζs_mag, Hs_π2, expansion_order, lamb_helmholtz)
            end
        end
    end
end

function upward_pass_singlethread!(tree::Tree, systems, expansion_order, lamb_helmholtz)
    upward_pass_singlethread_1!(tree, systems, expansion_order)
    upward_pass_singlethread_2!(tree, expansion_order, lamb_helmholtz)
end


# function upward_pass_singlethread!(branches::Vector{Branch{TF,N}}, systems, expansion_order, lamb_helmholtz, leaf_index) where {TF,N}
#
#     # try preallocating one container to be reused
#     Ts = zeros(TF, length_Ts(expansion_order))
#     eimϕs = zeros(TF, 2, expansion_order+1)
#     weights_tmp_1 = initialize_expansion(expansion_order, TF)
#     weights_tmp_2 = initialize_expansion(expansion_order, TF)
#
#     # multipole expansions
#     for i_branch in leaf_index
#         branch = branches[i_branch]
#         branch.source && body_to_multipole!(branch, systems, branch.harmonics, expansion_order)
#     end
#
#     # loop over branches
#     for i_branch in length(branches):-1:1 # no need to create a multipole expansion at the very top level
#         branch = branches[i_branch]
#
#         if branch.source && branch.n_branches !== 0 # branch is a source and is not a leaf
#             # iterate over children
#             for i_child in branch.branch_index
#                 child_branch = branches[i_child]
#                 multipole_to_multipole!(branch, child_branch, weights_tmp_1, weights_tmp_2, Ts, eimϕs, ζs_mag, Hs_π2, expansion_order, lamb_helmholtz)
#             end
#         end
#     end
# end

function body_to_multipole_multithread!(branches::Vector{<:Branch}, systems::Tuple, expansion_order, leaf_index, n_threads)
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

function translate_multipoles_multithread!(branches::Vector{Branch{TF,N}}, expansion_order, lamb_helmholtz, levels_index, n_threads) where {TF,N}

    # try preallocating one set of containers to be reused
    Ts = [zeros(TF, length_Ts(expansion_order)) for _ in 1:n_threads]
    eimϕs = [zeros(TF, 2, expansion_order+1) for _ in 1:n_threads]
    weights_tmp_1 = [initialize_expansion(expansion_order, TF) for _ in 1:n_threads]
    weights_tmp_2 = [initialize_expansion(expansion_order, TF) for _ in 1:n_threads]

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
                    multipole_to_multipole!(parent_expansion, parent_branch, child_expansion, child_branch, weights_tmp_1[i_task], weights_tmp_2[i_task], Ts[i_task], eimϕs[i_task], ζs_mag, Hs_π2, expansion_order, lamb_helmholtz)
                end
            end
        end
    end
end

function upward_pass_multithread!(branches, systems, expansion_order, lamb_helmholtz, levels_index, leaf_index, ::Val{n_threads}) where n_threads

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

# function horizontal_pass_debug!(target_system, target_branches::Vector{<:Branch{TF}}, source_system, source_branches, m2l_list, lamb_helmholtz, expansion_order, ε_abs; store="local", printstuff=74121) where TF
#
#     @assert isnothing(ε_abs) "ε_abs must be turned off for this function"
#
#     # increment the expansion order if ε_abs !== nothing
#     error_check = !(isnothing(ε_abs))
#
#     # preallocate containers to be reused
#     weights_tmp_1 = initialize_expansion(expansion_order + error_check, TF)
#     weights_tmp_2 = initialize_expansion(expansion_order + error_check, TF)
#     velocity_n_m = zeros(TF, 2, 3, size(weights_tmp_1, 3))
#     Ts = zeros(TF, length_Ts(expansion_order + error_check))
#     eimϕs = zeros(TF, 2, expansion_order+1+error_check)
#
#     # zero potential
#     for i in 1:get_n_bodies(target_system)
#         target_system[i,ScalarPotential()] = zero(TF)
#         target_system[i,Velocity()] = zero(SVector{3,TF})
#     end
#
#     # store observed and predicted errors
#     εs_obs = TF[]
#     εs_pred = TF[]
#
#     i_worst = 0
#     val_worst = 0.0
#
#     for (i_printstuff, (i_target, j_source)) in enumerate(m2l_list)
#         # get local expansion due to just this multipole expansion
#         target_branch = deepcopy(target_branches[i_target])
#         target_branch.local_expansion .= zero(TF)
#         source_branch = source_branches[j_source]
#         multipole_to_local!(target_branch, source_branch, weights_tmp_1, weights_tmp_2, Ts, eimϕs, ζs_mag, ηs_mag, Hs_π2, expansion_order, lamb_helmholtz, ε_abs)
#
#         # store previous potential/velocity
#         old_potential = [target_system[i, ScalarPotential()] for i in target_branch.bodies_index]
#         old_velocity = [target_system[i, Velocity()] for i in target_branch.bodies_index]
#
#         # calculate direct velocity
#         for i in target_branch.bodies_index
#             target_system[i, ScalarPotential()] = zero(TF)
#             target_system[i, Velocity()] = zero(SVector{3,TF})
#         end
#         derivatives_switch = DerivativesSwitch(true, true, false, target_system)
#         nearfield_singlethread!(target_system, target_branches, derivatives_switch, source_system, source_branches, [SVector{2,Int32}(i_target, j_source)])
#         direct_potential = [target_system[i, ScalarPotential()] for i in target_branch.bodies_index]
#         direct_velocity = [target_system[i, Velocity()] for i in target_branch.bodies_index]
#
#         # calculate multipole influence
#         for i in target_branch.bodies_index
#             target_system[i, ScalarPotential()] = zero(TF)
#             target_system[i, Velocity()] = zero(SVector{3,TF})
#         end
#         evaluate_multipole!(target_system, target_branch, source_branch, source_branch.harmonics, expansion_order, lamb_helmholtz, derivatives_switch)
#         multipole_potential = [target_system[i, ScalarPotential()] for i in target_branch.bodies_index]
#         multipole_velocity = [target_system[i, Velocity()] for i in target_branch.bodies_index]
#
#         # calculate expansion velocity
#         for i in target_branch.bodies_index
#             target_system[i, ScalarPotential()] = zero(TF)
#             target_system[i, Velocity()] = zero(SVector{3,TF})
#         end
#         evaluate_local!(target_system, target_branch, target_branch.harmonics, velocity_n_m, expansion_order, lamb_helmholtz, derivatives_switch)
#         expansion_potential = [target_system[i, ScalarPotential()] for i in target_branch.bodies_index]
#         expansion_velocity = [target_system[i, Velocity()] for i in target_branch.bodies_index]
#
#         # restore previous velocity
#         if store == "direct"
#             for (j,i) in enumerate(target_branch.bodies_index)
#                 target_system[i, ScalarPotential()] = old_potential[j] + direct_potential[j]
#                 target_system[i, Velocity()] = old_velocity[j] + direct_velocity[j]
#             end
#         elseif store == "multipole"
#             for (j,i) in enumerate(target_branch.bodies_index)
#                 target_system[i, ScalarPotential()] = old_potential[j] + multipole_potential[j]
#                 target_system[i, Velocity()] = old_velocity[j] + multipole_velocity[j]
#             end
#         elseif store == "local"
#             for (j,i) in enumerate(target_branch.bodies_index)
#                 target_system[i, ScalarPotential()] = old_potential[j] + expansion_potential[j]
#                 target_system[i, Velocity()] = old_velocity[j] + expansion_velocity[j]
#             end
#         end
#
#         # actual error
#         ε_obs = [norm(v_direct - v_expansion) for (v_direct, v_expansion) in zip(direct_velocity, expansion_velocity)]
#         ε_obs_mp = [norm(v_direct - v_expansion) for (v_direct, v_expansion) in zip(direct_velocity, multipole_velocity)]
#         ε_obs_max = maximum(abs.(ε_obs))
#
#         if ε_obs_max > val_worst
#             val_worst = ε_obs_max
#             i_worst = i_printstuff
#         end
#
#         if i_printstuff == printstuff
#             println("\nM2L DEBUG:")
#             @show mean(abs.(ε_obs_mp)) maximum(abs.(ε_obs_mp))
#             @show mean(abs.(ε_obs)) maximum(abs.(ε_obs))
#             @show target_branch.target_center source_branch.source_center target_branch.bodies_index source_branch.bodies_index target_branch.target_radius source_branch.source_radius target_branch.target_box source_branch.source_box
#             println()
#             println("for reproducibility:")
#             source_x = [source_system[i,Position()] for i in source_branch.bodies_index]
#             target_x = [target_system[i,Position()] for i in target_branch.bodies_index]
#             source_m = [source_system[i,Strength()] for i in source_branch.bodies_index]
#             target_m = [target_system[i,Strength()] for i in target_branch.bodies_index]
#             source_r = [source_system[i,Radius()] for i in source_branch.bodies_index]
#             target_r = [target_system[i,Radius()] for i in target_branch.bodies_index]
#
#             @show source_x source_m source_r target_x target_m target_r
#         end
#
#         # predicted error
#         ε_pred = error(target_branch, source_branch, expansion_order, RotatedCoefficients(), lamb_helmholtz)
#
#         # store results
#         push!(εs_obs, ε_obs_max)
#         push!(εs_pred, ε_pred)
#
#         if ε_obs_max > 0.1 && false
#             println("\nERROR > 0.1 ")
#             @show i_target j_source
#             @show ε_obs_max ε_pred
#             @show target_branch.target_radius
#             @show source_branch.source_radius
#             @show source_branch.bodies_index
#             dx = target_branch.target_center - source_branch.source_center
#             @show dx
#             @show (source_branch.source_radius + target_branch.target_radius) / norm(dx)
#
#             # ensure all bodies are contained inside
#             for i in target_branch.bodies_index
#                 dx = target_system[i,Position()] - target_branch.target_center
#                 @assert norm(dx) <= target_branch.target_radius "found a body outside: dx=$dx, r=$(target_branch.target_radius)"
#             end
#
#             for i in source_branch.bodies_index
#                 dx = source_system[i,Position()] - source_branch.source_center
#                 @assert norm(dx) <= source_branch.source_radius
#             end
#
#         end
#
#     end
#
#     println("\nWorst interaction:\n\ti=$i_worst\n")
#     return εs_obs, εs_pred
# end

function horizontal_pass_singlethread!(target_tree::Tree{TF1,<:Any,<:Any}, source_tree::Tree{TF2,<:Any,<:Any}, m2l_list, lamb_helmholtz, expansion_order, ε_abs; verbose=false) where {TF1,TF2}

    TF = promote_type(TF1, TF2)

    # increment the expansion order if ε_abs !== nothing
    error_check = !(isnothing(ε_abs))

    # preallocate containers to be reused
    weights_tmp_1 = initialize_expansion(expansion_order + error_check, TF)
    weights_tmp_2 = initialize_expansion(expansion_order + error_check, TF)
    Ts = zeros(TF, length_Ts(expansion_order + error_check))
    eimϕs = zeros(TF, 2, expansion_order + 1 + error_check)

    Pmax = 0
    error_success = true
    Ps = zeros(length(m2l_list))
    for (i,(i_target, j_source)) in enumerate(m2l_list)
        target_branch = target_tree.branches[i_target]
        target_expansion = view(target_tree.expansions, :, :, :, i_target)
        source_branch = source_tree.branches[j_source]
        source_expansion = view(source_tree.expansions, :, :, :, j_source)
        P, this_error_success = multipole_to_local!(target_expansion, target_branch, source_expansion, source_branch, weights_tmp_1, weights_tmp_2, Ts, eimϕs, ζs_mag, ηs_mag, Hs_π2, expansion_order, lamb_helmholtz, ε_abs)
        Pmax = max(P, Pmax)
        Ps[i] = P
        error_success = error_success && this_error_success
    end

    if verbose
        println("\n------- M2L Stats: -------")
        println("\n\tmean: ", mean(Ps))
        println("\tstd:  ", std(Ps))
        println("\tmax:  ", maximum(Ps))
        println("\tmin:  ", minimum(Ps))
        println("\n--------------------------\n")
    end

    return Pmax, error_success
end

function horizontal_pass_multithread!(target_tree::Tree{TF1,<:Any,<:Any}, source_tree::Tree{TF2,<:Any,<:Any}, m2l_list, lamb_helmholtz, expansion_order, ε_abs, n_threads) where {TF1, TF2}

    TF = promote_type(TF1, TF2)

    if n_threads == 1 # single thread

        horizontal_pass_singlethread!(target_branches, source_branches, m2l_list, lamb_helmholtz, expansion_order, ε_abs)

    else # multithread

        # number of translations per thread
        n_per_thread, rem = divrem(length(m2l_list),n_threads)
        rem > 0 && (n_per_thread += 1)
        assignments = 1:n_per_thread:length(m2l_list)

        # increment the expansion order if ε_abs !== nothing
        error_check = !(isnothing(ε_abs))

        # try preallocating one set of containers to be reused
        Ts = [zeros(TF, length_Ts(expansion_order + error_check)) for _ in 1:n_threads]
        eimϕs = [zeros(TF, 2, expansion_order + 1 + error_check) for _ in 1:n_threads]
        weights_tmp_1 = [initialize_expansion(expansion_order + error_check, TF) for _ in 1:n_threads]
        weights_tmp_2 = [initialize_expansion(expansion_order + error_check, TF) for _ in 1:n_threads]

        # execute tasks
        Threads.@threads for i_thread in 1:length(assignments)
            i_start = assignments[i_thread]
            i_stop = min(i_start+n_per_thread-1, length(m2l_list))
            for (i_target, j_source) in m2l_list[i_start:i_stop]
                Threads.lock(target_tree.branches[i_target].lock) do
                    multipole_to_local!(target_tree.branches[i_target], source_tree.branches[j_source], weights_tmp_1[i_thread], weights_tmp_2[i_thread], Ts[i_thread], eimϕs[i_thread], ζs_mag, ηs_mag, Hs_π2, expansion_order, lamb_helmholtz, ε_abs)
                end
            end
        end
    end

    return nothing
end

#------- DOWNWARD PASS -------#

function downward_pass_singlethread_1!(tree::Tree{TF,<:Any,<:Any}, expansion_order, lamb_helmholtz) where TF

    # try preallocating one container to be reused
    Ts = zeros(TF, length_Ts(expansion_order))
    eimϕs = zeros(TF, 2, expansion_order+1)
    weights_tmp_1 = initialize_expansion(expansion_order, TF)
    weights_tmp_2 = initialize_expansion(expansion_order, TF)

    # loop over branches
    for i_branch in 1:length(tree.branches)
        branch = tree.branches[i_branch]
        if branch.n_branches > 0 # if branch is a non-leaf target
            for i_child_branch in branch.branch_index
                child_branch = tree.branches[i_child_branch]
                child_expansion = view(tree.expansions, :, :, :, i_child_branch)
                branch_expansion = view(tree.expansions, :, :, :, i_branch)
                local_to_local!(child_expansion, child_branch, branch_expansion, branch, weights_tmp_1, weights_tmp_2, Ts, eimϕs, ηs_mag, Hs_π2, expansion_order, lamb_helmholtz)
            end
        end
    end
end

function downward_pass_singlethread_2!(tree::Tree{TF,<:Any,<:Any}, systems, expansion_order, lamb_helmholtz, derivatives_switches, velocity_n_m) where TF

    harmonics = initialize_harmonics(expansion_order)
    # loop over systems
    for (i_system, system) in enumerate(systems)
        evaluate_local!(system, i_system, tree, harmonics, velocity_n_m, expansion_order, lamb_helmholtz, derivatives_switches)
    end

end

@inline function downward_pass_singlethread!(tree, systems, expansion_order, lamb_helmholtz, derivatives_switches, velocity_n_m)

    downward_pass_singlethread_1!(tree, expansion_order, lamb_helmholtz)

    downward_pass_singlethread_2!(tree, systems, expansion_order, lamb_helmholtz, derivatives_switches, velocity_n_m)

end

# @inline function downward_pass_singlethread!(branches::Vector{<:Branch{TF,<:Any}}, systems, expansion_order, lamb_helmholtz, derivatives_switches, is_target) where TF
#     # try preallocating one container to be reused
#     Ts = zeros(TF, length_Ts(expansion_order))
#     eimϕs = zeros(TF, 2, expansion_order+1)
#     weights_tmp_1 = initialize_expansion(expansion_order, TF)
#     weights_tmp_2 = initialize_expansion(expansion_order, TF)
#     velocity_n_m = zeros(eltype(branches[1]), 2, 3, size(weights_tmp_1, 3))
#
#     # loop over branches
#     for branch in branches
#         if branch.target # if branch is a target
#             if branch.n_branches == 0 # leaf level
#                 evaluate_local!(systems, branch, branch.harmonics, velocity_n_m, expansion_order, lamb_helmholtz, derivatives_switches, is_target)
#             else
#                 for i_child_branch in branch.branch_index
#                     child_branch = branches[i_child_branch]
#                     child_branch.target && local_to_local!(child_branch, branch, weights_tmp_1, weights_tmp_2, Ts, eimϕs, ηs_mag, Hs_π2, expansion_order, lamb_helmholtz)
#                 end
#             end
#         end
#     end
# end

function translate_locals_multithread!(tree::Tree{TF,<:Any}, expansion_order, lamb_helmholtz, levels_index, n_threads) where TF

    # try preallocating one set of containers to be reused
    Ts = [zeros(TF, length_Ts(expansion_order)) for _ in 1:n_threads]
    eimϕs = [zeros(TF, 2, expansion_order+1) for _ in 1:n_threads]
    weights_tmp_1 = [initialize_expansion(expansion_order, TF) for _ in 1:n_threads]
    weights_tmp_2 = [initialize_expansion(expansion_order, TF) for _ in 1:n_threads]

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
                child_expansion = view(tree.expansions, :, :, :, i_child)
                branch_expansion = view(tree.expansions, :, :, :, child_branch.i_parent)
                local_to_local!(child_expansion, child_branch, branch_expansion, branches[child_branch.i_parent], weights_tmp_1[i_task], weights_tmp_2[i_task], Ts[i_task], eimϕs[i_task], ηs_mag, Hs_π2, expansion_order, lamb_helmholtz)
            end
        end
    end
    return nothing
end

function local_to_body_multithread!(branches::Vector{Branch{TF,<:Any}}, systems, derivatives_switches, expansion_order, lamb_helmholtz, leaf_index, n_threads) where TF
    # try preallocating one set of containers to be reused
    velocity_n_m = [initialize_velocity_n_m(expansion_order, TF) for _ in 1:n_threads]

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

# warn if lamb_helmholtz = true and scalar_potential=true

function warn_scalar_potential_with_lh(switch::DerivativesSwitch{PS,<:Any,<:Any}, lamb_helmholtz) where PS
    success = !(PS && lamb_helmholtz)
end

function warn_scalar_potential_with_lh(derivatives_switches::Tuple, lamb_helmholtz)
    success = true
    for switch in derivatives_switches
        success = success && warn_scalar_potential_with_lh(switch, lamb_helmholtz)
    end
    if !success
        @warn "\nScalar potential was requested with lamb_helmholtz=true; this may result in nonsensical potential predictions.\n"
    end
end

@inline function to_tuple(input::Tuple)
    return input
end

@inline function to_tuple(input)
    return (input,)
end

@inline function to_vector(input::AbstractVector, n)
    return SVector{n}(input)
end

@inline function to_vector(input, n)
    return SVector{n}(input for _ in 1:n)
end

# combine into tuple of unique systems
function unique_tuple(target_systems, source_systems, leaf_size_target::AbstractVector, leaf_size_source::AbstractVector)

    # unique tuple containing all systems
    systems = (target_systems..., (source_system for source_system in source_systems if !(source_system in target_systems))...)

    # tag telling which systems are targets
    is_target = SVector{length(systems),Bool}((true for _ in eachindex(target_systems))..., (false for _ in length(target_systems)+1:length(systems))...)

    # tag telling which systems are sources
    is_source = SVector{length(systems),Bool}(system in source_systems for system in systems)

    return systems, is_target, is_source
end

function match_unique(v_target::AbstractVector, v_source::AbstractVector, systems, target_systems, source_systems)
    # convert vector to match unique systems tuple
    v_unique = SVector{length(systems),eltype(v_target)}((val for val in v_target)..., (v_source[i] for i in eachindex(v_source) if !(source_systems[i] in target_systems))...)

    return v_unique
end

function system_to_target!(target, system)
    for i in 1:get_n_bodies(system)
        target[i,Position()] = system[i,Position()]
    end
end

function target_to_system!(system, target, ::DerivativesSwitch{PS,VS,GS}) where {PS,VS,GS}
    for i in 1:get_n_bodies(system)
        if PS
            system[i,ScalarPotential()] += target[i,ScalarPotential()]
        end
        if VS
            system[i,Velocity()] += target[i,Velocity()]
        end
        if GS
            system[i,VelocityGradient()] += target[i,VelocityGradient()]
        end
    end
end

fmm!(system; optargs...) = fmm!(system, system; optargs...)

function fmm!(target_systems, source_systems; optargs...)
    # promote arguments to Tuples
    target_systems = to_tuple(target_systems)
    source_systems = to_tuple(source_systems)

    return fmm!(target_systems, source_systems; optargs...)
end

function fmm!(target_systems::Tuple, source_systems::Tuple;
    leaf_size_target=nothing,
    leaf_size_source=default_leaf_size(source_systems),
    expansion_order=5,
    ε_abs=nothing,
    shrink_recenter=true,
    optargs...
)

    # promote leaf_size to vector
    leaf_size_source = to_vector(leaf_size_source, length(source_systems))
    leaf_size_target = to_vector(isnothing(leaf_size_target) ? minimum(leaf_size_source) : leaf_size_target, length(target_systems))

    # create trees
    target_tree = Tree(target_systems; expansion_order=expansion_order+!(isnothing(ε_abs)), leaf_size=leaf_size_target, shrink_recenter)
    source_tree = Tree(source_systems; expansion_order=expansion_order+!(isnothing(ε_abs)), leaf_size=leaf_size_source, shrink_recenter)

    # unsort systems
    unsort!(source_systems, source_tree)
    unsort!(target_systems, target_tree)

    return fmm!(target_systems, target_tree, source_systems, source_tree; expansion_order, leaf_size_source, ε_abs, optargs...)
end

function fmm!(target_systems::Tuple, target_tree::Tree, source_systems::Tuple, source_tree::Tree;
    leaf_size_source=default_leaf_size(source_systems), multipole_threshold=0.4,
    scalar_potential=true, velocity=true, velocity_gradient=true,
    farfield=true, nearfield=true, self_induced=true,
    optargs...
)

    # promote derivative arguments to a vector
    scalar_potential = to_vector(scalar_potential, length(target_systems))
    velocity = to_vector(velocity, length(target_systems))
    velocity_gradient = to_vector(velocity_gradient, length(target_systems))

    # assemble derivatives switch
    derivatives_switches = DerivativesSwitch(scalar_potential, velocity, velocity_gradient, target_systems)

    # create interaction lists
    m2l_list, direct_list = build_interaction_lists(target_tree.branches, source_tree.branches, leaf_size_source, multipole_threshold, farfield, nearfield, self_induced)

    # run fmm
    return fmm!(target_systems, target_tree, source_systems, source_tree, leaf_size_source, m2l_list, direct_list, derivatives_switches; optargs...)
end

"""
    fmm!(target_tree, target_systems, source_tree, source_systems, m2l_list, direct_list, derivatives_switches; kwargs...)

Dispatches `fmm!` using existing `::Tree` objects. Note that systems should be unsorted, as they will be resorted into their trees here.

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
function fmm!(target_systems::Tuple, target_tree::Tree, source_systems::Tuple, source_tree::Tree,  leaf_size_source, m2l_list, direct_list, derivatives_switches::Tuple;
    expansion_order=5, ε_abs=nothing, lamb_helmholtz::Bool=true,
    upward_pass::Bool=true, horizontal_pass::Bool=true, downward_pass::Bool=true,
    horizontal_pass_verbose::Bool=false,
    reset_target_tree::Bool=true, reset_source_tree::Bool=true,
    nearfield_device::Bool=false,
    tune=false,
)

    # check if systems are empty
    n_target_bodies = get_n_bodies(target_systems)
    n_source_bodies = get_n_bodies(source_systems)

    if n_target_bodies > 0 && n_source_bodies > 0

        # reset flags
        WARNING_FLAG_ERROR[] = true

        # check that lamb_helmholtz and ScalarPotential are not both true
        warn_scalar_potential_with_lh(derivatives_switches, lamb_helmholtz)

        # wrap lamb_helmholtz in Val
        lamb_helmholtz = Val(lamb_helmholtz)

        # increment the expansion order if ε_abs !== nothing
        error_check = !(isnothing(ε_abs))

        # precompute y-axis rotation by π/2 matrices (if not already done)
        update_Hs_π2!(Hs_π2, expansion_order + error_check)

        # precompute y-axis Wigner matrix normalization (if not already done)
        update_ζs_mag!(ζs_mag, expansion_order + error_check)
        update_ηs_mag!(ηs_mag, expansion_order + error_check)

        # available threads
        n_threads = Threads.nthreads()

        # reset trees
        reset_target_tree && (reset_expansions!(target_tree))
        reset_source_tree && (reset_expansions!(source_tree))

        # begin FMM
        if nearfield_device # use GPU

            # allow nearfield_device! to be called concurrently with upward and horizontal passes
            t1 = Threads.@spawn nearfield && nearfield_device!(target_systems, target_tree, derivatives_switches, source_systems, source_tree, direct_list)
            n_threads_multipole = n_threads == 1 ? n_threads : n_threads - 1
            t2 = Threads.@spawn begin
                    upward_pass && upward_pass_multithread!(source_tree.branches, source_systems, expansion_order + error_check, lamb_helmholtz, source_tree.levels_index, source_tree.leaf_index, n_threads_multipole)
                    horizontal_pass && length(m2l_list) > 0 && horizontal_pass_multithread!(target_tree.branches, source_tree.branches, m2l_list, lamb_helmholtz, expansion_order, ε_abs, n_threads_multipole)
	                downward_pass && translate_locals_multithread!(target_tree.branches, expansion_order, lamb_helmholtz, target_tree.levels_index, n_threads_multipole)
                end

            fetch(t1)
            fetch(t2)

            # local to body interaction
            downward_pass && local_to_body_multithread!(target_tree.branches, target_systems, derivatives_switches, Pmax, lamb_helmholtz, tree.leaf_index, n_threads)

        else # use CPU

            if n_threads == 1

                # sort into target tree
                resort!(target_systems, target_tree)

                # allocate target matrices
                targets = [zeros(i_VELOCITY_GRADIENT[end], get_n_bodies(system)) for system in target_systems]

                # copy system information to targets
                for (target, system) in zip(targets, target_systems)
                    system_to_target!(target, system)
                end

                # unsort targets NO NEED TO DO THIS YET?
                # unsort!(target_systems, target_tree)

                # sort into source tree
                resort!(source_systems, source_tree)

                # perform nearfield calculations
                t_direct = nearfield_singlethread!(targets, target_tree.branches, source_systems, source_tree.branches, derivatives_switches, direct_list)

                # check number of interactions
                if tune
                    n_interactions = 0
                    for i_source_system in eachindex(source_systems)
                        for (i_target, i_source) in direct_list
                            source_branch = source_tree.branches[i_source]
                            target_branch = target_tree.branches[i_target]
                            n_interactions += source_branch.n_bodies[i_source_system] * sum(target_branch.n_bodies)
                        end
                        t_direct[i_source_system] /= n_interactions
                    end
                end

                # farfield computations
                upward_pass && upward_pass_singlethread!(source_tree, source_systems, expansion_order + error_check, lamb_helmholtz)

                t_m2l = 0.0
                Pmax = 0
                error_success = true
                if horizontal_pass
                    t_m2l = @elapsed Pmax, error_success = horizontal_pass_singlethread!(target_tree, source_tree, m2l_list, lamb_helmholtz, expansion_order, ε_abs; verbose=horizontal_pass_verbose)
                end
                if !error_success
                    Pmax += 1
                end

                # @time downward_pass && downward_pass_singlethread!(tree.branches, tree.leaf_index, systems, expansion_order, lamb_helmholtz, derivatives_switches)
                if downward_pass
                    downward_pass_singlethread_1!(target_tree, expansion_order, lamb_helmholtz)
                    velocity_n_m = initialize_velocity_n_m(expansion_order, eltype(target_tree.branches[1]))
                    downward_pass_singlethread_2!(target_tree, targets, expansion_order, lamb_helmholtz, derivatives_switches, velocity_n_m)
                end

                # unsort sources (back in target tree permutation in case some systems are both sources and targets)
                unsort!(source_systems, source_tree)

                # copy results to target systems
                for (system, target, derivatives_switch) in zip(target_systems, targets, derivatives_switches)
                    target_to_system!(system, target, derivatives_switch)
                end

                # unsort targets
                unsort!(target_systems, target_tree)

                if tune
                    #=
                    println("\n#------- Tuning FastMultipole -------#\n\tDirect Cost:    $t_nf\n\tExpansion Cost: $t_exp")
                    suggestion = t_nf > t_exp * 1.25 ? "decreasing" : t_nf < t_exp * 0.75 ? "increasing" : "keeping current"
                    println("\tSuggest " * suggestion * " leaf size")
                    println("#------------------------------------#\n")
                    =#
                    if length(m2l_list) > 0
                        #--- compute optimal leaf_size_source ---#

                        # t per m2l transformation
                        t_m2l /= length(m2l_list)

                        # t_per_interaction * LS^2 = t_per_m2l
                        leaf_size_source = SVector{length(source_systems),Int}(Int(ceil(sqrt(t_m2l / t_direct[i]))) for i in eachindex(source_systems))
                    else
                        # make leaf size smaller so that some m2l operations exist
                        leaf_size_source = max.(leaf_size_source .>> 1, Ref(1))
                    end

                    expansion_order = Pmax
                end

            else

                nearfield && nearfield_multithread!(systems, tree.branches, derivatives_switches, direct_list, n_threads)
                upward_pass && upward_pass_multithread!(tree.branches, systems, expansion_order + error_check, lamb_helmholtz, tree.levels_index, tree.leaf_index, n_threads_2)
                horizontal_pass && length(m2l_list) > 0 && horizontal_pass_multithread!(tree.branches, tree.branches, m2l_list, lamb_helmholtz, expansion_order, ε_abs, n_threads)
                downward_pass && downward_pass_multithread!(tree.branches, systems, derivatives_switches, expansion_order, lamb_helmholtz, tree.levels_index, tree.leaf_index, n_threads)
            end

        end

    else

        @warn "fmm! called but either sources or targets are empty; foregoing calculation"

    end

    # unsort bodies
    # n_target_bodies > 0 && unsort_target_bodies && unsort!(target_systems, source_tree)
    # n_source_bodies > 0 && unsort_source_bodies && unsort!(source_systems, source_tree)

    # visualize tree
    # save_target_tree && (visualize(save_name, target_systems, target_tree))
    # save_source_tree && (visualize(save_name, source_systems, source_tree))

    return target_tree, source_tree, m2l_list, direct_list, derivatives_switches, leaf_size_source, expansion_order, error_success
end
