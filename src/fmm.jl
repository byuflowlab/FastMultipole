#####
##### direct interactions
#####
@inline function nearfield_singlethread!(target_systems::Tuple, direct_target_bodies::Tuple, derivatives_switches, source_system, direct_source_bodies, gpu=Val(false))
    for (target_system, target_index, derivatives_switch) in zip(target_systems, direct_target_bodies, derivatives_switches)
        nearfield_singlethread!(target_system, target_index, derivatives_switch, source_system, direct_source_bodies, gpu)
    end
end

@inline function nearfield_singlethread!(target_system, target_index::Vector, derivatives_switch, source_systems::Tuple, direct_source_bodies::Tuple, gpu=Val(false))
    for (source_system, source_index) in zip(source_systems, direct_source_bodies)
        nearfield_singlethread!(target_system, target_index, derivatives_switch, source_system, source_index, gpu)
    end
end

function nearfield_singlethread!(target_system, target_index::Vector, derivatives_switch, source_system, source_index::Vector, gpu=Val(false))
    _direct!(target_system, target_index, derivatives_switch, source_system, source_index, gpu)
end

function get_n_interactions(target_systems, direct_target_bodies::Tuple, source_system, source_index::Vector)
    n_interactions = 0
    for (target_system, target_index) in zip(target_systems, direct_target_bodies)
        n_interactions += get_n_interactions(target_system, target_index, source_system, source_index)
    end
    return n_interactions
end

function get_n_interactions(target_system, target_index::Vector, source_system, source_index::Vector)
    n_interactions = 0
    for (targets, sources) in zip(target_index, source_index)
        n_interactions += length(targets) * length(sources)
    end
    return n_interactions
end

function make_assignments!(assignments, direct_target_bodies::Tuple, source_index, n_threads, n_per_thread)
    i_start = 1
    i_thread = 1
    n_interactions = 0
    i_targets = 1
    for (i_end, sources) in enumerate(source_index)
        for i_target_bodies in eachindex(direct_target_bodies)
            n_interactions += length(direct_target_bodies[i_target_bodies][i_targets]) * length(sources)
        end
        if n_interactions >= n_per_thread
            assignments[i_thread] = i_start:i_end
            i_start = i_end + 1
            i_thread += 1
            n_interactions = 0
        end
        i_targets += 1
    end

    i_thread <= n_threads && (assignments[i_thread] = i_start:length(source_index))
end

function make_assignments!(assignments, target_index::Vector, source_index, n_threads, n_per_thread)
    i_start = 1
    i_thread = 1
    n_interactions = 0
    for (i_end,(targets, sources)) in enumerate(zip(target_index, source_index))
        n_interactions += length(targets) * length(sources)
        if n_interactions >= n_per_thread
            assignments[i_thread] = i_start:i_end
            i_start = i_end + 1
            i_thread += 1
            n_interactions = 0
        end
    end

    i_thread <= n_threads && (assignments[i_thread] = i_start:length(target_index))
end

@inline function execute_assignment!(target_system, target_index, derivatives_switch, source_system, source_index, assignment)
    _direct!(target_system, view(target_index, assignment), derivatives_switch, source_system, view(source_index, assignment), Val{false}())
end

@inline function execute_assignment!(target_systems::Tuple, direct_target_bodies::Tuple, derivatives_switches, source_system, source_index, assignment)
    for (target_system, target_index, derivatives_switch) in zip(target_systems, direct_target_bodies, derivatives_switches)
        execute_assignment!(target_system, target_index, derivatives_switch, source_system, source_index, assignment)
    end
end

function nearfield_multithread!(target_system, target_index, derivatives_switch, source_system, source_index::Vector, n_threads)
    ## load balance

    # total number of interactions
    n_interactions = get_n_interactions(target_system, target_index, source_system, source_index)

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
    make_assignments!(assignments, target_index, source_index, n_threads, n_per_thread)

    # execute tasks
    Threads.@threads for assignment in assignments
        execute_assignment!(target_system, target_index, derivatives_switch, source_system, source_index, assignment)
    end

    return nothing
end

function nearfield_multithread!(target_system, target_index, derivatives_switch, source_systems::Tuple, direct_source_bodies::Tuple, n_threads)
    for (source_system, source_index) in zip(source_systems, direct_source_bodies)
        nearfield_multithread!(target_system, target_index, derivatives_switch, source_system, source_index, n_threads)
    end
end

#------- UPWARD PASS -------#

function upward_pass_singlethread!(branches::AbstractVector{<:Branch{TF}}, systems, expansion_order::Val{P}, lamb_helmholtz) where {TF,P}

    # try preallocating one container to be reused
    Ts = zeros(length_Ts(P))
    eimϕs = zeros(2, P+1)
    weights_tmp_1 = initialize_expansion(P, TF)
    weights_tmp_2 = initialize_expansion(P, TF)

    # loop over branches
    for branch in view(branches,length(branches):-1:1) # no need to create a multipole expansion at the very top level
        if branch.n_branches == 0 # branch is a leaf
            body_to_multipole!(branch, systems, branch.harmonics, expansion_order)
        else # not a leaf
            # iterate over children
            for child_branch in view(branches, branch.branch_index)

                # try using redundant storage on each branch
                #weights_tmp_1 = child_branch.expansion_storage
                #weights_tmp_2 = branch.expansion_storage

                multipole_to_multipole!(branch, child_branch, weights_tmp_1, weights_tmp_2, Ts, eimϕs, ζs_mag, Hs_π2, expansion_order, lamb_helmholtz)
            end
        end
    end
end

function body_2_multipole_multithread!(branches, systems::Tuple, expansion_order::Val{P}, leaf_index, n_threads) where P
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
            for i_branch in view(leaf_index, leaf_assignment)
                branch = branches[i_branch]
                Threads.lock(branch.lock) do
                    body_to_multipole!(system, branch, branch.bodies_index[i_system], branch.harmonics, expansion_order)
                end
            end
        end
    end
end

function body_2_multipole_multithread!(branches, system, expansion_order::Val{P}, leaf_index, n_threads) where P
    ## load balance
    leaf_assignments = fill(1:0, n_threads)

    # total number of bodies
    n_bodies = 0
    for i_leaf in leaf_index
        n_bodies += length(branches[i_leaf].bodies_index)
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
        n_bodies += length(branches[i_leaf].bodies_index)
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
        for i_branch in view(leaf_index, leaf_assignments[i_thread])
            branch = branches[i_branch]
            body_to_multipole!(system, branch, branch.bodies_index, branch.harmonics, expansion_order)
        end
    end
end

function translate_multipoles_multithread!(branches, expansion_order::Val{P}, levels_index, n_threads) where P

    # iterate over levels
    for level_index in view(levels_index,length(levels_index):-1:2)

        # load balance
        n_branches = length(level_index)
        n_per_thread, rem = divrem(n_branches, n_threads)
        rem > 0 && (n_per_thread += 1)

        # if there are too many threads, we'll actually hurt performance
        n_per_thread < MIN_NPT_M2M && (n_per_thread = MIN_NPT_M2M)

        # assign thread start branches
        Threads.@threads for i_start in 1:n_per_thread:n_branches
            # get final branch
            i_end = min(n_branches, i_start+n_per_thread-1)

            # loop over branches
            for i_branch in view(level_index,i_start:i_end)
                child_branch = branches[i_branch]
                parent_branch = branches[child_branch.i_parent]
                Threads.lock(parent_branch.lock) do
                    M2M!(parent_branch, child_branch, child_branch.harmonics, child_branch.ML, expansion_order)
                end
            end
        end
    end
end

function upward_pass_multithread!(branches, systems, expansion_order, levels_index, leaf_index, n_threads)
    # create multipole expansions
    body_2_multipole_multithread!(branches, systems, expansion_order, leaf_index, n_threads)

    # m2m translation
    translate_multipoles_multithread!(branches, expansion_order, levels_index, n_threads)
end

#------- direct interaction matrix -------#

# TODO: add influence matrix approach to direct interactions

#------- horizontal pass -------#

function horizontal_pass_singlethread!(target_branches::Vector{<:Branch{TF}}, source_branches, m2l_list, expansion_order::Val{P}, lamb_helmholtz) where {TF,P}
    # preallocate containers to be reused
    weights_tmp_1 = initialize_expansion(P, TF)
    weights_tmp_2 = initialize_expansion(P, TF)
    Ts = zeros(length_Ts(P))
    eimϕs = zeros(2, P+1)

    for (i_target, j_source) in m2l_list
        target_branch = target_branches[i_target]
        source_branch = source_branches[j_source]
        #weights_tmp_1 = target_branch.expansion_storage
        #weights_tmp_2 = source_branch.expansion_storage
        multipole_to_local!(target_branch, source_branch, weights_tmp_1, weights_tmp_2, Ts, eimϕs, ζs_mag, ηs_mag, Hs_π2, expansion_order, lamb_helmholtz)
    end

end

# function horizontal_pass_singlethread!(target_branches, source_branches, m2l_list, expansion_order, harmonics, L)
#     for (i_target, j_source) in m2l_list
#         @lock target_branches[i_target].lock M2L!(target_branches[i_target], source_branches[j_source], harmonics, L, expansion_order)
#     end
# end

# @inline function preallocate_horizontal_pass(expansion_type, expansion_order)
#     harmonics = zeros(expansion_type, (expansion_order<<1 + 1)*(expansion_order<<1 + 1))
#     L = zeros(expansion_type, 4)
#     return harmonics, L
# end

# @inline function preallocate_horizontal_pass(expansion_type, expansion_order, n)
#     containers = [preallocate_horizontal_pass(expansion_type, expansion_order) for _ in 1:n]
# end

function horizontal_pass_multithread!(target_branches, source_branches::Vector{<:Branch{TF}}, m2l_list, expansion_order::Val{P}, n_threads) where {TF,P}
    # number of translations per thread
    n_per_thread, rem = divrem(length(m2l_list),n_threads)
    rem > 0 && (n_per_thread += 1)
    assignments = 1:n_per_thread:length(m2l_list)

    # preallocate memory
    # harmonics_preallocated = [initialize_harmonics(P,TF) for _ in 1:length(assignments)]
    # ML_preallocated = [initialize_ML(P,TF) for _ in 1:length(assignments)]

    # execute tasks
    Threads.@threads for i_thread in 1:length(assignments)
	# Threads.@threads for i_start in 1:n_per_thread:length(m2l_list)
        i_start = assignments[i_thread]
        i_stop = min(i_start+n_per_thread-1, length(m2l_list))
        # harmonics = harmonics_preallocated[i_thread]
        # ML = ML_preallocated[i_thread]
        # harmonics = initialize_harmonics(P,TF)
        # ML = initialize_ML(P,TF)
        for (i_target, j_source) in m2l_list[i_start:i_stop]
            Threads.lock(target_branches[i_target].lock) do
                M2L!(target_branches[i_target], source_branches[j_source], expansion_order)
                # M2L!(target_branches[i_target], source_branches[j_source], harmonics, ML, expansion_order)
            end
            # target_branch = target_branches[i_target]
            # Threads.@lock target_branch.lock M2L!(target_branch, source_branches[j_source], expansion_order)
        end
    end

    return nothing
end

#------- DOWNWARD PASS -------#

function preallocate_l2b(float_type, expansion_type, expansion_order::Val{P}, n_threads) where P
    containers = [preallocate_l2b(float_type, expansion_type, expansion_order) for _ in 1:n_threads]
    return containers
end

function downward_pass_singlethread!(branches::AbstractVector{<:Branch{TF}}, systems, expansion_order::Val{P}, lamb_helmholtz, derivatives_switches) where {TF,P}
    # try preallocating one container to be reused
    Ts = zeros(length_Ts(P))
    eimϕs = zeros(2, P+1)
    weights_tmp_1 = initialize_expansion(P, TF)
    weights_tmp_2 = initialize_expansion(P, TF)
    velocity_n_m = zeros(eltype(branches[1]), 2, 3, size(weights_tmp_1, 3))

    # loop over branches
    for branch in branches
        if branch.n_branches == 0 # leaf level
            evaluate_local!(systems, branch, branch.harmonics, velocity_n_m, expansion_order, lamb_helmholtz, derivatives_switches)
        else
            for i_child_branch in branch.branch_index
                local_to_local!(branch, branches[i_child_branch], weights_tmp_1, weights_tmp_2, Ts, eimϕs, ηs_mag, Hs_π2, expansion_order, lamb_helmholtz)
            end
        end
    end
end

function translate_locals_multithread!(branches, expansion_order::Val{P}, levels_index, n_threads) where P

    # iterate over levels
    for level_index in view(levels_index,2:length(levels_index))

        # divide chunks
        n_per_thread, rem = divrem(length(level_index),n_threads)
        rem > 0 && (n_per_thread += 1)

        # if there are too many threads, we'll actually hurt performance
        n_per_thread < MIN_NPT_L2L && (n_per_thread = MIN_NPT_L2L)

        # loop over branches
        Threads.@threads for i_start in 1:n_per_thread:length(level_index)
            i_stop = min(i_start+n_per_thread-1,length(level_index))

            # loop over branches
            for child_branch in view(branches,view(level_index,i_start:i_stop))
                L2L!(branches[child_branch.i_parent], child_branch, child_branch.harmonics, child_branch.ML, expansion_order)
            end
        end
    end
    return nothing
end

function local_2_body_multithread!(branches, systems, derivatives_switches, expansion_order, leaf_index, n_threads)
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
        for i_leaf in view(leaf_index,assignment)
            leaf = branches[i_leaf]
			L2B!(systems, leaf, derivatives_switches, expansion_order)
        end
    end
end

function downward_pass_multithread!(branches, systems, derivatives_switch, expansion_order, levels_index, leaf_index, n_threads)
    # m2m translation
	translate_locals_multithread!(branches, expansion_order, levels_index, n_threads)

    # local to body interaction
    local_2_body_multithread!(branches, systems, derivatives_switch, expansion_order, leaf_index, n_threads)
end

#####
##### create interaction lists
#####
function build_interaction_lists(target_branches, source_branches, source_leaf_index, multipole_threshold, farfield, nearfield, self_induced)

    # prepare containers
    m2l_list = Vector{SVector{2,Int32}}(undef,0)
    direct_list = Vector{SVector{2,Int32}}(undef,0)

    # populate lists
    build_interaction_lists!(m2l_list, direct_list, Int32(1), Int32(1), target_branches, source_branches, source_leaf_index, multipole_threshold, Val(farfield), Val(nearfield), Val(self_induced))

    # sort by target branch (helps with GPU computation)
    direct_target_bodies, direct_source_bodies = sort_list_by_target(direct_list, target_branches, source_branches, length(source_leaf_index))


    return m2l_list, direct_target_bodies, direct_source_bodies
end

@inline preallocate_bodies_index(T::Type{<:MultiBranch{<:Any,NT}}, n) where NT = Tuple(Vector{UnitRange{Int64}}(undef, n) for _ in 1:NT)
@inline preallocate_bodies_index(T::Type{<:SingleBranch}, n) = Vector{UnitRange{Int64}}(undef, n)

function build_interaction_lists!(m2l_list, direct_list, i_target, j_source, target_branches, source_branches, source_leaf_index, multipole_threshold, farfield::Val{ff}, nearfield::Val{nf}, self_induced::Val{si}) where {ff,nf,si}
    # unpack
    source_branch = source_branches[j_source]
    target_branch = target_branches[i_target]

    # determine multipole criterion
    spacing = source_branch.center - target_branch.center
    center_spacing_squared = spacing[1]*spacing[1] + spacing[2]*spacing[2] + spacing[3]*spacing[3]
    summed_radii_squared = target_branch.radius + source_branch.radius
    summed_radii_squared *= summed_radii_squared
    if center_spacing_squared * multipole_threshold * multipole_threshold >= summed_radii_squared # meet M2L criteria
        ff && push!(m2l_list, SVector{2}(i_target, j_source))
    elseif source_branch.n_branches == target_branch.n_branches == 0 # both leaves
        nf && (i_target!=j_source || si) && push!(direct_list, SVector{2}(i_target, j_source))
    elseif source_branch.n_branches == 0 || (target_branch.radius >= source_branch.radius && target_branch.n_branches != 0) # source is a leaf OR target is not a leaf and is bigger or the same size
        for i_child in target_branch.branch_index
            build_interaction_lists!(m2l_list, direct_list, i_child, j_source, target_branches, source_branches, source_leaf_index, multipole_threshold, farfield, nearfield, self_induced)
        end
    else # source is not a leaf AND target is a leaf or is smaller
        for j_child in source_branch.branch_index
            build_interaction_lists!(m2l_list, direct_list, i_target, j_child, target_branches, source_branches, source_leaf_index, multipole_threshold, farfield, nearfield, self_induced)
        end
    end
end

function sort_list_by_target(direct_list, target_branches::Vector{TT}, source_branches::Vector{TS}, n_leaves) where {TT,TS}
    target_counter = zeros(Int32, n_leaves)
    place_counter = zeros(Int32, n_leaves)
    direct_target_bodies = preallocate_bodies_index(TT, length(direct_list))
    direct_source_bodies = preallocate_bodies_index(TS, length(direct_list))

    # tally the contributions of each source
    for (i_target, j_source) in direct_list
        i_leaf = source_branches[i_target].i_leaf
        target_counter[i_leaf] += 1
    end

    # prepare place counter
    i_cum = 1
    for (i,n) in enumerate(target_counter)
        place_counter[i] = i_cum
        i_cum += n
    end

    # place interactions
    for (i,(i_target, j_source)) in enumerate(direct_list)
        i_leaf = target_branches[i_target].i_leaf
        update_direct_bodies!(direct_target_bodies, place_counter[i_leaf], target_branches[i_target].bodies_index)
        update_direct_bodies!(direct_source_bodies, place_counter[i_leaf], source_branches[j_source].bodies_index)
        place_counter[i_leaf] += 1
    end

    return direct_target_bodies, direct_source_bodies
end

function sort_list_by_source(direct_list, target_branches::Vector{TT}, source_branches::Vector{TS}, n_leaves) where {TT,TS}
    source_counter = zeros(Int32, n_leaves)
    place_counter = zeros(Int32, n_leaves)
    direct_target_bodies = preallocate_bodies_index(TT, length(direct_list))
    direct_source_bodies = preallocate_bodies_index(TS, length(direct_list))

    # tally the contributions of each source
    for (i_target, j_source) in direct_list
        j_leaf = source_branches[j_source].i_leaf
        source_counter[j_leaf] += 1
    end

    # prepare place counter
    i_cum = 1
    for (i,n) in enumerate(source_counter)
        place_counter[i] = i_cum
        i_cum += n
    end

    # place interactions
    for (i, (i_target, j_source)) in enumerate(direct_list)
        i_leaf = target_branches[i_target].i_leaf
        j_leaf = source_branches[j_source].i_leaf
        update_direct_bodies!(direct_target_bodies, place_counter[j_leaf], target_branches[i_target].bodies_index)
        update_direct_bodies!(direct_source_bodies, place_counter[j_leaf], source_branches[j_source].bodies_index)
        place_counter[j_leaf] += 1
    end

    return direct_target_bodies, direct_source_bodies
end


@inline function update_direct_bodies!(direct_bodies::Vector{<:UnitRange}, leaf_index, bodies_index::UnitRange)
    direct_bodies[leaf_index] = bodies_index
end

@inline function update_direct_bodies!(direct_bodies_list, leaf_index, bodies_indices::AbstractVector{<:UnitRange})
    for (direct_bodies, bodies_index) in zip(direct_bodies_list, bodies_indices)
        update_direct_bodies!(direct_bodies, leaf_index, bodies_index)
    end
end

function InteractionList(direct_list, target_systems, target_tree::Tree, source_systems, source_tree::Tree{TF,<:Any}, derivatives_switches) where TF
    # unpack tree
    leaf_index = source_tree.leaf_index

    # preallocate containers
    influence_matrices = Vector{Matrix{TF}}(undef, length(leaf_index))

    # determine strength dimensions
    strength_dims = get_strength_dims(source_systems)

    # add influence matrices
    for (i_matrix,i_source_branch) in enumerate(leaf_index)
        add_influence_matrix!(influence_matrices, i_matrix, target_systems, target_tree.branches, source_systems, source_tree.branches, i_source_branch, strength_dims, direct_list, derivatives_switches)
    end

    # create largest needed storage strength and influence vectors
    n_cols_max = 0
    n_rows_max = 0
    for influence_matrix in influence_matrices
        n_rows, n_cols = size(influence_matrix)
        n_cols_max = max(n_cols_max, n_cols)
        n_rows_max = max(n_rows_max, n_rows)
    end
    strengths = zeros(TF,n_cols_max)
    influence = zeros(TF,n_rows_max)

    return InteractionList{TF}(influence_matrices, strengths, influence, direct_list)
end

#####
##### running FMM
#####

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
    expansion_order=5, leaf_size_source=50, leaf_size_target=50, multipole_threshold=0.4,
    lamb_helmholtz::Bool=false,
    scalar_potential::Bool=true, velocity::Bool=true, velocity_gradient::Bool=true,
    upward_pass::Bool=true, horizontal_pass::Bool=true, downward_pass::Bool=true,
    nearfield::Bool=true, farfield::Bool=true, self_induced::Bool=true,
    unsort_source_bodies::Bool=true, unsort_target_bodies::Bool=true,
    source_shrink_recenter::Bool=true, target_shrink_recenter::Bool=true,
    save_tree_source::Bool=false, save_tree_target::Bool=false, save_name_source="source_tree", save_name_target="target_tree", gpu::Bool=false
)
    # check for duplicate systems
    target_systems = wrap_duplicates(target_systems, source_systems)

    # create trees
    source_tree = Tree(source_systems; expansion_order, leaf_size=leaf_size_source, shrink_recenter=source_shrink_recenter)
    target_tree = Tree(target_systems; expansion_order, leaf_size=leaf_size_target, shrink_recenter=target_shrink_recenter)

    # perform fmm
    m2l_list, direct_list, derivatives_switches = fmm!(target_tree, target_systems, source_tree, source_systems;
        scalar_potential, velocity, velocity_gradient,
        multipole_threshold,
        reset_source_tree=false, reset_target_tree=false,
        upward_pass, horizontal_pass, downward_pass,
        nearfield, farfield, self_induced,
        unsort_source_bodies, unsort_target_bodies, gpu
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
- `gpu::Bool`: indicates whether or not GPU is to be used for direct interactions

"""
function fmm!(systems;
    expansion_order=5, leaf_size=50, multipole_threshold=0.4,
    lamb_helmholtz::Bool=false,
    scalar_potential::Bool=true, velocity::Bool=true, velocity_gradient::Bool=true,
    upward_pass::Bool=true, horizontal_pass::Bool=true, downward_pass::Bool=true,
    nearfield::Bool=true, farfield::Bool=true, self_induced::Bool=true,
    unsort_bodies::Bool=true, shrink_recenter::Bool=true,
    save_tree::Bool=false, save_name="tree", gpu::Bool=false
)
    # create tree
    tree = Tree(systems; expansion_order, leaf_size, shrink_recenter)

    # perform fmm
    m2l_list, direct_target_bodies, direct_source_bodies, derivatives_switches = fmm!(tree, systems;
        scalar_potential, velocity, velocity_gradient,
        multipole_threshold, reset_tree=false,
        upward_pass, horizontal_pass, downward_pass,
        nearfield, farfield, self_induced,
        unsort_bodies, gpu
    )

    # visualize
    save_tree && (visualize(save_name, systems, tree))

    return tree, m2l_list, direct_target_bodies, direct_source_bodies, derivatives_switches
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
- `gpu::Bool`: indicates whether or not GPU is to be used for direct interactions

"""
function fmm!(tree::Tree, systems;
    multipole_threshold=0.4, reset_tree::Bool=true,
    lamb_helmholtz::Bool=false,
    scalar_potential::Bool=true, velocity::Bool=true, velocity_gradient::Bool=true,
    upward_pass::Bool=true, horizontal_pass::Bool=true, downward_pass::Bool=true,
    nearfield::Bool=true, farfield::Bool=true, self_induced::Bool=true,
    unsort_bodies::Bool=true, gpu::Bool=false
)

    # assemble derivatives switch
    derivatives_switches = DerivativesSwitch(scalar_potential, velocity, velocity_gradient, systems)

    # create interaction lists
    m2l_list, direct_target_bodies, direct_source_bodies = build_interaction_lists(tree.branches, tree.branches, tree.leaf_index, multipole_threshold, farfield, nearfield, self_induced)

    # run fmm
    fmm!(tree, systems, m2l_list, direct_target_bodies, direct_source_bodies,  derivatives_switches;
        lamb_helmholtz,
        reset_tree,
        nearfield, upward_pass, horizontal_pass, downward_pass,
        unsort_bodies, gpu
    )

    return m2l_list, direct_target_bodies, direct_source_bodies, derivatives_switches
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
- `gpu::Bool`: indicates whether or not GPU is to be used for direct interactions

"""
function fmm!(target_tree::Tree, target_systems, source_tree::Tree, source_systems;
    multipole_threshold=0.4,
    scalar_potential::Bool=true, velocity::Bool=true, velocity_gradient::Bool=true,
    lamb_helmholtz::Bool=false,
    reset_source_tree::Bool=true, reset_target_tree::Bool=true,
    upward_pass::Bool=true, horizontal_pass::Bool=true, downward_pass::Bool=true,
    nearfield::Bool=true, farfield::Bool=true, self_induced::Bool=true,
    unsort_source_bodies::Bool=true, unsort_target_bodies::Bool=true,
    gpu::Bool=false
)

    # assemble derivatives switch
    derivatives_switches = DerivativesSwitch(scalar_potential, velocity, velocity_gradient, target_systems)

    # create interaction lists
    m2l_list, direct_target_bodies, direct_source_bodies = build_interaction_lists(target_tree.branches, source_tree.branches, source_tree.leaf_index, multipole_threshold, farfield, nearfield, self_induced)

    # run fmm
    fmm!(target_tree, target_systems, source_tree, source_systems, m2l_list, direct_target_bodies, direct_source_bodies, derivatives_switches;
        lamb_helmholtz,
        reset_source_tree, reset_target_tree,
        nearfield, upward_pass, horizontal_pass, downward_pass,
        unsort_source_bodies, unsort_target_bodies, gpu
    )

    return m2l_list, direct_target_bodies, direct_source_bodies, derivatives_switches
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
function fmm!(tree::Tree, systems, m2l_list, direct_target_bodies, direct_source_bodies, derivatives_switches;
    lamb_helmholtz::Bool=false,
    reset_tree::Bool=true,
    nearfield::Bool=true, upward_pass::Bool=true, horizontal_pass::Bool=true, downward_pass::Bool=true,
    unsort_bodies::Bool=true, gpu::Bool=false
)

    fmm!(tree, systems, tree, systems, m2l_list, direct_target_bodies, direct_source_bodies, derivatives_switches;
        lamb_helmholtz,
        reset_source_tree=reset_tree, reset_target_tree=false,
        nearfield, upward_pass, horizontal_pass, downward_pass,
        unsort_source_bodies=unsort_bodies, unsort_target_bodies=false,
        gpu
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
function fmm!(target_tree::Tree, target_systems, source_tree::Tree, source_systems, m2l_list, direct_target_bodies, direct_source_bodies, derivatives_switches;
    lamb_helmholtz::Bool=false,
    reset_source_tree::Bool=true, reset_target_tree::Bool=true,
    nearfield::Bool=true, upward_pass::Bool=true, horizontal_pass::Bool=true, downward_pass::Bool=true,
    unsort_source_bodies::Bool=true, unsort_target_bodies::Bool=true,
    gpu::Bool=false
)
    # check if systems are empty
    n_sources = get_n_bodies(source_systems)
    n_targets = get_n_bodies(target_systems)

    if n_sources > 0 && n_targets > 0

        # precompute y-axis rotation by π/2 matrices (if not already done)
        update_Hs_π2!(Hs_π2, source_tree.expansion_order)

        # precompute y-axis Wigner matrix normalization (if not already done)
        update_ζs_mag!(ζs_mag, source_tree.expansion_order)
        update_ηs_mag!(ηs_mag, source_tree.expansion_order)

        # reset multipole/local expansions
        reset_target_tree && (reset_expansions!(source_tree))
        reset_source_tree && (reset_expansions!(source_tree))

        # available threads
        n_threads = Threads.nthreads()

        # gpu
        if gpu && n_threads > 1

            @sync begin

                Threads.@spawn nearfield_singlethread!(target_systems, direct_target_bodies, derivatives_switches, source_systems, direct_source_bodies, Val(gpu))
                upward_pass && upward_pass_multithread!(source_tree.branches, source_systems, source_tree.expansion_order, source_tree.levels_index, source_tree.leaf_index, n_threads-1)
                horizontal_pass && length(m2l_list) > 0 && horizontal_pass_multithread!(target_tree.branches, source_tree.branches, m2l_list, source_tree.expansion_order, n_threads-1)

            end
            downward_pass && downward_pass_multithread!(target_tree.branches, target_systems, derivatives_switches, target_tree.expansion_order, target_tree.levels_index, target_tree.leaf_index, n_threads)

        else

            # nearfield interactions
            if n_threads == 1 # && !gpu || gpu

                nearfield_singlethread!(target_systems, direct_target_bodies, derivatives_switches, source_systems, direct_source_bodies, Val(gpu))
                upward_pass && upward_pass_singlethread!(source_tree.branches, source_systems, source_tree.expansion_order, Val(lamb_helmholtz))
                horizontal_pass && length(m2l_list) > 0 && horizontal_pass_singlethread!(target_tree.branches, source_tree.branches, m2l_list, source_tree.expansion_order, Val(lamb_helmholtz))
                downward_pass && downward_pass_singlethread!(target_tree.branches, target_systems, target_tree.expansion_order, Val(lamb_helmholtz), derivatives_switches)

            else # n_threads > 1 && !gpu

                nearfield_multithread!(target_systems, direct_target_bodies, derivatives_switches, source_systems, direct_source_bodies, n_threads)
                upward_pass && upward_pass_multithread!(source_tree.branches, source_systems, source_tree.expansion_order, source_tree.levels_index, source_tree.leaf_index, n_threads)
                horizontal_pass && length(m2l_list) > 0 && horizontal_pass_multithread!(target_tree.branches, source_tree.branches, m2l_list, source_tree.expansion_order, n_threads)
                downward_pass && downward_pass_multithread!(target_tree.branches, target_systems, derivatives_switches, target_tree.expansion_order, target_tree.levels_index, target_tree.leaf_index, n_threads)

            end

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
