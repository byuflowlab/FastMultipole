#####
##### direct interactions
#####
function P2P!(target_system, target_branch::SingleBranch, derivatives_switch, source_system, source_branch::SingleBranch)
    _direct!(target_system, target_branch.bodies_index, derivatives_switch, source_system, source_branch.bodies_index)
end

function P2P!(target_system, target_branch::SingleBranch, derivatives_switch, source_systems, source_branch::MultiBranch)
    for (source_system, source_bodies_index) in zip(source_systems, source_branch.bodies_index)
        _direct!(target_system, target_branch.bodies_index, derivatives_switch, source_system, source_bodies_index)
    end
end

function P2P!(target_systems, target_branch::MultiBranch, derivatives_switches, source_system, source_branch::SingleBranch)
    for (target_system, target_bodies_index, derivatives_switch) in zip(target_systems, target_branch.bodies_index, derivatives_switches)
        _direct!(target_system, target_bodies_index, derivatives_switch, source_system, source_branch.bodies_index)
    end
end

function P2P!(target_systems, target_branch::MultiBranch, derivatives_switches, source_systems, source_branch::MultiBranch)
    for (source_system, source_bodies_index) in zip(source_systems, source_branch.bodies_index)
        for (target_system, target_bodies_index, derivatives_switch) in zip(target_systems, target_branch.bodies_index, derivatives_switches)
            _direct!(target_system, target_bodies_index, derivatives_switch, source_system, source_bodies_index)
        end
    end
end

function P2P!(target_systems, target_branch::MultiBranch, derivatives_switches, source_system, source_bodies_index::UnitRange)
    for (target_system, target_bodies_index, derivatives_switch) in zip(target_systems, target_branch.bodies_index, derivatives_switches)
        _direct!(target_system, target_bodies_index, derivatives_switch, source_system, source_bodies_index)
    end
end

function P2P!(target_system, target_branch::SingleBranch, derivatives_switch, source_system, source_bodies_index::UnitRange)
    _direct!(target_system, target_branch.bodies_index, derivatives_switch, source_system, source_bodies_index)
end

#####
##### upward pass
#####
function upward_pass_singlethread!(branches, systems, expansion_order::Val{P}) where P

    # loop over branches
    for branch in view(branches,length(branches):-1:1) # no need to create a multipole expansion at the very top level
        if branch.n_branches == 0 # branch is a leaf
            B2M!(branch, systems, branch.harmonics, expansion_order)
        else # not a leaf
            # iterate over children
            for child_branch in view(branches, branch.branch_index)
                M2M!(branch, child_branch, child_branch.harmonics, child_branch.ML, expansion_order)
            end
        end
    end
end

function body_2_multipole_multithread!(branches, systems::Tuple, expansion_order::Val{P}, leaf_index, concurrent_direct) where P
    ## load balance
    n_threads = Threads.nthreads() - 0^!concurrent_direct
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
                    B2M!(system, branch, branch.bodies_index[i_system], branch.harmonics, expansion_order)
                end
            end
        end
    end
end

function body_2_multipole_multithread!(branches, system, expansion_order::Val{P}, leaf_index, concurrent_direct) where P
    ## load balance
    n_threads = Threads.nthreads() - 0^!concurrent_direct
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
            B2M!(system, branch, branch.bodies_index, branch.harmonics, expansion_order)
        end
    end
end

function translate_multipoles_multithread!(branches, expansion_order::Val{P}, levels_index, concurrent_direct) where P
    # initialize memory
    n_threads = Threads.nthreads() - 0^!concurrent_direct

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

function upward_pass_multithread!(branches, systems, expansion_order, levels_index, leaf_index, concurrent_direct)
    # create multipole expansions
    body_2_multipole_multithread!(branches, systems, expansion_order, leaf_index, concurrent_direct)

    # m2m translation
    translate_multipoles_multithread!(branches, expansion_order, levels_index, concurrent_direct)
end

#####
##### horizontal pass
#####
function nearfield_singlethread!(target_system, target_tree::Tree, derivatives_switch, source_system, source_tree::Tree, direct_list)
    for (i_target, j_source) in direct_list
        P2P!(target_system, target_tree.branches[i_target], derivatives_switch, source_system, source_tree.branches[j_source])
    end
end

function update_strengths!(strengths, source_branch::MultiBranch, source_systems::Tuple)
    i_strength = 0
    for (source_system, bodies_index) in zip(source_systems, source_branch.bodies_index)
        i_strength = update_strengths!(strengths, i_strength, bodies_index, source_system)
    end
    return i_strength
end

function update_strengths!(strengths, source_branch, source_system)
    update_strengths!(strengths, 0, source_branch.bodies_index, source_system)
end

function update_strengths!(strengths, i_strength, bodies_index, source_system)
    strength_dims = get_strength_dims(source_system)
    for i_source_body in bodies_index
        σ = source_system[i_source_body, STRENGTH]
        for σi in σ
            i_strength += 1
            strengths[i_strength] = σi
        end
    end
    return i_strength
end

function update_influence!(target_systems::Tuple, this_influence::AbstractVector{TF}, i_row_prev, bodies_indices, derivatives_switches::Tuple) where TF
    for (target_system, bodies_index, switch) in zip(target_systems, bodies_indices, derivatives_switches)
        i_row_prev = update_influence!(target_system, this_influence, i_row_prev, bodies_index, switch)
    end
    return i_row_prev
end

function update_influence!(target_system, this_influence::AbstractVector{TF}, i_row_prev, bodies_index::UnitRange, derivatives_switch::DerivativesSwitch{PS,VPS,VS,GS}) where {TF,PS,VPS,VS,GS}

    for i_body in bodies_index
        if PS
            i_row_prev += 1
            ϕ = target_system[i_body,SCALAR_POTENTIAL]
            target_system[i_body,SCALAR_POTENTIAL] = ϕ + this_influence[i_row_prev]
        end

        if VPS
            i_row_prev += 3
            ψ = target_system[i_body,VECTOR_POTENTIAL]
            target_system[i_body,VECTOR_POTENTIAL] = ψ + SVector{3,TF}(this_influence[i_row_prev-2], this_influence[i_row_prev-1], this_influence[i_row_prev])
        end

        if VS
            i_row_prev += 3
            v = target_system[i_body,VELOCITY]
            target_system[i_body,VELOCITY] = v + SVector{3,TF}(this_influence[i_row_prev-2], this_influence[i_row_prev-1], this_influence[i_row_prev])
        end

        if GS
            i_row_prev += 9
            ∇v = target_system[i_body,VELOCITY_GRADIENT]
            target_system[i_body,VELOCITY_GRADIENT] = ∇v + SMatrix{3,3,TF,9}(this_influence[i_row_prev-8], this_influence[i_row_prev-7], this_influence[i_row_prev-6], this_influence[i_row_prev-5], this_influence[i_row_prev-4], this_influence[i_row_prev-3], this_influence[i_row_prev-2], this_influence[i_row_prev-1], this_influence[i_row_prev])
        end
    end

    return i_row_prev
end

function nearfield_singlethread!(target_system, target_tree::Tree, derivatives_switch, source_system, source_tree::Tree, interaction_list::InteractionList)
    @assert length(source_tree.leaf_index) == length(interaction_list.influence_matrices)

    # unpack
    strengths = interaction_list.strengths
    influence = interaction_list.influence

    # loop over source leaves
    for (i_source_branch, matrix) in zip(source_tree.leaf_index, interaction_list.influence_matrices)

        # update strengths
        source_branch = source_tree.branches[i_source_branch]
        i_strength = update_strengths!(strengths, source_branch, source_system)

        # obtain influence
        n_rows, n_cols = size(matrix)
        this_influence = view(influence, 1:n_rows)
        @assert i_strength == n_cols
        this_strength = view(strengths, 1:n_cols)
        mul!(this_influence, matrix, this_strength)

        # apply influence to targets
        i_row_prev = 0
        for (i_target, j_source) in interaction_list.direct_list
            if j_source == i_source_branch # found a target
                i_row_prev = update_influence!(target_system, this_influence, i_row_prev, target_tree.branches[i_target].bodies_index, derivatives_switch)
            end
        end
        @assert i_row_prev == n_rows

    end
end

function nearfield_multithread!(target_system, target_branches, derivatives_switch, source_systems::Tuple, source_branches, direct_list)
    ## load balance
    n_threads = Threads.nthreads()
    assignments = Vector{UnitRange{Int64}}(undef,n_threads)

    for (i_source_system, source_system) in enumerate(source_systems)
        # total number of interactions
        n_interactions = 0
        for (i_target, i_source) in direct_list
            target_leaf = view(target_branches,i_target)
            source_leaf = view(source_branches,i_source)
            n_interactions += get_n_bodies(target_leaf[].bodies_index) * get_n_bodies(source_leaf[].bodies_index[i_source_system])
        end

        # interactions per thread
        n_per_thread, rem = divrem(n_interactions, n_threads)
        rem > 0 && (n_per_thread += 1)

        # if there are too many threads, we'll actually hurt performance
        n_per_thread < MIN_NPT_NF && (n_per_thread = MIN_NPT_NF)

        # create assignments
        for i in eachindex(assignments)
            assignments[i] = 1:0
        end
        i_start = 1
        i_thread = 1
        n_interactions = 0
        for (i_end,(i_target, j_source)) in enumerate(direct_list)
            target_leaf = view(target_branches,i_target)
            source_leaf = view(source_branches,j_source)
            n_interactions += get_n_bodies(target_leaf[].bodies_index) * get_n_bodies(source_leaf[].bodies_index[i_source_system])
            if n_interactions >= n_per_thread
                assignments[i_thread] = i_start:i_end
                i_start = i_end+1
                i_thread += 1
                n_interactions = 0
            end
        end
        i_thread <= n_threads && (assignments[i_thread] = i_start:length(direct_list))

        # execute tasks
        Threads.@threads for i_thread in eachindex(assignments)
            assignment = assignments[i_thread]
            for (i_target, j_source) in view(direct_list, assignment)
                target_branch = target_branches[i_target]
                source_bodies_index = source_branches[j_source].bodies_index[i_source_system]
                Threads.lock(target_branch.lock) do
                    P2P!(target_system, target_branch, derivatives_switch, source_system, source_bodies_index)
                end
            end
        end
    end

    return nothing
end

function nearfield_multithread!(target_system, target_branches, derivatives_switch, source_system, source_branches, direct_list)
    ## load balance
    n_threads = Threads.nthreads()
    assignments = Vector{UnitRange{Int64}}(undef,n_threads)

    # total number of interactions
    n_interactions = 0
    for (i_target, i_source) in direct_list
        target_leaf = view(target_branches,i_target)
        source_leaf = view(source_branches,i_source)
        n_interactions += get_n_bodies(target_leaf[].bodies_index) * get_n_bodies(source_leaf[].bodies_index)
    end

    # interactions per thread
    n_per_thread, rem = divrem(n_interactions, n_threads)
    rem > 0 && (n_per_thread += 1)

    # if there are too many threads, we'll actually hurt performance
    n_per_thread < MIN_NPT_NF && (n_per_thread = MIN_NPT_NF)

    # create assignments
    for i in eachindex(assignments)
        assignments[i] = 1:0
    end
    i_start = 1
    i_thread = 1
    n_interactions = 0
    for (i_end,(i_target, i_source)) in enumerate(direct_list)
        target_leaf = view(target_branches,i_target)
        source_leaf = view(source_branches,i_source)
        n_interactions += get_n_bodies(target_leaf[].bodies_index) * get_n_bodies(source_leaf[].bodies_index)
        if n_interactions >= n_per_thread
            assignments[i_thread] = i_start:i_end
            i_start = i_end+1
            i_thread += 1
            n_interactions = 0
        end
    end
    i_thread <= n_threads && (assignments[i_thread] = i_start:length(direct_list))

    # execute tasks
    Threads.@threads for assignment in assignments
        for (i_target, j_source) in view(direct_list, assignment)
            target_branch = target_branches[i_target]
            Threads.lock(target_branch.lock) do
                P2P!(target_system, target_branch, derivatives_switch, source_system, source_branches[j_source])
            end
        end
    end

    return nothing
end

function horizontal_pass_singlethread!(target_branches, source_branches, m2l_list, expansion_order)
    for (i_target, j_source) in m2l_list
        M2L!(target_branches[i_target], source_branches[j_source], expansion_order)
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

function horizontal_pass_multithread!(target_branches, source_branches::Vector{<:Branch{TF}}, m2l_list, expansion_order::Val{P}, concurrent_direct) where {TF,P}
    # number of translations per thread
    n_threads = Threads.nthreads() - 0^!concurrent_direct
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

#####
##### downward pass
#####
function preallocate_l2b(float_type, expansion_type, expansion_order::Val{P}, n_threads) where P
    containers = [preallocate_l2b(float_type, expansion_type, expansion_order) for _ in 1:n_threads]
    return containers
end

function downward_pass_singlethread!(branches, systems, derivatives_switches, expansion_order::Val{P}) where P
    regular_harmonics = zeros(eltype(branches[1].multipole_expansion), 2, (P+1)*(P+1))
    for branch in branches
        if branch.n_branches == 0 # leaf level
			L2B!(systems, branch, derivatives_switches, expansion_order)
        else
            for child_branch in view(branches,branch.branch_index)
                L2L!(branch, child_branch, regular_harmonics, branch.ML, expansion_order)
            end
        end
    end
end

function translate_locals_multithread!(branches, expansion_order::Val{P}, levels_index) where P
    # initialize memory
    n_threads = Threads.nthreads()

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

function local_2_body_multithread!(branches, systems, derivatives_switches, expansion_order, leaf_index)
    # create assignments
    n_threads = Threads.nthreads()

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

function downward_pass_multithread!(branches, systems, derivatives_switch, expansion_order, levels_index, leaf_index)
    # m2m translation
	translate_locals_multithread!(branches, expansion_order, levels_index)

    # local to body interaction
    local_2_body_multithread!(branches, systems, derivatives_switch, expansion_order, leaf_index)
end

#####
##### create interaction lists
#####
function build_interaction_lists(target_branches, source_branches, source_leaf_index, multipole_threshold, farfield, nearfield, self_induced)
    m2l_list = Vector{SVector{2,Int32}}(undef,0)
    direct_list = Vector{SVector{2,Int32}}(undef,0)
    build_interaction_lists!(m2l_list, direct_list, Int32(1), Int32(1), target_branches, source_branches, source_leaf_index, multipole_threshold, Val(farfield), Val(nearfield), Val(self_induced))
    return m2l_list, direct_list
end

function build_interaction_lists!(m2l_list, direct_list, i_target, j_source, target_branches, source_branches, source_leaf_index, multipole_threshold, farfield::Val{ff}, nearfield::Val{nf}, self_induced::Val{si}, sort_direct=false) where {ff,nf,si}
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

    sort_direct && sort_list!(direct_list, length(source_tree.leaf_index))
end

function sort_list!(direct_list, n_leaves)
    source_counter = zeros(Int32, n_leaves)
    place_counter = zeros(Int32, n_leaves)
    buffer = similar(direct_list)

    # tally the contributions of each source
    for (i_target, j_source) in direct_list
        source_counter[j_source] += 1
    end

    # place interactions
    for interaction in direct_list
        j_source = interaction[2]
    end
end

@inline function get_strength_dims(systems::Tuple)
    return SVector{length(systems),Int}(get_strength_dims(system) for system in systems)
end

@inline function get_strength_dims(system)
    return length(system[1,STRENGTH])
end

@inline function get_one(TF,::Val{n},i) where {n}
    return SVector{n,TF}(0.0^(i!=j) for j in 1:n)
end

@inline function get_one(TF,::Val{1},i)
    return one(TF)
end

@inline function get_n_rows(bodies_index::UnitRange, ::DerivativesSwitch{PS,VPS,VS,GS}) where {PS,VPS,VS,GS}
    return length(bodies_index) * (0^!PS + 3*0^!VPS + 3*0^!VS + 9*0^!GS)
end

@inline function get_n_rows(bodies_indices, derivatives_switches::Tuple)
    return sum(get_n_rows(bodies_index, switch) for (bodies_index, switch) in zip(bodies_indices, derivatives_switches))
end

"multiple source systems"
function populate_influence_matrix!(matrix, target_systems, direct_list, target_branches, derivatives_switches, source_systems, i_source_branch, source_bodies_indices, strength_dims)
    i_col_prev = 0
    for (source_system, source_bodies_index, strength_dim) in zip(source_systems, source_bodies_indices, strength_dims)
        i_col_prev = populate_influence_matrix!(matrix, i_col_prev, target_systems, direct_list, target_branches, derivatives_switches, source_system, i_source_branch, source_bodies_index, strength_dim)
    end
end

"single source system (without specifying i_col_prev)"
function populate_influence_matrix!(matrix::Matrix{TF}, target_systems, direct_list, target_branches, derivatives_switches, source_system, i_source_branch, source_bodies_index::UnitRange, strength_dims::Int) where TF
    populate_influence_matrix!(matrix, 0, target_systems, direct_list, target_branches, derivatives_switches, source_system, i_source_branch, source_bodies_index, strength_dims)
end

"single source system"
function populate_influence_matrix!(matrix::Matrix{TF}, i_col_prev, target_systems, direct_list, target_branches, derivatives_switches, source_system, i_source_branch, source_bodies_index::UnitRange, strength_dims::Int) where TF
    strength_dims_val = Val(strength_dims)
    for i_source_body in source_bodies_index
        strength = source_system[i_source_body,STRENGTH]
        for i_strength_component in 1:strength_dims
            i_col_prev += 1
            source_system[i_source_body,STRENGTH] = get_one(TF,strength_dims_val,i_strength_component)
            i_row_prev = 0

            # search for targets
            for (i_target, j_source) in direct_list
                if j_source == i_source_branch # found a target
                    target_bodies_index = target_branches[i_target].bodies_index
                    i_row_prev = populate_influence_matrix!(matrix, i_row_prev, i_col_prev, target_systems, target_bodies_index, derivatives_switches, source_system, i_source_body)
                end
            end

        end
        source_system[i_source_body,STRENGTH] = strength
    end
    return i_col_prev
end

"single source body, multiple target systems"
function populate_influence_matrix!(matrix, i_row_prev, i_col, target_systems::Tuple, target_bodies_indices, derivatives_switches::Tuple, source_system, i_source_body)
    for (target_system, target_bodies_index, switch) in zip(target_systems, target_bodies_indices, derivatives_switches)
        i_row_prev = populate_influence_matrix!(matrix, i_row_prev, i_col, target_system, target_bodies_index, switch, source_system, i_source_body)
    end
    return i_row_prev
end

"single source body, single target system"
function populate_influence_matrix!(matrix::Matrix{TF}, i_row_prev, i_col, target_system, target_bodies_index::UnitRange, derivatives_switch::DerivativesSwitch{PS,VPS,VS,GS}, source_system, i_source_body) where {TF,PS,VPS,VS,GS}

    # loop over targets
    for i_target_body in target_bodies_index

        # save old influence and reset
        if PS
            ϕ_old = target_system[i_target_body, SCALAR_POTENTIAL]
            target_system[i_target_body, SCALAR_POTENTIAL] = zero(TF)
        end
        if VPS
            ψ_old = SVector{3,TF}(target_system[i_target_body, VECTOR_POTENTIAL])
            target_system[i_target_body, VECTOR_POTENTIAL] = zero(SVector{3,TF})
        end
        if VS
            v_old = SVector{3,TF}(target_system[i_target_body, VELOCITY])
            target_system[i_target_body, VELOCITY] = zero(SVector{3,TF})
        end
        if GS
            ∇v_old = SMatrix{3,3,TF,9}(target_system[i_target_body, VELOCITY_GRADIENT])
            target_system[i_target_body, VELOCITY_GRADIENT] = zero(SMatrix{3,3,TF,9})
        end

        # compute unit influence of source
        direct!(target_system, i_target_body, derivatives_switch, source_system, i_source_body)

        # update influence matrix
        if PS
            i_row_prev += 1
            matrix[i_row_prev, i_col] = target_system[i_target_body, SCALAR_POTENTIAL]
        end

        if VPS
            ψ = target_system[i_target_body, VECTOR_POTENTIAL]
            for i in eachindex(ψ)
                i_row_prev += 1
                matrix[i_row_prev, i_col] = ψ[i]
            end
        end

        if VS
            v = target_system[i_target_body, VELOCITY]
            for i in eachindex(v)
                i_row_prev += 1
                matrix[i_row_prev, i_col] = v[i]
            end
        end

        if GS
            ∇v = target_system[i_target_body, VELOCITY_GRADIENT]
            for i in eachindex(∇v)
                i_row_prev += 1
                matrix[i_row_prev, i_col] = ∇v[i]
            end
        end

        # restore old influence
        if PS
            target_system[i_target_body, SCALAR_POTENTIAL] = ϕ_old
        end
        if VPS
            target_system[i_target_body, VECTOR_POTENTIAL] = ψ_old
        end
        if VS
            target_system[i_target_body, VELOCITY] = v_old
        end
        if GS
            target_system[i_target_body, VELOCITY_GRADIENT] = ∇v_old
        end

    end

    return i_row_prev
end

@inline function get_n_strengths(bodies_index::UnitRange, strength_dim::Int)
    return length(bodies_index) * strength_dim
end

@inline function get_n_strengths(bodies_indices, strength_dims)
    n = 0
    for (bodies_index, strength_dim) in zip(bodies_indices, strength_dims)
        n += get_n_strengths(bodies_index, strength_dim)
    end
    return n
end

function add_influence_matrix!(influence_matrices::Vector{Matrix{TF}}, i_matrix, target_systems, target_branches, source_systems, source_branches, i_source_branch, strength_dims, direct_list, derivatives_switches) where TF

    # unpack
    source_branch = source_branches[i_source_branch]

    # number of source strength values
    n_cols = get_n_strengths(source_branch.bodies_index, strength_dims)

    # number of target locations
    n_rows = 0
    for (i_target, j_source) in direct_list
        j_source == i_source_branch && (n_rows += get_n_rows(target_branches[i_target].bodies_index, derivatives_switches))
    end

    # preallocate influence matrix
    matrix = Matrix{TF}(undef, n_rows, n_cols)

    # populate influence matrix
    populate_influence_matrix!(matrix, target_systems, direct_list, target_branches, derivatives_switches, source_systems, i_source_branch, source_branch.bodies_index, strength_dims)

    # update matrices
    influence_matrices[i_matrix] = matrix

end

@inline function get_influence_storage(::DerivativesSwitch{PS,VPS,VS,GS}, tree::Tree{TF,<:Any}) where {PS,VPS,VS,GS,TF}
    n_rows = 1^PS + 3^VPS + 3^VS + 9^GS
    n_cols = tree.leaf_size
    return Matrix{TF}(undef, n_rows, n_cols)
end

@inline function get_influence_storage(derivatives_switches::Tuple, tree::Tree)
    return Tuple(get_influence_storage(switch,tree) for switch in derivatives_switches)
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

    - a system object for which compatibility functions have been overloaded, or
    - a tuple of system objects for which compatibility functions have been overloaded

- `source_systems`: either

    - a system object for which compatibility functions have been overloaded, or
    - a tuple of system objects for which compatibility functions have been overloaded

# Optional Arguments

- `expansion_order::Int`: the expansion order to be used
- `leaf_size_source::Int`: maximum number of bodies from `source_systems` allowed in a leaf-level branch
- `leaf_size_target::Int`: maximum number of bodies from `target_systems` allowed in a leaf-level branch
- `multipole_threshold::Float64`: number between 0 and 1 (often denoted theta in [0,1]) controls the accuracy by determining the non-dimensional distance after which multipoles are used; 0 means an infinite distance (no error, high cost), and 1 means barely convergent (high error, low cost)
- `scalar_potential::Bool`: either a `::Bool` or a `::AbstractVector{Bool}` of length `length(target_systems)` indicating whether each system should receive a scalar potential from `source_systems`
- `vector_potential::Bool`: either a `::Bool` or a `::AbstractVector{Bool}` of length `length(target_systems)` indicating whether each system should receive a vector potential from `source_systems`
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
- `save_tree::Bool`: indicates whether or not to save a VTK file for visualizing the octree
- `save_name::String`: name and path of the octree visualization if `save_tree == true`

"""
function fmm!(target_systems, source_systems;
    scalar_potential=true, vector_potential=true, velocity=true, velocity_gradient=true,
    expansion_order=5, leaf_size_source=50, leaf_size_target=50, multipole_threshold=0.4,
    upward_pass=true, horizontal_pass=true, downward_pass=true,
    nearfield=true, farfield=true, self_induced=true,
    influence_matrices=false,
    unsort_source_bodies=true, unsort_target_bodies=true,
    source_shrink_recenter=true, target_shrink_recenter=true,
    save_tree=false, save_name="tree",
    concurrent_direct=false
)
    # check for duplicate systems
    target_systems = wrap_duplicates(target_systems, source_systems)

    # create trees
    source_tree = Tree(source_systems; expansion_order, leaf_size=leaf_size_source, shrink_recenter=source_shrink_recenter)
    target_tree = Tree(target_systems; expansion_order, leaf_size=leaf_size_target, shrink_recenter=target_shrink_recenter)

    # perform fmm
    m2l_list, direct_list, derivatives_switches = fmm!(target_tree, target_systems, source_tree, source_systems;
        scalar_potential, vector_potential, velocity, velocity_gradient,
        multipole_threshold,
        reset_source_tree=false, reset_target_tree=false,
        upward_pass, horizontal_pass, downward_pass,
        influence_matrices,
        nearfield, farfield, self_induced,
        unsort_source_bodies, unsort_target_bodies,
        concurrent_direct
    )

    # visualize
    save_tree && (visualize(save_name, systems, tree))

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
- `vector_potential::Bool`: either a `::Bool` or a `::AbstractVector{Bool}` of length `length(systems)` indicating whether each system should receive a vector potential from `source_systems`
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

"""
function fmm!(systems;
    scalar_potential=true, vector_potential=true, velocity=true, velocity_gradient=true,
    expansion_order=5, leaf_size=50, multipole_threshold=0.4,
    upward_pass=true, horizontal_pass=true, downward_pass=true,
    nearfield=true, farfield=true, self_induced=true,
    influence_matrices=false,
    unsort_bodies=true, shrink_recenter=true,
    save_tree=false, save_name="tree",
    concurrent_direct=false
)
    # create tree
    tree = Tree(systems; expansion_order, leaf_size, shrink_recenter)

    # perform fmm
    m2l_list, direct_list, derivatives_switches = fmm!(tree, systems;
        scalar_potential, vector_potential, velocity, velocity_gradient,
        multipole_threshold, reset_tree=false,
        upward_pass, horizontal_pass, downward_pass,
        nearfield, farfield, self_induced,
        influence_matrices,
        unsort_bodies, concurrent_direct
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
- `vector_potential::Bool`: either a `::Bool` or a `::AbstractVector{Bool}` of length `length(systems)` indicating whether each system should receive a vector potential from `source_systems`
- `velocity::Bool`: either a `::Bool` or a `::AbstractVector{Bool}` of length `length(systems)` indicating whether each system should receive a velocity from `source_systems`
- `velocity_gradient::Bool`: either a `::Bool` or a `::AbstractVector{Bool}` of length `length(systems)` indicating whether each system should receive a velocity gradient from `source_systems`
- `upward_pass::Bool`: whether or not to form the multipole expansions from source bodies and translate them upward in the source tree
- `horizontal_pass::Bool`: whether or not to transform multipole expansions from the source tree into local expansions in the target tree
- `downward_pass::Bool`: whether or not to translate local expansions down to the leaf level of the target tree and evaluate them
- `nearfield::Bool`: indicates whether near-field (computed without multipoles) interactions should be included in the direct_list
- `farfield::Bool`: indicates whether far-field (computed with multipoles) interactions should be included in the m2l_list
- `self_induced::Bool`: indicates whether to include the interactions of each leaf-level branch on itself in the direct_list
- `unsort_bodies::Bool`: indicates whether or not to undo the sort operation used to generate the octree for `systems`

"""
function fmm!(tree::Tree, systems;
    scalar_potential=true, vector_potential=true, velocity=true, velocity_gradient=true,
    multipole_threshold=0.4, reset_tree=true,
    upward_pass=true, horizontal_pass=true, downward_pass=true,
    nearfield=true, farfield=true, self_induced=true,
    influence_matrices=false,
    unsort_bodies=true,
    concurrent_direct=false
)

    # assemble derivatives switch
    derivatives_switches = DerivativesSwitch(scalar_potential, vector_potential, velocity, velocity_gradient, systems)

    # create interaction lists
    m2l_list, direct_list = build_interaction_lists(tree.branches, tree.branches, tree.leaf_index, multipole_threshold, farfield, nearfield, self_induced)
    influence_matrices && ( direct_list = InteractionList(direct_list, systems, tree, systems, tree, derivatives_switches) )

    # run fmm
    fmm!(tree, systems, m2l_list, direct_list, derivatives_switches;
        reset_tree,
        nearfield, upward_pass, horizontal_pass, downward_pass,
        unsort_bodies, concurrent_direct
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
- `vector_potential::Bool`: either a `::Bool` or a `::AbstractVector{Bool}` of length `length(target_systems)` indicating whether each system should receive a vector potential from `source_systems`
- `velocity::Bool`: either a `::Bool` or a `::AbstractVector{Bool}` of length `length(target_systems)` indicating whether each system should receive a velocity from `source_systems`
- `velocity_gradient::Bool`: either a `::Bool` or a `::AbstractVector{Bool}` of length `length(target_systems)` indicating whether each system should receive a velocity gradient from `source_systems`
- `upward_pass::Bool`: whether or not to form the multipole expansions from source bodies and translate them upward in the source tree
- `horizontal_pass::Bool`: whether or not to transform multipole expansions from the source tree into local expansions in the target tree
- `downward_pass::Bool`: whether or not to translate local expansions down to the leaf level of the target tree and evaluate them
- `nearfield::Bool`: indicates whether near-field (computed without multipoles) interactions should be included in the direct_list
- `farfield::Bool`: indicates whether far-field (computed with multipoles) interactions should be included in the m2l_list
- `self_induced::Bool`: indicates whether to include the interactions of each leaf-level branch on itself in the direct_list
- `unsort_source_bodies::Bool`: indicates whether or not to undo the sort operation used to generate the octree for `source_systems`
- `unsort_target_bodies::Bool`: indicates whether or not to undo the sort operation used to generate the octree for `target_systems`

"""
function fmm!(target_tree::Tree, target_systems, source_tree::Tree, source_systems;
    scalar_potential=true, vector_potential=true, velocity=true, velocity_gradient=true,
    multipole_threshold=0.4,
    reset_source_tree=true, reset_target_tree=true,
    upward_pass=true, horizontal_pass=true, downward_pass=true,
    nearfield=true, farfield=true, self_induced=true,
    influence_matrices=false,
    unsort_source_bodies=true, unsort_target_bodies=true,
    concurrent_direct=false
)

    # assemble derivatives switch
    derivatives_switches = DerivativesSwitch(scalar_potential, vector_potential, velocity, velocity_gradient, target_systems)

    # create interaction lists
    m2l_list, direct_list = build_interaction_lists(target_tree.branches, source_tree.branches, source_tree.leaf_index, multipole_threshold, farfield, nearfield, self_induced)
    influence_matrices && ( direct_list = InteractionList(direct_list, target_systems, target_tree, source_systems, source_tree, derivatives_switches) )

    # run fmm
    fmm!(target_tree, target_systems, source_tree, source_systems, m2l_list, direct_list, derivatives_switches;
        reset_source_tree, reset_target_tree,
        nearfield, upward_pass, horizontal_pass, downward_pass,
        unsort_source_bodies, unsort_target_bodies, concurrent_direct
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
function fmm!(tree::Tree, systems, m2l_list, direct_list, derivatives_switches;
    reset_tree=true,
    nearfield=true, upward_pass=true, horizontal_pass=true, downward_pass=true,
    unsort_bodies=true, concurrent_direct=false
)

    fmm!(tree, systems, tree, systems, m2l_list, direct_list, derivatives_switches;
        reset_source_tree=reset_tree, reset_target_tree=false,
        nearfield, upward_pass, horizontal_pass, downward_pass,
        unsort_source_bodies=unsort_bodies, unsort_target_bodies=false,
        concurrent_direct
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
- `vector_potential::Bool`: either a `::Bool` or a `::AbstractVector{Bool}` of length `length(target_systems)` indicating whether each system should receive a vector potential from `source_systems`
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
function fmm!(target_tree::Tree, target_systems, source_tree::Tree, source_systems, m2l_list, direct_list, derivatives_switches;
    reset_source_tree=true, reset_target_tree=true,
    nearfield=true, upward_pass=true, horizontal_pass=true, downward_pass=true,
    unsort_source_bodies=true, unsort_target_bodies=true,
    concurrent_direct=false
)
    # check if systems are empty
    n_sources = get_n_bodies(source_systems)
    n_targets = get_n_bodies(target_systems)

    if n_sources > 0 && n_targets > 0

        # reset multipole/local expansions
        reset_target_tree && (reset_expansions!(source_tree))
        reset_source_tree && (reset_expansions!(source_tree))

        # run FMM
        if Threads.nthreads() == 1
            nearfield && nearfield_singlethread!(target_systems, target_tree, derivatives_switches, source_systems, source_tree, direct_list)
            upward_pass && upward_pass_singlethread!(source_tree.branches, source_systems, source_tree.expansion_order)
            horizontal_pass && horizontal_pass_singlethread!(target_tree.branches, source_tree.branches, m2l_list, source_tree.expansion_order)
            downward_pass && downward_pass_singlethread!(target_tree.branches, target_systems, derivatives_switches, target_tree.expansion_order)
        else # multithread
            # @assert !(typeof(direct_list) <: InteractionList) "`InteractionList` objects do not yet support multithreading"
            @sync begin
                if nearfield && length(direct_list) > 0
                    if concurrent_direct || typeof(direct_list) <: InteractionList
                        Threads.@spawn nearfield_singlethread!(target_systems, target_tree, derivatives_switches, source_systems, source_tree, direct_list)
                    else
                        nearfield_multithread!(target_systems, target_tree.branches, derivatives_switches, source_systems, source_tree.branches, direct_list)
                    end
                end
                upward_pass && upward_pass_multithread!(source_tree.branches, source_systems, source_tree.expansion_order, source_tree.levels_index, source_tree.leaf_index, concurrent_direct)
                horizontal_pass && length(m2l_list) > 0 && horizontal_pass_multithread!(target_tree.branches, source_tree.branches, m2l_list, source_tree.expansion_order, concurrent_direct)
            end
            downward_pass && downward_pass_multithread!(target_tree.branches, target_systems, derivatives_switches, target_tree.expansion_order, target_tree.levels_index, target_tree.leaf_index)
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
