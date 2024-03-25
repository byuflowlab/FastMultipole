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

function body_2_multipole_multithread!(branches, systems::Tuple, expansion_order::Val{P}, leaf_index) where P
    ## load balance
    n_threads = Threads.nthreads()
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

function body_2_multipole_multithread!(branches, system, expansion_order::Val{P}, leaf_index) where P
    ## load balance
    n_threads = Threads.nthreads()
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

function translate_multipoles_multithread!(branches, expansion_order::Val{P}, levels_index) where P
    # initialize memory
    n_threads = Threads.nthreads()
    
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

function upward_pass_multithread!(branches, systems, expansion_order, levels_index, leaf_index)
    # create multipole expansions
    body_2_multipole_multithread!(branches, systems, expansion_order, leaf_index)

    # m2m translation
    translate_multipoles_multithread!(branches, expansion_order, levels_index)
end

#####
##### horizontal pass
#####
function nearfield_singlethread!(target_system, target_branches, derivatives_switch, source_system, source_branches, direct_list)
    for (i_target, j_source) in direct_list
        P2P!(target_system, target_branches[i_target], derivatives_switch, source_system, source_branches[j_source])
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

function horizontal_pass_multithread!(target_branches, source_branches::Vector{<:Branch{TF}}, m2l_list, expansion_order::Val{P}) where {TF,P}
    # number of translations per thread
    n_threads = Threads.nthreads()
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

@inline function preallocate_l2b(float_type, expansion_type, expansion_order::Val{P}) where P
    vector_potential = zeros(float_type,3)
    potential_jacobian = zeros(float_type,3,4)
    potential_hessian = zeros(float_type,3,3,4)
    derivative_harmonics = zeros(expansion_type, 2, ((P+1) * (P+2)) >> 1)
    derivative_harmonics_theta = zeros(expansion_type, 2, ((P+1) * (P+2)) >> 1)
    derivative_harmonics_theta_2 = zeros(expansion_type, 2, ((P+1) * (P+2)) >> 1)
    workspace = zeros(float_type,3,4)
    return vector_potential, potential_jacobian, potential_hessian, derivative_harmonics, derivative_harmonics_theta, derivative_harmonics_theta_2, workspace
end

function downward_pass_singlethread!(branches, systems, derivatives_switches, expansion_order::Val{P}) where P
    regular_harmonics = zeros(eltype(branches[1].multipole_expansion), 2, (P+1)*(P+1))
    # L = zeros(eltype(branches[1].multipole_expansion),4)
    vector_potential, potential_jacobian, potential_hessian, derivative_harmonics, derivative_harmonics_theta, derivative_harmonics_theta_2, workspace = preallocate_l2b(eltype(branches[1]), eltype(branches[1].multipole_expansion), expansion_order)
    for branch in branches
        if branch.n_branches == 0 # leaf level
            L2B!(systems, branch, derivatives_switches, expansion_order, vector_potential, potential_jacobian, potential_hessian, 
                derivative_harmonics, derivative_harmonics_theta, derivative_harmonics_theta_2, workspace)
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

    # preallocate containers
    containers = preallocate_l2b(eltype(branches[1]), eltype(branches[1].multipole_expansion), expansion_order, n_threads)

    # spread remainder across rem chunks
    Threads.@threads for i_thread in eachindex(assignments)
        assignment = assignments[i_thread]
        vector_potential, potential_jacobian, potential_hessian, derivative_harmonics, derivative_harmonics_theta, derivative_harmonics_theta_2, workspace = containers[i_thread]
        for i_leaf in view(leaf_index,assignment)
            leaf = branches[i_leaf]
            L2B!(systems, leaf, derivatives_switches, expansion_order, vector_potential, potential_jacobian, potential_hessian, derivative_harmonics, derivative_harmonics_theta, derivative_harmonics_theta_2, workspace)
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
function build_interaction_lists(target_branches, source_branches, multipole_acceptance_criterion, farfield, nearfield, self_induced)
    m2l_list = Vector{SVector{2,Int32}}(undef,0)
    direct_list = Vector{SVector{2,Int32}}(undef,0)
    build_interaction_lists!(m2l_list, direct_list, 1, 1, target_branches, source_branches, multipole_acceptance_criterion, farfield, nearfield, self_induced)
    return m2l_list, direct_list
end

function build_interaction_lists!(m2l_list, direct_list, i_target, j_source, target_branches, source_branches, multipole_acceptance_criterion, farfield, nearfield, self_induced)
    source_branch = source_branches[j_source]
    target_branch = target_branches[i_target]

    spacing = source_branch.center - target_branch.center
    center_spacing_squared = spacing[1]*spacing[1] + spacing[2]*spacing[2] + spacing[3]*spacing[3]
    summed_radii_squared = target_branch.radius + source_branch.radius
    summed_radii_squared *= summed_radii_squared
    if center_spacing_squared * multipole_acceptance_criterion * multipole_acceptance_criterion >= summed_radii_squared && farfield # meet M2L criteria
        push!(m2l_list, SVector{2}(i_target, j_source))
    elseif source_branch.n_branches == target_branch.n_branches == 0 && nearfield && (i_target!=j_source || self_induced) # both leaves
        push!(direct_list, SVector{2}(i_target, j_source))
    elseif source_branch.n_branches == 0 || (target_branch.radius >= source_branch.radius && target_branch.n_branches != 0) # source is a leaf OR target is not a leaf and is bigger or the same size
        for i_child in target_branch.branch_index
            build_interaction_lists!(m2l_list, direct_list, i_child, j_source, target_branches, source_branches, multipole_acceptance_criterion, farfield, nearfield, self_induced)
        end
    else # source is not a leaf AND target is a leaf or is smaller
        for j_child in source_branch.branch_index
            build_interaction_lists!(m2l_list, direct_list, i_target, j_child, target_branches, source_branches, multipole_acceptance_criterion, farfield, nearfield, self_induced)
        end
    end
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
- `n_per_branch_source::Int`: maximum number of bodies from `source_systems` allowed in a leaf-level branch
- `n_per_branch_target::Int`: maximum number of bodies from `target_systems` allowed in a leaf-level branch
- `multipole_acceptance_criterion::Float64`: number between 0 and 1 (often denoted theta in [0,1]) controls the accuracy by determining the non-dimensional distance after which multipoles are used; 0 means an infinite distance (no error, high cost), and 1 means barely convergent (high error, low cost)
- `scalar_potential::Bool`: either a `::Bool` or a `::AbstractVector{Bool}` of length `length(target_systems)` indicating whether each system should receive a scalar potential from `source_systems`
- `vector_potential::Bool`: either a `::Bool` or a `::AbstractVector{Bool}` of length `length(target_systems)` indicating whether each system should receive a vector potential from `source_systems`
- `velocity::Bool`: either a `::Bool` or a `::AbstractVector{Bool}` of length `length(target_systems)` indicating whether each system should receive a velocity from `source_systems`
- `velocity_gradient::Bool`: either a `::Bool` or a `::AbstractVector{Bool}` of length `length(target_systems)` indicating whether each system should receive a velocity gradient from `source_systems`
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
    expansion_order=5, n_per_branch_source=50, n_per_branch_target=50, multipole_acceptance_criterion=0.4, 
    nearfield=true, farfield=true, self_induced=true, 
    unsort_source_bodies=true, unsort_target_bodies=true, 
    source_shrink_recenter=true, target_shrink_recenter=true, 
    save_tree=false, save_name="tree"
)
    # check for duplicate systems
    target_systems = wrap_duplicates(target_systems, source_systems)

    # create trees
    source_tree = Tree(source_systems; expansion_order, n_per_branch=n_per_branch_source, shrink_recenter=source_shrink_recenter)
    target_tree = Tree(target_systems; expansion_order, n_per_branch=n_per_branch_target, shrink_recenter=target_shrink_recenter)

    # perform fmm
    fmm!(target_tree, target_systems, source_tree, source_systems; 
        scalar_potential, vector_potential, velocity, velocity_gradient, 
        multipole_acceptance_criterion, 
        reset_source_tree=false, reset_target_tree=false, 
        nearfield, farfield, self_induced, 
        unsort_source_bodies, unsort_target_bodies
    )

    # visualize
    save_tree && (visualize(save_name, systems, tree))
    
    return source_tree, target_tree
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
- `n_per_branch::Int`: maximum number of bodies from `systems` allowed in a leaf-level branch
- `multipole_acceptance_criterion::Float64`: number between 0 and 1 (often denoted theta in [0,1]) controls the accuracy by determining the non-dimensional distance after which multipoles are used; 0 means an infinite distance (no error, high cost), and 1 means barely convergent (high error, low cost)
- `scalar_potential::Bool`: either a `::Bool` or a `::AbstractVector{Bool}` of length `length(systems)` indicating whether each system should receive a scalar potential from `source_systems`
- `vector_potential::Bool`: either a `::Bool` or a `::AbstractVector{Bool}` of length `length(systems)` indicating whether each system should receive a vector potential from `source_systems`
- `velocity::Bool`: either a `::Bool` or a `::AbstractVector{Bool}` of length `length(systems)` indicating whether each system should receive a velocity from `source_systems`
- `velocity_gradient::Bool`: either a `::Bool` or a `::AbstractVector{Bool}` of length `length(systems)` indicating whether each system should receive a velocity gradient from `source_systems`
- `nearfield::Bool`: indicates whether near-field (comuted without multipoles) interactions should be included
- `farfield::Bool`: indicates whether far-field (comuted with multipoles) interactions should be included
- `self_induced::Bool`: indicates whether to include the interactions of each leaf-level branch on itself
- `unsort_bodies::Bool`: indicates whether or not to undo the sort operation used to generate the octree for `systems`
- `shink_recenter::Bool`: indicates whether or not to resize branches for the octree after it is created to increase computational efficiency
- `save_tree::Bool`: indicates whether or not to save a VTK file for visualizing the octree
- `save_name::String`: name and path of the octree visualization if `save_tree == true`
    
"""
function fmm!(systems; 
    scalar_potential=true, vector_potential=true, velocity=true, velocity_gradient=true,
    expansion_order=5, n_per_branch=50, multipole_acceptance_criterion=0.4, 
    nearfield=true, farfield=true, self_induced=true, 
    unsort_bodies=true, shrink_recenter=true, 
    save_tree=false, save_name="tree"
)
    # create tree
    tree = Tree(systems; expansion_order, n_per_branch, shrink_recenter)
    
    # perform fmm
    fmm!(tree, systems;
        scalar_potential, vector_potential, velocity, velocity_gradient,
        multipole_acceptance_criterion, reset_tree=false, 
        nearfield, farfield, self_induced, 
        unsort_bodies
    )
    
    # visualize
    save_tree && (visualize(save_name, systems, tree))

    return tree
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

- `multipole_acceptance_criterion::Float64`: number between 0 and 1 (often denoted theta in [0,1]) controls the accuracy by determining the non-dimensional distance after which multipoles are used; 0 means an infinite distance (no error, high cost), and 1 means barely convergent (high error, low cost)
- `scalar_potential::Bool`: either a `::Bool` or a `::AbstractVector{Bool}` of length `length(systems)` indicating whether each system should receive a scalar potential from `source_systems`
- `vector_potential::Bool`: either a `::Bool` or a `::AbstractVector{Bool}` of length `length(systems)` indicating whether each system should receive a vector potential from `source_systems`
- `velocity::Bool`: either a `::Bool` or a `::AbstractVector{Bool}` of length `length(systems)` indicating whether each system should receive a velocity from `source_systems`
- `velocity_gradient::Bool`: either a `::Bool` or a `::AbstractVector{Bool}` of length `length(systems)` indicating whether each system should receive a velocity gradient from `source_systems`
- `nearfield::Bool`: indicates whether near-field (comuted without multipoles) interactions should be included
- `farfield::Bool`: indicates whether far-field (comuted with multipoles) interactions should be included
- `self_induced::Bool`: indicates whether to include the interactions of each leaf-level branch on itself
- `unsort_bodies::Bool`: indicates whether or not to undo the sort operation used to generate the octree for `systems`
    
"""
function fmm!(tree::Tree, systems; 
    scalar_potential=true, vector_potential=true, velocity=true, velocity_gradient=true,
    multipole_acceptance_criterion=0.4, reset_tree=true, 
    nearfield=true, farfield=true, self_induced=true, 
    unsort_bodies=true
)
    fmm!(tree, systems, tree, systems;
        scalar_potential, vector_potential, velocity, velocity_gradient,
        multipole_acceptance_criterion, reset_source_tree=reset_tree, reset_target_tree=false, 
        nearfield, farfield, self_induced, 
        unsort_source_bodies=unsort_bodies, unsort_target_bodies=false
    )
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

- `multipole_acceptance_criterion::Float64`: number between 0 and 1 (often denoted theta in [0,1]) controls the accuracy by determining the non-dimensional distance after which multipoles are used; 0 means an infinite distance (no error, high cost), and 1 means barely convergent (high error, low cost)
- `scalar_potential::Bool`: either a `::Bool` or a `::AbstractVector{Bool}` of length `length(target_systems)` indicating whether each system should receive a scalar potential from `source_systems`
- `vector_potential::Bool`: either a `::Bool` or a `::AbstractVector{Bool}` of length `length(target_systems)` indicating whether each system should receive a vector potential from `source_systems`
- `velocity::Bool`: either a `::Bool` or a `::AbstractVector{Bool}` of length `length(target_systems)` indicating whether each system should receive a velocity from `source_systems`
- `velocity_gradient::Bool`: either a `::Bool` or a `::AbstractVector{Bool}` of length `length(target_systems)` indicating whether each system should receive a velocity gradient from `source_systems`
- `nearfield::Bool`: indicates whether near-field (comuted without multipoles) interactions should be included
- `farfield::Bool`: indicates whether far-field (comuted with multipoles) interactions should be included
- `self_induced::Bool`: indicates whether to include the interactions of each leaf-level branch on itself
- `unsort_source_bodies::Bool`: indicates whether or not to undo the sort operation used to generate the octree for `source_systems`
- `unsort_target_bodies::Bool`: indicates whether or not to undo the sort operation used to generate the octree for `target_systems`

"""
function fmm!(target_tree::Tree, target_systems, source_tree::Tree, source_systems;
    scalar_potential=true, vector_potential=true, velocity=true, velocity_gradient=true,
    multipole_acceptance_criterion=0.4, 
    reset_source_tree=true, reset_target_tree=true, 
    nearfield=true, farfield=true, self_induced=true, 
    unsort_source_bodies=true, unsort_target_bodies=true
)
    # reset multipole/local expansions
    reset_target_tree && (reset_expansions!(source_tree))
    reset_source_tree && (reset_expansions!(source_tree))

    # create interaction lists
    m2l_list, direct_list = build_interaction_lists(target_tree.branches, source_tree.branches, multipole_acceptance_criterion, farfield, nearfield, self_induced)

    # assemble derivatives switch
    derivatives_switch = DerivativesSwitch(scalar_potential, vector_potential, velocity, velocity_gradient, target_systems)

    # run FMM
    if Threads.nthreads() == 1
        nearfield && (nearfield_singlethread!(target_systems, target_tree.branches, derivatives_switch, source_systems, source_tree.branches, direct_list))
        if farfield
            upward_pass_singlethread!(source_tree.branches, source_systems, source_tree.expansion_order)
            horizontal_pass_singlethread!(target_tree.branches, source_tree.branches, m2l_list, source_tree.expansion_order)
            downward_pass_singlethread!(target_tree.branches, target_systems, derivatives_switch, target_tree.expansion_order)
        end
    else # multithread
        # println("nearfield")
        nearfield && length(direct_list) > 0 && (nearfield_multithread!(target_systems, target_tree.branches, derivatives_switch, source_systems, source_tree.branches, direct_list))
        # @time nearfield && length(direct_list) > 0 && (nearfield_multithread!(target_systems, target_tree.branches, source_systems, source_tree.branches, direct_list))
        if farfield
            # println("upward pass")
            upward_pass_multithread!(source_tree.branches, source_systems, source_tree.expansion_order, source_tree.levels_index, source_tree.leaf_index)
            # @time upward_pass_multithread!(source_tree.branches, source_systems, source_tree.expansion_order, source_tree.levels_index, source_tree.leaf_index)
            # println("horizontal pass")
            length(m2l_list) > 0 && (horizontal_pass_multithread!(target_tree.branches, source_tree.branches, m2l_list, source_tree.expansion_order))
            # @time length(m2l_list) > 0 && (horizontal_pass_multithread!(target_tree.branches, source_tree.branches, m2l_list, source_tree.expansion_order))
            # println("downward pass")
            downward_pass_multithread!(target_tree.branches, target_systems, derivatives_switch, target_tree.expansion_order, target_tree.levels_index, target_tree.leaf_index)
            # @time downward_pass_multithread!(target_tree.branches, target_systems, target_tree.expansion_order, target_tree.levels_index, target_tree.leaf_index)
        end
    end

    # unsort bodies
    unsort_target_bodies && (unsort!(target_systems, target_tree))
    unsort_source_bodies && (unsort!(source_systems, source_tree))
end
