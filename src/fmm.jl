function P2P!(tree, systems::Tuple, i_target, j_source, targets_index, sources_index)
    target_branch = tree.branches[i_target]
    source_branch = tree.branches[j_source]
    for (ii_source,source_system) in enumerate(systems[sources_index])
        i_source = sources_index[ii_source]
        for (ii_target,target_system) in enumerate(systems[targets_index])
            i_target = targets_index[ii_target]
            target_indices = target_branch.first_body[i_target]:target_branch.first_body[i_target]+target_branch.n_bodies[i_target]-1
            source_indices = source_branch.first_body[i_source]:source_branch.first_body[i_source]+source_branch.n_bodies[i_source]-1
            direct!(target_system, target_indices, source_system, source_indices)
        end
    end
end

function upward_pass!(tree::Tree, systems, sources_index)
    upward_pass!(tree, systems, 1, sources_index)
end

function upward_pass!(tree, systems, i_branch, sources_index)
    # recursively iterate through branches
    branch = tree.branches[i_branch]
    Threads.@threads for i_child in branch.first_branch:branch.first_branch + branch.n_branches-1
    # for i_child in branch.first_branch:branch.first_branch + branch.n_branches-1
        # @inbounds upward_pass!(tree, systems, i_child, sources_index)
        upward_pass!(tree, systems, i_child, sources_index)
    end
    
    contains_sources = false
    for i in sources_index; branch.n_bodies[i] > 0 && (contains_sources = true); end

    if contains_sources # contains sources
        # perform P2M (leaf level) or M2M (not leaf level) translations
        @lock branch.lock begin
            if branch.first_branch == -1 # leaf level
                # println("B2M!")
                B2M!(tree, systems, i_branch, sources_index)
            else
                # println("M2M!")
                M2M!(tree, i_branch) # not leaf level
            end
        end
    end
end

function horizontal_pass!(tree, systems, theta, targets_index, sources_index, local_P2P=true)
    horizontal_pass!(tree, systems, 1, 1, theta, targets_index, sources_index, local_P2P)
end

function horizontal_pass!(tree, systems, i_target, j_source, theta, targets_index, sources_index, local_P2P=true)
    branches = tree.branches
    source_branch = branches[j_source]
    target_branch = branches[i_target]
    # check if target branch contains targets AND source branch contains sources
    contains_targets = false
    for i in targets_index; target_branch.n_bodies[i] > 0 && (contains_targets = true); end
    if contains_targets
        contains_sources = false
        for i in sources_index; source_branch.n_bodies[i] > 0 && (contains_sources = true); end
        if contains_sources
            spacing = source_branch.center - target_branch.center
            spacing_squared = spacing' * spacing
            threshold_squared = (target_branch.radius[1] + source_branch.radius[1])^2 * theta # theta is the number of radii squared
            if spacing_squared >= threshold_squared # meet separation criteria
                @lock target_branch.lock M2L!(tree, i_target, j_source)
            elseif source_branch.first_branch == target_branch.first_branch == -1 && (local_P2P || i_target != j_source) # both leaves
                @lock target_branch.child_lock P2P!(tree, systems, i_target, j_source, targets_index, sources_index)
            elseif source_branch.first_branch == -1 || (target_branch.radius >= source_branch.radius && target_branch.first_branch != -1)
                Threads.@threads for i_child in target_branch.first_branch:target_branch.first_branch + target_branch.n_branches - 1
                # for i_child in target_branch.first_branch:target_branch.first_branch + target_branch.n_branches - 1
                    # @inbounds horizontal_pass!(tree, systems, i_child, j_source, theta, targets_index, sources_index, local_P2P)
                    horizontal_pass!(tree, systems, i_child, j_source, theta, targets_index, sources_index, local_P2P)
                end
            else
                Threads.@threads for j_child in source_branch.first_branch:source_branch.first_branch + source_branch.n_branches - 1
                # for j_child in source_branch.first_branch:source_branch.first_branch + source_branch.n_branches - 1
                    # @inbounds horizontal_pass!(tree, systems, i_target, j_child, theta, targets_index, sources_index, local_P2P)
                    horizontal_pass!(tree, systems, i_target, j_child, theta, targets_index, sources_index, local_P2P)
                end
            end
        end
    end
end

function downward_pass!(tree, systems, targets_index)
    downward_pass!(tree, systems, 1, targets_index)
end

function downward_pass!(tree, systems, j_source, targets_index)
    # expose branch
    branch = tree.branches[j_source]

    contains_targets = false
    for i in targets_index; branch.n_bodies[i] > 0 && (contains_targets = true); end

    if contains_targets
        # if a leaf, perform L2P on
        if branch.first_branch == -1 # leaf
            L2B!(tree, systems, j_source, targets_index)
        else # not a leaf, so perform L2L! and recurse
            L2L!(tree, j_source)
            Threads.@threads for i_child in branch.first_branch:branch.first_branch + branch.n_branches - 1
            # for i_child in branch.first_branch:branch.first_branch + branch.n_branches - 1
                # @inbounds downward_pass!(tree, systems, i_child, targets_index)
                downward_pass!(tree, systems, i_child, targets_index)
            end
        end
    end
end

#####
##### running FMM
#####
function fmm!(tree::Tree, systems, options::Options, targets_index, sources_index; reset_tree=true, local_P2P=true, unsort_bodies=true)
    reset_tree && (reset_expansions!(tree))
    upward_pass!(tree, systems, sources_index)
    horizontal_pass!(tree, systems, options.theta, targets_index, sources_index, local_P2P)
    downward_pass!(tree, systems, targets_index)
    unsort_bodies && (unsort!(systems, tree))
end

"""
    fmm!(elements::Tuple{Element1, Element2,...}, options::Options)

Calculates the influence of all scalar and/or vector potential elements using the fast multipole method in spherical harmonics.

# Inputs

- `elements::Tuple{Element1, Element2, ...}`- a tuple of structs, each containing the following members:

    * `bodies::Array{Float64,2}`- 3+4+mxN array containing element positions, strengths, and m other values that must be sorted into the octree
    * `index::Vector{Int32}`- length N; `index[i]` contains the new index of the element that used to occupy the ith position before being sorted into the octree
    * `potential::Array{Float64,2}`- 4+12+36xN array of the potential, Jacobian, and Hessian that are reset every iteration (and don't require sorting)
    * `direct!::Function`- function calculates the direct influence of the body at the specified location
    * `B2M!::Function`- function converts the body's influence into a multipole expansion

Note: this function merely adds to existing potential of its elements to avoid overwriting potential calculated due to other sources.

The user must reset the potential manually.
"""
function fmm!(systems::Tuple, options::Options; optargs...)
    full_index = collect(1:length(systems))
    fmm!(systems, options, full_index, full_index; optargs...)
    return tree
end

function fmm!(tree::Tree, systems::Tuple, options::Options; optargs...)
    full_index = collect(1:length(systems))
    fmm!(tree, systems, options, full_index, full_index; optargs...)
    return tree
end

function fmm!(systems::Tuple, options::Options, target_index, source_index; sort_in_place = Tuple(true for _ in systems), optargs...)
    # if not sorting in place, wrap systems in SortWrapper (default behavior does nothing)
    if false in sort_in_place
        systems_2 = []
        for (i,(in_place,system)) in enumerate(zip(sort_in_place,systems))
            in_place ? push!(systems_2, system) : push!(systems_2, SortWrapper(system,index_list[i]))
        end
        systems = Tuple(systems_2)
    end

    tree = Tree(systems, options, Val{eltype(systems[1][1,POSITION])}())
    fmm!(tree, systems, options, target_index, source_index; optargs...)
    return tree
end