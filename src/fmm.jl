function P2P!(system, target_branch::SingleBranch, source_branch::SingleBranch)
    target_indices = target_branch.first_body:target_branch.first_body+target_branch.n_bodies-1
    source_indices = source_branch.first_body:source_branch.first_body+source_branch.n_bodies-1
    direct!(system, target_indices, system, source_indices)
end

function P2P!(systems, target_branch::MultiBranch, source_branch::MultiBranch)
    for (i,source_system) in enumerate(systems)
        source_indices = source_branch.first_body[i]:source_branch.first_body[i]+source_branch.n_bodies[i]-1
        for (j,target_system) in enumerate(systems)
            target_indices = target_branch.first_body[j]:target_branch.first_body[j]+target_branch.n_bodies[j]-1
            direct!(target_system, target_indices, source_system, source_indices)
        end
    end
end

function P2P!(target_system, target_branch::SingleBranch, source_system, source_branch::SingleBranch)
    source_indices = source_branch.first_body:source_branch.first_body+source_branch.n_bodies-1
    target_indices = target_branch.first_body:target_branch.first_body+target_branch.n_bodies-1
    direct!(target_system, target_indices, source_system, source_indices)
end

function P2P!(target_system, target_branch::SingleBranch, source_systems, source_branch::MultiBranch)
    target_indices = target_branch.first_body[j]:target_branch.first_body[j]+target_branch.n_bodies[j]-1
    for (i,source_system) in enumerate(source_systems)
        source_indices = source_branch.first_body[i]:source_branch.first_body[i]+source_branch.n_bodies[i]-1
        direct!(target_system, target_indices, source_system, source_indices)
    end
end

function P2P!(target_systems, target_branch::MultiBranch, source_system, source_branch::SingleBranch)
    source_indices = source_branch.first_body:source_branch.first_body+source_branch.n_bodies-1
    for (j,target_system) in enumerate(target_systems)
        target_indices = target_branch.first_body[j]:target_branch.first_body[j]+target_branch.n_bodies[j]-1
        direct!(target_system, target_indices, source_system, source_indices)
    end
end

function P2P!(target_systems, target_branch::MultiBranch, source_systems, source_branch::MultiBranch)
    for (i,source_system) in enumerate(source_systems)
        source_indices = source_branch.first_body[i]:source_branch.first_body[i]+source_branch.n_bodies[i]-1
        for (j,target_system) in enumerate(target_systems)
            target_indices = target_branch.first_body[j]:target_branch.first_body[j]+target_branch.n_bodies[j]-1
            direct!(target_system, target_indices, source_system, source_indices)
        end
    end
end

function upward_pass!(tree, system)
    upward_pass!(tree.branches, system, 1, tree.expansion_order)
end

function upward_pass!(branches, system, i_branch, expansion_order)
    branch = branches[i_branch]
    Threads.@threads for i_child in branch.first_branch:branch.first_branch + branch.n_branches-1
        upward_pass!(branches, system, i_child, expansion_order)
    end

    if branch.first_branch == -1
        B2M!(branch, system, expansion_order)
    else
        M2M!(branches, i_branch, expansion_order) # no locks needed, since nothing else will be modifying this branch
    end
end

function horizontal_pass!(tree, system, theta, farfield, nearfield)
    horizontal_pass!(tree.branches, system, 1, 1, theta, farfield, nearfield, tree.expansion_order)
end

function horizontal_pass!(branches, system, i_target, j_source, theta, farfield, nearfield, expansion_order)
    source_branch = branches[j_source]
    target_branch = branches[i_target]

    spacing = source_branch.center - target_branch.center
    center_spacing_squared = spacing[1]*spacing[1] + spacing[2]*spacing[2] + spacing[3]*spacing[3]
    summed_radii_squared = target_branch.radius + source_branch.radius
    summed_radii_squared *= summed_radii_squared
    if farfield && center_spacing_squared * theta * theta >= summed_radii_squared
        Threads.lock(target_branch.lock) do
            M2L!(target_branch, source_branch, expansion_order)
        end
    elseif source_branch.first_branch == target_branch.first_branch == -1 && nearfield # both leaves
        Threads.lock(target_branch.child_lock) do
            P2P!(system, target_branch, source_branch)
        end
    elseif source_branch.first_branch == -1 || (target_branch.radius >= source_branch.radius && target_branch.first_branch != -1)
        Threads.@threads for i_child in target_branch.first_branch:target_branch.first_branch + target_branch.n_branches - 1
            horizontal_pass!(branches, system, i_child, j_source, theta, farfield, nearfield, expansion_order)
        end
    else
        Threads.@threads for j_child in source_branch.first_branch:source_branch.first_branch + source_branch.n_branches - 1
            horizontal_pass!(branches, system, i_target, j_child, theta, farfield, nearfield, expansion_order)
        end
    end
end

function horizontal_pass!(target_tree, target_systems, source_tree, source_systems, theta, farfield, nearfield)
    @assert target_tree.expansion_order == source_tree.expansion_order "source and target trees must use the same expansion order"
    horizontal_pass!(target_tree.branches, target_systems, source_tree.branches, source_systems, 1, 1, theta, farfield, nearfield, target_tree.expansion_order)
end

function horizontal_pass!(target_branches, target_systems, source_branches, source_systems, i_target, j_source, theta, farfield, nearfield, expansion_order)
    source_branch = source_branches[j_source]
    target_branch = target_branches[i_target]

    # determine whether to perform M2L
    spacing = source_branch.center - target_branch.center
    spacing_squared = spacing[1]*spacing[1] + spacing[2]*spacing[2] + spacing[3]*spacing[3]
    threshold_squared = (target_branch.radius + source_branch.radius)
    threshold_squared *= threshold_squared
    if farfield && spacing_squared * theta * theta >= threshold_squared # meet separation criteria; theta is the spacing parameter
        Threads.lock(target_branch.lock) do
            M2L!(target_branch, source_branch, expansion_order)
        end
    elseif source_branch.first_branch == -1 == target_branch.first_branch && nearfield # both leaves
        Threads.lock(target_branch.child_lock) do
            P2P!(target_systems, target_branch, source_systems, source_branch)
        end
    elseif source_branch.first_branch == -1 || (target_branch.radius >= source_branch.radius && target_branch.first_branch != -1)
        Threads.@threads for i_child in target_branch.first_branch:target_branch.first_branch + target_branch.n_branches - 1
            horizontal_pass!(target_branches, target_systems, source_branches, source_systems, i_child, j_source, theta, farfield, nearfield, expansion_order)
        end
    else
        Threads.@threads for j_child in source_branch.first_branch:source_branch.first_branch + source_branch.n_branches - 1
            horizontal_pass!(target_branches, target_systems, source_branches, source_systems, i_target, j_child, theta, farfield, nearfield, expansion_order)
        end
    end
end

function downward_pass!(tree, systems)
    downward_pass!(tree.branches, systems, 1, tree.expansion_order)
end

function downward_pass!(branches, systems, j_source, expansion_order)
    branch = branches[j_source]

    if branch.first_branch == -1 # branch is a leaf
        L2B!(systems, branch, expansion_order)
    else # recurse to child branches until we hit a leaf
        println("L2L!")
        L2L!(branches, j_source, expansion_order)
        Threads.@threads for i_child in branch.first_branch:branch.first_branch + branch.n_branches - 1
            downward_pass!(branches, systems, i_child, expansion_order)
        end
    end
end

#####
##### running FMM
#####
function fmm!(tree::Tree, systems; theta=0.4, reset_tree=true, nearfield=true, farfield=true, unsort_bodies=true)
    # reset multipole/local expansions
    reset_tree && (reset_expansions!(tree))
    
    # run FMM
    farfield && (upward_pass!(tree, systems))
    horizontal_pass!(tree, systems, theta, farfield, nearfield)
    farfield && (downward_pass!(tree, systems))
    
    # unsort bodies
    unsort_bodies && (unsort!(systems, tree))
end

function fmm!(target_tree::Tree, target_systems, source_tree::Tree, source_systems; theta=0.4, reset_source_tree=true, reset_target_tree=true, nearfield=true, farfield=true, unsort_source_bodies=true, unsort_target_bodies=true)
    # reset multipole/local expansions
    reset_target_tree && (reset_expansion!(source_tree))
    reset_source_tree && (reset_expansion!(source_tree))

    # run FMM
    farfield && (upward_pass!(source_tree, source_systems))
    horizontal_pass!(target_tree, target_systems, source_tree, source_systems, theta, farfield, nearfield)
    farfield && (downward_pass!(target_tree, target_systems))

    # unsort bodies
    unsort_source_bodies && (unsort!(source_systems, source_tree))
    unsort_target_bodies && (unsort!(target_systems, target_tree))
end

# function fmm!(tree::Tree, systems, options::Options, targets_index, sources_index; reset_tree=true, nearfield=true, farfield=true, unsort_bodies=true)
#     reset_tree && (reset_expansions!(tree))
#     farfield && (upward_pass!(tree, systems, sources_index))
#     horizontal_pass!(tree, systems, options.theta, targets_index, sources_index, farfield, nearfield)
#     farfield && (downward_pass!(tree, systems, targets_index))
#     unsort_bodies && (unsort!(systems, tree))
# end

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
function fmm!(systems; expansion_order=5, n_per_branch=50, theta=0.4, nearfield=true, farfield=true, unsort_bodies=true, shrinking=true)
    tree = Tree(systems, expansion_order, n_per_branch, shrinking=shrinking)
    fmm!(tree, systems; theta=theta, reset_tree=false, nearfield=nearfield, farfield=farfield, unsort_bodies=unsort_bodies)
    return tree
end

function fmm!(target_systems, source_systems; expansion_order=5, n_per_branch_source=50, n_per_branch_target=50, theta=0.4, nearfield=true, farfield=true, unsort_source_bodies=true, unsort_target_bodies=true, source_shrinking=true, target_shrinking=true)
    source_tree = Tree(source_systems, expansion_order, n_per_branch_source; shrinking=source_shrinking)
    target_tree = Tree(target_systems, expansion_order, n_per_branch_target; shrinking=target_shrinking)
    fmm!(target_tree, target_systems, source_tree, source_systems; theta=theta, reset_source_tree=false, reset_target_tree=false, nearfield=nearfield, farfield=farfield, unsort_source_bodies=unsort_source_bodies, unsort_target_bodies=unsort_target_bodies)
    return source_tree, target_tree
end
