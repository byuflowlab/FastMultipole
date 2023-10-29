function P2P!(system, target_branch::SingleBranch, source_branch::SingleBranch)
    direct!(system, target_branch.bodies_index, system, source_branch.bodies_index)
end

function P2P!(systems, target_branch::MultiBranch, source_branch::MultiBranch)
    for (i,source_system) in enumerate(systems)
        source_indices = source_branch.bodies_index[i]
        for (j,target_system) in enumerate(systems)
            target_indices = target_branch.bodies_index[j]
            direct!(target_system, target_indices, source_system, source_indices)
        end
    end
end

function P2P!(target_system, target_branch::SingleBranch, source_system, source_branch::SingleBranch)
    direct!(target_system, target_branch.bodies_index, source_system, source_branch.bodies_index)
end

function P2P!(target_system, target_branch::SingleBranch, source_systems, source_branch::MultiBranch)
    target_indices = target_branch.bodies_index
    for (i,source_system) in enumerate(source_systems)
        source_indices = source_branch.bodies_index[i]
        direct!(target_system, target_indices, source_system, source_indices)
    end
end

function P2P!(target_systems, target_branch::MultiBranch, source_system, source_branch::SingleBranch)
    source_indices = source_branch.bodies_index
    for (j,target_system) in enumerate(target_systems)
        target_indices = target_branch.bodies_index[j]
        direct!(target_system, target_indices, source_system, source_indices)
    end
end

function P2P!(target_systems, target_branch::MultiBranch, source_systems, source_branch::MultiBranch)
    for (i,source_system) in enumerate(source_systems)
        source_indices = source_branch.bodies_index[i]
        for (j,target_system) in enumerate(target_systems)
            target_indices = target_branch.bodies_index[j]
            direct!(target_system, target_indices, source_system, source_indices)
        end
    end
end

#####
##### upward pass
#####
function upward_pass!(tree, system)
    upward_pass!(tree.branches, system, 1, tree.expansion_order)
end

function upward_pass!(branches, system, i_branch, expansion_order)
    branch = branches[i_branch]
    for i_child in branch.branch_index
        upward_pass!(branches, system, i_child, expansion_order)
    end

    if branch.n_branches == 0
        B2M!(branch, system, expansion_order)
    else
        M2M!(branches, i_branch, expansion_order) # no locks needed, since nothing else will be modifying this branch
    end
end

#####
##### horizontal pass
#####
function nearfield!(system, branches, direct_list)
    for (i_target, j_source) in direct_list
        P2P!(system, branches[i_target], branches[j_source])
    end
end

function nearfield!(target_system, target_branches, source_system, source_branches, direct_list)
    for (i_target, j_source) in direct_list
        P2P!(target_system, target_branches[i_target], source_system, source_branches[j_source])
    end
end

function horizontal_pass!(branches, m2l_list, expansion_order)
    for (i_target, j_source) in m2l_list
        M2L!(branches[i_target], branches[j_source], expansion_order)
    end
end

function horizontal_pass!(target_branches, source_branches, m2l_list, expansion_order)
    for (i_target, j_source) in m2l_list
        M2L!(target_branches[i_target], source_branches[j_source], expansion_order)
    end
end

#####
##### downward pass
#####
function downward_pass!(tree, systems)
    downward_pass!(tree.branches, systems, 1, tree.expansion_order)
end

function downward_pass!(branches, systems, j_source, expansion_order)
    branch = branches[j_source]

    if branch.n_branches == 0 # branch is a leaf
        L2B!(systems, branch, expansion_order)
    else # recurse to child branches until we hit a leaf
        L2L!(branches, j_source, expansion_order)
        for i_child in branch.branch_index
            downward_pass!(branches, systems, i_child, expansion_order)
        end
    end
end

#####
##### create interaction lists
#####
function build_interaction_lists(branches, theta, farfield, nearfield)
    m2l_list = Vector{SVector{2,Int32}}(undef,0)
    direct_list = Vector{SVector{2,Int32}}(undef,0)
    build_interaction_lists!(m2l_list, direct_list, 1, 1, branches, theta, farfield, nearfield)
    
    return m2l_list, direct_list
end

function build_interaction_lists!(m2l_list, direct_list, i_target, j_source, branches, theta, farfield, nearfield)
    source_branch = branches[j_source]
    target_branch = branches[i_target]

    spacing = source_branch.center - target_branch.center
    center_spacing_squared = spacing[1]*spacing[1] + spacing[2]*spacing[2] + spacing[3]*spacing[3]
    summed_radii_squared = target_branch.radius + source_branch.radius
    summed_radii_squared *= summed_radii_squared
    if center_spacing_squared * theta * theta >= summed_radii_squared && farfield# meet M2L criteria
        push!(m2l_list, SVector{2}(i_target, j_source))
    elseif source_branch.n_branches == target_branch.n_branches == 0 && nearfield # both leaves
        push!(direct_list, SVector{2}(i_target, j_source))
    elseif source_branch.n_branches == 0 || (target_branch.radius >= source_branch.radius && target_branch.n_branches != 0)
        for i_child in target_branch.branch_index
            build_interaction_lists!(m2l_list, direct_list, i_child, j_source, branches, theta, farfield, nearfield)
        end
    else
        for j_child in source_branch.branch_index
            build_interaction_lists!(m2l_list, direct_list, i_target, j_child, branches, theta, farfield, nearfield)
        end
    end
end

function build_interaction_lists(target_branches, source_branches, theta, farfield, nearfield)
    m2l_list = Vector{SVector{2,Int32}}(undef,0)
    direct_list = Vector{SVector{2,Int32}}(undef,0)
    build_interaction_lists!(m2l_list, direct_list, 1, 1, target_branches, source_branches, theta, farfield, nearfield)
    
    return m2l_list, direct_list
end

function build_interaction_lists!(m2l_list, direct_list, i_target, j_source, target_branches, source_branches, theta, farfield, nearfield)
    source_branch = source_branches[j_source]
    target_branch = target_branches[i_target]

    spacing = source_branch.center - target_branch.center
    center_spacing_squared = spacing[1]*spacing[1] + spacing[2]*spacing[2] + spacing[3]*spacing[3]
    summed_radii_squared = target_branch.radius + source_branch.radius
    summed_radii_squared *= summed_radii_squared
    if center_spacing_squared * theta * theta >= summed_radii_squared && farfield# meet M2L criteria
        push!(m2l_list, SVector{2}(i_target, j_source))
    elseif source_branch.n_branches == target_branch.n_branches == 0 && nearfield # both leaves
        push!(direct_list, SVector{2}(i_target, j_source))
    elseif source_branch.n_branches == 0 || (target_branch.radius >= source_branch.radius && target_branch.n_branches != 0)
        for i_child in target_branch.branch_index
            build_interaction_lists!(m2l_list, direct_list, i_child, j_source, target_branches, source_branches, theta, farfield, nearfield)
        end
    else
        for j_child in source_branch.branch_index
            build_interaction_lists!(m2l_list, direct_list, i_target, j_child, target_branches, source_branches, theta, farfield, nearfield)
        end
    end
end

#####
##### running FMM
#####
function fmm!(tree::Tree, systems; theta=0.4, reset_tree=true, nearfield=true, farfield=true, unsort_bodies=true)
    # reset multipole/local expansions
    reset_tree && (reset_expansions!(tree))

    # create interaction lists
    m2l_list, direct_list = build_interaction_lists(tree.branches, theta, farfield, nearfield)
    
    # run FMM
    println("nearfield")
    @time nearfield && (nearfield!(systems, tree.branches, direct_list))
    if farfield
        println("upward pass")
        @time upward_pass!(tree, systems)
        println("horizontal pass")
        @time horizontal_pass!(tree.branches, m2l_list, tree.expansion_order)
        println("downward pass")
        @time downward_pass!(tree, systems)
    end
    
    # unsort bodies
    println("unsort")
    @time unsort_bodies && (unsort!(systems, tree))
end

function fmm!(target_tree::Tree, target_systems, source_tree::Tree, source_systems; theta=0.4, reset_source_tree=true, reset_target_tree=true, nearfield=true, farfield=true, unsort_source_bodies=true, unsort_target_bodies=true)
    # reset multipole/local expansions
    reset_target_tree && (reset_expansions!(source_tree))
    reset_source_tree && (reset_expansions!(source_tree))

    # create interaction lists
    m2l_list, direct_list = build_interaction_lists(target_tree.branches, source_tree.branches, theta, farfield, nearfield)

    # run FMM
    nearfield && (nearfield!(target_systems, target_tree.branches, source_systems, source_tree.branches, direct_list))
    if farfield
        upward_pass!(source_tree, source_systems)
        horizontal_pass!(target_tree.branches, source_tree.branches, m2l_list, source_tree.expansion_order)
        downward_pass!(target_tree, target_systems)
    end

    # unsort bodies
    unsort_source_bodies && (unsort!(source_systems, source_tree))
    unsort_target_bodies && (unsort!(target_systems, target_tree))
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
function fmm!(systems; expansion_order=5, n_per_branch=50, theta=0.4, nearfield=true, farfield=true, unsort_bodies=true, shrink_recenter=true)
    tree = Tree(systems; expansion_order, n_per_branch, shrink_recenter=shrink_recenter)
    fmm!(tree, systems; theta=theta, reset_tree=false, nearfield=nearfield, farfield=farfield, unsort_bodies=unsort_bodies)
    return tree
end

function fmm!(target_systems, source_systems; expansion_order=5, n_per_branch_source=50, n_per_branch_target=50, theta=0.4, nearfield=true, farfield=true, unsort_source_bodies=true, unsort_target_bodies=true, source_shrink_recenter=true, target_shrink_recenter=true)
    source_tree = Tree(source_systems; expansion_order, n_per_branch_source, shrink_recenter=source_shrink_recenter)
    target_tree = Tree(target_systems; expansion_order, n_per_branch_target, shrink_recenter=target_shrink_recenter)
    fmm!(target_tree, target_systems, source_tree, source_systems; theta=theta, reset_source_tree=false, reset_target_tree=false, nearfield=nearfield, farfield=farfield, unsort_source_bodies=unsort_source_bodies, unsort_target_bodies=unsort_target_bodies)
    return source_tree, target_tree
end
