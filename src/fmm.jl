function P2P!(tree, elements_tuple::Tuple, i_target, j_source)
    target_branch = tree.branches[i_target]
    source_branch = tree.branches[j_source]
    for (i_source,source_elements) in enumerate(elements_tuple)
        for (i_target,target_elements) in enumerate(elements_tuple)
            is_target = target_branch.first_body[i_target]:target_branch.first_body[i_target]+target_branch.n_bodies[i_target]-1
            target_positions = target_elements.bodies[i_POSITION,is_target]
            target_potential = view(target_elements.potential, :, is_target)
            is_source = source_branch.first_body[i_source]:source_branch.first_body[i_source]+source_branch.n_bodies[i_source]-1
            sources = source_elements.bodies[:,is_source]
            source_elements.direct!(target_potential, target_positions, sources)
        end
    end
end

function upward_pass!(tree::Tree, elements)
    upward_pass!(tree, elements, 1)
end

function upward_pass!(tree, elements, i_branch)
    # recursively iterate through branches
    branch = tree.branches[i_branch]
    for i_child in branch.first_branch:branch.first_branch + branch.n_branches-1
        upward_pass!(tree, elements, i_child)
    end

    # perform P2M (leaf level) or M2M (not leaf level) translations
    if branch.first_branch == -1 # no child branches
        # println("UP: B2M: $i_branch")
        B2M!(tree, elements, i_branch)
    else
        # println("UP: M2M: $i_branch")
        M2M!(tree, i_branch)
    end
end

function horizontal_pass!(tree, elements::Tuple, theta)
    horizontal_pass!(tree, elements, 1, 1, theta)
end

function horizontal_pass!(tree, elements, i_target, j_source, theta)
    branches = tree.branches
    source_branch = branches[j_source]
    target_branch = branches[i_target]
    spacing = source_branch.center - target_branch.center
    spacing_squared = spacing' * spacing
    threshold_squared = (target_branch.radius + source_branch.radius)^2 * theta # theta is the number of radii squared
    if spacing_squared >= threshold_squared # meet separation criteria
        M2L!(tree, i_target, j_source)
        # println("HP: M2L: $i_target $j_source")
    elseif source_branch.first_branch == target_branch.first_branch == -1 # both leaves
        # println("HP: P2P: $i_target $j_source")
        P2P!(tree, elements, i_target, j_source)
    elseif source_branch.first_branch == -1 || (target_branch.radius >= source_branch.radius && target_branch.first_branch != -1)
        for i_child in target_branch.first_branch:target_branch.first_branch + target_branch.n_branches - 1
            horizontal_pass!(tree, elements, i_child, j_source, theta)
        end
    else
        for j_child in source_branch.first_branch:source_branch.first_branch + source_branch.n_branches - 1
            horizontal_pass!(tree, elements, i_target, j_child, theta)
        end
    end
end

function downward_pass!(tree, elements)
    downward_pass!(tree, elements, 1)
end

function downward_pass!(tree, elements, j_source)
    # expose branch
    branch = tree.branches[j_source]

    # if a leaf, perform L2P on
    if branch.first_branch == -1 # leaf
        # println("DP: L2B: $j_source")
        L2B!(tree, elements, j_source)
    else # not a leaf, so perform L2L! and recurse
        # println("DP: L2L: $j_source")
        L2L!(tree, j_source)
        for i_child in branch.first_branch:branch.first_branch + branch.n_branches - 1
            downward_pass!(tree, elements, i_child)
        end
    end
end

function fmm!(tree::Tree, elements, theta; reset_tree=true)
    if reset_tree; reset_expansions!(tree); end
    upward_pass!(tree, elements)
    horizontal_pass!(tree, elements, theta)
    downward_pass!(tree, elements)
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

Note: this function merely adds to existing potential of its elements to avoid overwriting potential calcualted due to other sources.

The user must reset the potential manually.
"""
function fmm!(elements::Tuple, options::Options)
    tree = Tree(elements, options.expansion_order, options.n_per_branch)
    fmm!(tree, elements, options.theta)
    return tree
end
