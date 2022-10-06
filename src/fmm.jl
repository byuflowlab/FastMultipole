# struct Options{TF} where TF
#     expansion_order::Int8
#     n_per_branch::Int16
#     threshold_squared::TF
# end

# function Options(expansion_order, n_per_branch, threshold_squared)
#     Options(Int8(expansion_order), Int16(n_per_branch), threshold_squared)
# end

include("cartesian.jl")
include("spherical.jl")

function P2P!(i_target, j_source, tree, elements)
    target_branch = tree.branches[i_target]
    source_branch = tree.branches[j_source]
    for source_element in elements[tree.indices[source_branch.first_element:source_branch.first_element+source_branch.n_elements-1]]
        for target_element in elements[tree.indices[target_branch.first_element:target_branch.first_element+target_branch.n_elements-1]]
            dx = get_x(target_element) - get_x(source_element)
            kernel!(target_element, source_element)
        end
    end
end

function upward_pass!(tree::Tree, elements, basis::Basis)
    upward_pass!(1, tree, elements, basis)
end

function upward_pass!(i_branch, tree, elements, basis::Basis)
    # recursively iterate through branches
    branch = tree.branches[i_branch]
    for i_child in branch.first_branch:branch.first_branch + branch.n_branches-1
        upward_pass!(i_child, tree, elements, basis)
    end

    # perform P2M (leaf level) or M2M (not leaf level) translations
    if branch.first_branch == -1 # no child branches
        P2M!(i_branch, tree, elements, basis)
    else
        M2M!(i_branch, tree, basis)
    end
end

function horizontal_pass!(tree, elements, derivatives, theta, basis::Basis)
    horizontal_pass!(1, 1, tree, elements, derivatives, theta, basis)
end

function horizontal_pass!(i_target, j_source, tree, elements, derivatives, theta, basis::Basis)
    branches = tree.branches
    source_branch = branches[j_source]
    target_branch = branches[i_target]
    spacing = source_branch.center - target_branch.center
    spacing_squared = spacing' * spacing
    threshold_squared = (target_branch.radius + source_branch.radius) * (target_branch.radius + source_branch.radius) * theta # theta is the number of radii squared
    if spacing_squared >= threshold_squared # meet separation criteria
        M2L!(i_target, j_source, tree, elements, derivatives, basis)
    elseif source_branch.first_branch == target_branch.first_branch == -1 # both leaves
        P2P!(i_target, j_source, tree, elements)
    elseif source_branch.first_branch == -1 || (target_branch.radius >= source_branch.radius && target_branch.first_branch != -1)
        for i_child in target_branch.first_branch:target_branch.first_branch + target_branch.n_branches - 1
            horizontal_pass!(i_child, j_source, tree, elements, derivatives, theta, basis)
        end
    else
        for j_child in source_branch.first_branch:source_branch.first_branch + source_branch.n_branches - 1
            horizontal_pass!(i_target, j_child, tree, elements, derivatives, theta, basis)
        end
    end
end

function downward_pass!(tree, elements, basis::Basis)
    downward_pass!(1, tree, elements, basis)
end

function downward_pass!(j_source, tree, elements, basis::Basis)
    # expose branch
    branch = tree.branches[j_source]

    # if a leaf, perform L2P on
    if branch.first_branch == -1 # leaf
        L2P!(j_source, tree, elements, basis)
    else # not a leaf, so perform L2L! and recurse
        L2L!(j_source, tree, basis)
        for i_child in branch.first_branch:branch.first_branch + branch.n_branches - 1
            downward_pass!(i_child, tree, elements, basis)
        end
    end
end

function fmm!(tree::Tree, elements, derivatives, theta, basis::Basis; reset_tree=true)
    if reset_tree; reset_expansions!(tree); end
    upward_pass!(tree, elements, basis)
    horizontal_pass!(tree, elements, derivatives, theta, basis)
    downward_pass!(tree, elements, basis)
end

function fmm!(elements, derivatives, expansion_order, n_per_branch, theta, basis::Basis)
    tree = Tree(elements, basis; expansion_order, n_per_branch)
    fmm!(tree, elements, derivatives, theta, basis)
    return tree
end
