function P2P!(tree, elements, i_target, j_source)
    target_branch = tree.branches[i_target]
    is_target = target_branch.first_body:target_branch.first_body+target_branch.n_bodies-1
    target_positions = view(elements.bodies,i_POSITION,is_target)
    target_potential = view(elements.potential, :, is_target)
    source_branch = tree.branches[j_source]
    is_source = source_branch.first_body:source_branch.first_body+source_branch.n_bodies-1
    sources = view(elements.bodies,:,is_source)
    elements.direct!(target_potential, target_positions, sources)
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
        B2M!(tree, elements, i_branch)
    else
        M2M!(tree, i_branch)
    end
end

function horizontal_pass!(tree, elements, theta)
    horizontal_pass!(tree, elements, 1, 1, theta)
end

function horizontal_pass!(tree, elements, i_target, j_source, theta)
    branches = tree.branches
    source_branch = branches[j_source]
    target_branch = branches[i_target]
    spacing = source_branch.center - target_branch.center
    spacing_squared = spacing' * spacing
    threshold_squared = (target_branch.radius + source_branch.radius) * (target_branch.radius + source_branch.radius) * theta # theta is the number of radii squared
    if spacing_squared >= threshold_squared # meet separation criteria
        M2L!(tree, i_target, j_source)
    elseif source_branch.first_branch == target_branch.first_branch == -1 # both leaves
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
        L2B!(tree, elements, j_source)
    else # not a leaf, so perform L2L! and recurse
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

function fmm!(elements, expansion_order, n_per_branch, theta)
    tree = Tree(elements; expansion_order, n_per_branch)
    fmm!(tree, elements, theta)
    return tree
end
