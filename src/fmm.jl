# struct Options{TF} where TF
#     expansion_order::Int8
#     n_per_branch::Int16
#     threshold_squared::TF
# end

# function Options(expansion_order, n_per_branch, threshold_squared)
#     Options(Int8(expansion_order), Int16(n_per_branch), threshold_squared)
# end

function P2M!(i_branch, tree, elements)
    branch = tree.branches[i_branch]

    # iterate over coefficients
    i_coeff = 1
    for order in 0:get_expansion_order(tree)
        for i in order:-1:0
            for j in order-i:-1:0
                k = order - i - j # constrained by i and j
                # iterate over elements
                for i_element in tree.indices[branch.first_element:branch.first_element + branch.n_elements-1]
                    element = elements[i_element]
                    dx = branch.center .- get_x(element)
                    branch.multipole_expansion[i_coeff] += ^(dx,i,j,k) * get_q(element)
                end
                branch.multipole_expansion[i_coeff] /= factorial(i) * factorial(j) * factorial(k)
                i_coeff += 1
            end
        end
    end
end

function M2M!(i_branch, tree)
    # expose objects
    branch = tree.branches[i_branch]

    # iterate over children
    for i_child in branch.first_branch:branch.first_branch + branch.n_branches - 1
        child = tree.branches[i_child]
        # get distance vector
        dx = branch.center - child.center
        # iterate over coefficients
        i_coeff = 1
        for order in 0:get_expansion_order(tree)
            for i in order:-1:0
                for j in order-i:-1:0
                    k = order - i - j # constrained by order, i, and j
                    # iterate over child multipole expansion coefficients
                    i_coeff_c = 1
                    for order_c in 0:order
                        for i_c in order_c:-1:0
                            for j_c in order_c-i_c:-1:0
                                k_c = order_c - i_c - j_c
                                if i_c > i || j_c > j || k_c > k; i_coeff_c += 1; continue; end
                                # perform translation of this coefficient
                                branch.multipole_expansion[i_coeff] += ^(dx, i-i_c, j-j_c, k-k_c) * child.multipole_expansion[i_coeff_c] /
                                    factorial(i-i_c) / factorial(j-j_c) / factorial(k-k_c)
                                i_coeff_c += 1
                            end
                        end
                    end
                    i_coeff += 1
                end
            end
        end
    end
end

function M2L!(i_local, j_multipole, tree, elements, derivatives)
    local_branch = tree.branches[i_local]
    multipole_branch = tree.branches[j_multipole]
    dx = local_branch.center - multipole_branch.center
    local_coeff = 1 # iterate over local expansion coefficients
    for order_local in 0:get_expansion_order(tree)
        for i_local in order_local:-1:0
            for j_local in order_local-i_local:-1:0
                k_local = order_local - i_local - j_local
                multipole_coeff = 1
                for order_multipole in 0:get_expansion_order(tree) - order_local
                    for i_multipole in order_multipole:-1:0
                        for j_multipole in order_multipole-i_multipole:-1:0
                            k_multipole = order_multipole - i_multipole - j_multipole
                            # gradient = kernel.potential_derivatives[i_multipole+i_local+1, j_multipole+j_local+1, k_multipole+k_local+1](dx,1,1)
                            gradient_index = ijk_2_index(i_multipole+i_local, j_multipole+j_local, k_multipole+k_local)
                            gradient = derivatives[gradient_index](dx)
                            local_branch.local_expansion[local_coeff] += gradient * multipole_branch.multipole_expansion[multipole_coeff]
                            multipole_coeff += 1
                        end
                    end
                end
                local_coeff += 1
            end
        end
    end
end

# test_d10dz10(x,y,z) = (14175 * (5760 * z^8 * (x^2 + y^2) - 20160 * z^6 * (x^2 + y^2)^2 + 16800 * z^4 * (x^2 + y^2)^3 - 3150 * z^2 * (x^2 + y^2)^4 + 63 * (x^2 + y^2)^5 - 256 * z^10)) / (x^2 + y^2 + z^2)^(21/2)

# test_d6dz6(x,y,z) = (45 * (120 * z^4 * (x^2 + y^2) - 90 * z^2 * (x^2 + y^2)^2 + 5 * (x^2 + y^2)^3 - 16 * z^6))/(x^2 + y^2 + z^2)^(13/2)

# test_d6dx2dy2dz2(x,y,z) = (45* (-2* x^6 + 15* x^4* y^2 + 15* x^2* y^4 + 15* z^4 * (x^2 + y^2) + 15* z^2 * (x^4 - 12* x^2* y^2 + y^4) - 2* y^6 - 2* z^6))/(x^2 + y^2 + z^2)^(13/2)

function L2L!(j_source, tree)
    # expose branch
    branch = tree.branches[j_source]

    # iterate over children
    for i_child in branch.first_branch:branch.first_branch + branch.n_branches - 1
        child = tree.branches[i_child]
        dx = child.center - branch.center

        # iterate over coefficients
        i_coeff_child = 1
        for order in 0:get_expansion_order(tree)
            for i in order:-1:0
                for j in order-i:-1:0
                    k = order - i - j # constrained by order, i, and j

                    # inner summation
                    i_coeff_sum = 1
                    for order_sum in 0:get_expansion_order(tree)
                        for i_sum in order_sum:-1:0
                            for j_sum in order_sum-i_sum:-1:0
                                k_sum = order_sum - i_sum - j_sum

                                # if out of bounds, skip this one
                                if i_sum < i || j_sum < j || k_sum < k; i_coeff_sum += 1; continue; end

                                # translate local expansion from source to child
                                child.local_expansion[i_coeff_child] += ^(dx, i_sum-i, j_sum-j, k_sum-k) * branch.local_expansion[i_coeff_sum] /
                                    (factorial(i_sum - i) * factorial(j_sum - j) * factorial(k_sum - k))

                                i_coeff_sum += 1
                            end
                        end
                    end

                    i_coeff_child += 1
                end
            end
        end
    end
end

"Calculates the potential at all child elements of a branch."
function L2P!(i_branch, tree, elements)
    branch = tree.branches[i_branch]
    for i_element in branch.first_element:branch.first_element + branch.n_elements - 1
        element = elements[tree.indices[i_element]]
        dx = get_x(element) - branch.center
        i_coeff = 1
        for order in 0:get_expansion_order(tree)
            for i in order:-1:0
                for j in order-i:-1:0
                    k = order - i - j
                    element.potential .+= ^(dx, i, j, k) * branch.local_expansion[i_coeff] /
                        (factorial(i) * factorial(j) * factorial(k))
                    i_coeff += 1
                end
            end
        end
    end
end

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

function upward_pass!(tree::Tree, elements)
    upward_pass!(1, tree, elements)
end

function upward_pass!(i_branch, tree, elements)
    # recursively iterate through branches
    branch = tree.branches[i_branch]
    for i_child in branch.first_branch:branch.first_branch + branch.n_branches-1
        upward_pass!(i_child, tree, elements)
    end

    # perform P2M (leaf level) or M2M (not leaf level) translations
    if branch.first_branch == -1 # no child branches
        P2M!(i_branch, tree, elements)
    else
        M2M!(i_branch, tree)
    end
end

function horizontal_pass!(tree, elements, derivatives, theta)
    horizontal_pass!(1, 1, tree, elements, derivatives, theta)
end

function horizontal_pass!(i_target, j_source, tree, elements, derivatives, theta)
    branches = tree.branches
    source_branch = branches[j_source]
    target_branch = branches[i_target]
    spacing = source_branch.center - target_branch.center
    spacing_squared = spacing' * spacing
    threshold_squared = (target_branch.radius + source_branch.radius) * (target_branch.radius + source_branch.radius) * theta # theta is the number of radii squared
    if spacing_squared >= threshold_squared # meet separation criteria
        M2L!(i_target, j_source, tree, elements, derivatives)
    elseif source_branch.first_branch == target_branch.first_branch == -1 # both leaves
        P2P!(i_target, j_source, tree, elements)
    elseif source_branch.first_branch == -1 || (target_branch.radius >= source_branch.radius && target_branch.first_branch != -1)
        for i_child in target_branch.first_branch:target_branch.first_branch + target_branch.n_branches - 1
            horizontal_pass!(i_child, j_source, tree, elements, derivatives, theta)
        end
    else
        for j_child in source_branch.first_branch:source_branch.first_branch + source_branch.n_branches - 1
            horizontal_pass!(i_target, j_child, tree, elements, derivatives, theta)
        end
    end
end

function downward_pass!(tree, elements)
    downward_pass!(1, tree, elements)
end

function downward_pass!(j_source, tree, elements)
    # expose branch
    branch = tree.branches[j_source]

    # if a leaf, perform L2P on
    if branch.first_branch == -1 # leaf
        L2P!(j_source, tree, elements)
    else # not a leaf, so perform L2L! and recurse
        L2L!(j_source, tree)
        for i_child in branch.first_branch:branch.first_branch + branch.n_branches - 1
            downward_pass!(i_child, tree, elements)
        end
    end
end

function fmm!(tree::Tree, elements, derivatives, theta; reset_tree=true, reset_elements=true)
    if reset_elements; reset_elements!(elements); end
    if reset_tree; reset_expansions!(tree); end
    upward_pass!(tree, elements)
    horizontal_pass!(tree, elements, derivatives, theta)
    downward_pass!(tree, elements)
end

function fmm!(elements, derivatives, expansion_order, n_per_branch, theta=4; reset_elements=true)
    if reset_elements; reset_elements!(elements); end
    tree = Tree(elements; expansion_order, n_per_branch)
    fmm!(tree, elements, derivatives, theta)
    return tree
end
