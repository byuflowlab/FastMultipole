
function P2M!(i_branch, tree, elements, ::Cartesian)
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

function M2M!(i_branch, tree, ::Cartesian)
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

function M2L!(i_local, j_multipole, tree, elements, derivatives, ::Cartesian)
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

function L2L!(j_source, tree, ::Cartesian)
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
function L2P!(i_branch, tree, elements, ::Cartesian)
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
