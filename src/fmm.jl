"Sum charges at each level."
function upward_pass!(root::Root, elements::AbstractArray{e}) where e<:Element
    # reset coefficients
    reset_multipole!(root)

    # iterate over branches/leaves
    for (i,level) in enumerate(root.branches)
        for branch in level
            X_branch = get_X(branch)
            nursery = i > 1 ? root.branches[i-1][branch.children] : elements[branch.children]
            for child in nursery
                q = get_q(child)
                set_q(branch.element, get_q(branch.element) + q)
                dx_complex = complexify(get_X(child) - X_branch)
                for i_coefficient in 1:root.p_expansion
                    if i == 1
                        branch.multipole_coefficients[i_coefficient] -= q * dx_complex^i_coefficient / i_coefficient
                    else
                        branch.multipole_coefficients[i_coefficient] += -q / i_coefficient * dx_complex^i_coefficient + sum([child.multipole_coefficients[k] * dx_complex^(i_coefficient-k) * binomial(i_coefficient-1,k-1) for k in 1:i_coefficient])
                    end
                end
            end
        end
    end
end

function complexify(vec)
    @assert length(vec) == 2 "Cannot complexify vector of length $(length(vec))"
    return vec[1] + vec[2] * im
end


"returns the nearest neighbor list of the specified object"
function i_nearest_neighbors(root, level, ci)
    this_size = size(root.branches[level])[1]
    nn = CartesianIndices(Tuple([(ci[i] > 1 ? ci[i] - 1 : ci[i]):(ci[i] < this_size ? ci[i] + 1 : ci[i]) for i in 1:length(ci)]))
    # nearest_branches = root.branches[level][nn]
    # return nearest_branches
end

"Interaction list consists of all children of the nearest neighbors of the specified cell's parents which are not the cell's nearest neighbors."
function interaction_list(root, level, ci)
    dims = length(ci)
    if level >= length(root.branches) - 1 # empty interaction list at top 2 levels
        return CartesianIndex{dims}[]
    end

    # nearest neighbors of the cell's parent (including the parent)
    i_parent_nn = i_nearest_neighbors(root, level+1, root.branches[level][ci].parent[1])

    # children of the parents nearest neighbors
    candidates = vcat([root.branches[level+1][nn].children for nn in i_parent_nn]...)

    # which ones are well separated
    i_child_nn = i_nearest_neighbors(root, level, ci)
    i_well_separated = candidates[findall(x -> !in(x, i_child_nn), candidates)]

    return i_well_separated
end

"Approximate influence of all cells on all other cells"
function downward_pass!(root, elements, kernel)
    # reset local expansion coefficients
    reset_local!(root)

    # reset element potential
    reset_potential!(elements)

    for i_level in length(root.branches):-1:1
        # interactions caused by siblings in the interaction list
        for target_ci in CartesianIndices(root.branches[i_level])
            il = interaction_list(root, i_level, target_ci)
            for source_ci in il
                if abs(get_q(root.branches[i_level][source_ci])) > 0.0 # note that the local expansion is nonzero if there are 0 elements and must therefore be omitted
                    # get z_0
                    dx_complex = complexify(get_X(root.branches[i_level][source_ci]) - get_X(root.branches[i_level][target_ci]))
                    q_source = get_q(root.branches[i_level][source_ci])

                    # if source_ci == CartesianIndex(5,3) && target_ci == CartesianIndex(3,5)
                    #     println("SHERLOCK!\n\tBEFORE: $(root.branches[i_level][target_ci].local_coefficients)")
                    # end
                    # 0th coefficient of the power series
                    root.branches[i_level][target_ci].local_coefficients[1] += q_source * log(-dx_complex) +
                        sum([root.branches[i_level][source_ci].multipole_coefficients[k] / (-dx_complex)^k for k in 1:root.p_expansion])

                    # remainding coefficients
                    for i_coefficient in 1:root.p_expansion
                        root.branches[i_level][target_ci].local_coefficients[i_coefficient+1] += dx_complex^(-i_coefficient) * (
                            -q_source/i_coefficient + sum([root.branches[i_level][source_ci].multipole_coefficients[k] * (-dx_complex)^(-k) * binomial(i_coefficient+k-1, k-1) for k in 1:root.p_expansion])
                        )
                    end
                    # if source_ci == CartesianIndex(5,3) && target_ci == CartesianIndex(3,5)
                    #     println("SHERLOCK!\n\tAFTER: $(root.branches[i_level][target_ci].local_coefficients)")
                    # end
                end
            end
        end
        # translate local expansions to children; note that child local expansion coefficients should be zero at this point
        if i_level > 1 # translate local expansion to children
            for parent_ci in CartesianIndices(root.branches[i_level]) # select parent
                branch = root.branches[i_level][parent_ci]
                # for child_ci in CartesianIndices(root.branches[i_level-1][branch.children])
                for child_ci in branch.children
                    child = root.branches[i_level-1][child_ci]
                    dx_complex = complexify(get_X(branch) - get_X(child))
                    # if parent_ci == CartesianIndex(2,3) && child_ci == CartesianIndex(3,5)
                    #     println("Sherlock!\n\tBEFORE: $(child.local_coefficients)")
                    # end
                    child.local_coefficients .+= branch.local_coefficients
                    for j in 0:root.p_expansion-1
                        for k in root.p_expansion-j:root.p_expansion
                            child.local_coefficients[k] -= dx_complex * child.local_coefficients[k+1]
                        end
                    end
                    # if parent_ci == CartesianIndex(2,3) && child_ci == CartesianIndex(3,5)
                    #     println("\tAFTER: $(child.local_coefficients)")
                    # end
                end
            end
        end
    end
end

function evaluate!(root, elements, kernel)
    # evaluate local expansion at each particle location
    for leaf in root.branches[1]
        for element in elements[leaf.children]
            dx_complex = complexify(get_X(element) - get_X(leaf))
            add_V(element, real(sum([leaf.local_coefficients[l+1] * dx_complex^l for l in 0:root.p_expansion])))
        end
    end

    # direct calculation of nearest neighbors
    for target_ci in CartesianIndices(root.branches[1])
        target_i = root.branches[1][target_ci].children
        for source_ci in i_nearest_neighbors(root, 1, target_ci)
            source_i = root.branches[1][source_ci].children
            direct!(elements[source_i], elements[target_i], kernel)
        end
    end
end

function fmm!(root, elements, kernel)
    upward_pass!(root, elements)
    downward_pass!(root, elements, kernel)
    evaluate!(root, elements, kernel)
    return nothing
end

function reset_local!(root::Root)
    for level in root.branches
        for branch in level
            branch.local_coefficients[:] *= 0
        end
    end
end

function reset_multipole!(root::Root)
    for level in root.branches
        for branch in level
            branch.multipole_coefficients[:] *= 0
        end
    end
end

function reset_potential!(elements::AbstractArray{e}) where e<:Element
    for element in elements
        set_V(element, 0.0)
    end
end

#####
##### taylor series cartesian fmm
#####
function cfmm_upward_pass!(root::CRoot, elements::AbstractArray{e}) where e<:Element
    # reset coefficients
    reset_multipole!(root)

    # iterate over branches/leaves
    for (i,level) in enumerate(root.branches)
        for branch in level
            X_branch = get_X(branch)
            nursery = i > 1 ? root.branches[i-1][branch.children] : elements[branch.children]
            for child in nursery
                q = get_q(child)
                set_q(branch.element, get_q(branch.element) + q)
                dx_complex = complexify(get_X(child) - X_branch)
                for i_coefficient in 1:root.p_expansion
                    if i == 1
                        branch.multipole_coefficients[i_coefficient] -= q * dx_complex^i_coefficient / i_coefficient
                    else
                        branch.multipole_coefficients[i_coefficient] += -q / i_coefficient * dx_complex^i_coefficient + sum([child.multipole_coefficients[k] * dx_complex^(i_coefficient-k) * binomial(i_coefficient-1,k-1) for k in 1:i_coefficient])
                    end
                end
            end
        end
    end
end
