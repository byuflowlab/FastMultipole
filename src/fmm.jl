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

# function l2p(p, x_target, x_local, x_multipole, x_source, q_source, kernel::Kernel{TF,dims}; verbose=false, store_series=nothing) where {TF,dims}
#     Rho = x_target - x_local
#     value = 0.0
#     for k in VectorIndex(ntuple(_ -> 0, dims),p)
#         if verbose
#             val = 1/vectorial(k) * Rho^k
#             # @show k val
#             if !isnothing(store_series)
#                 push!(store_series[1], val)
#             end
#         end
#         value += 1/vectorial(k) * Rho^k * m2l(p, k, x_local, x_multipole, x_source, q_source, kernel; verbose, store_series)
#     end
#     return value
# end

# function m2m(x_M, x_mu, n, dims, x_source, q_source)
#     x_M_mu = x_M - x_mu
#     val = 0.0
#     for k in VectorIndex(ntuple(_ -> 0, dims), n)
#         val += 1 / vectorial(n-k) * x_M_mu ^ (n-k) * p2m(n, x_mu, x_source, q_source)
#     end
#     return val
# end

# function m2m(elements::AbstractArray{e}, n, mu::l1, M::l2, kernel::Kernel{TF,dims}) where {TF,dims,e<:Element, l1<:Tree, l2<:Tree}
#     x_mu = mu.center
#     x_M = M.center
#     x_M_mu = x_M - x_mu
#     val = 0.0
#     for k in VectorIndex(ntuple(_ -> 0, dims), n)
#         val += 1/vectorial(n-k) * x_M_mu ^ (n-k) * p2m(n, x_mu, x_source, q_source)
#     end
#     return val
# end

# function m2l(p, k, x_local, x_multipole, x_source, q_source, kernel::Kernel{TF,dims}; verbose=false, store_series=nothing) where {TF,dims}
#     sum_k = sum(k)
#     Rho = x_local - x_multipole
#     value = 0.0
#     one_vec = ones(UInt8,dims)
#     for n in VectorIndex(ntuple(_ -> 0, dims),p-sum_k)
#         if verbose
#             val = kernel(n + k + one_vec, x_multipole, 1.0, x_local ) * p2m(n, x_multipole, x_source, q_source)
#             # @show n val
#             if !isnothing(store_series)
#                 push!(store_series[2], val)
#             end
#         end
#         value += kernel(n + k + one_vec, x_multipole, 1.0, x_local ) * p2m(n, x_multipole, x_source, q_source)
#     end
#     return value
# end

# # how to make this package?
# # upward pass: add up p2m and m2m from the bottom up (zero first)
# # horizontal pass: m2l to the root/branch/leaf where the potential is desired
# # say we desire the potential at all points in a tree
# # then, m2l at the root level only, from each branch to each branch
# # then, l2l from each branch downward to each child
# # finally, l2p from each leaf to its points

# function p2m(n, x_multipole, x_source, q_source)
#     x_ms = x_multipole - x_source
#     return q_source / vectorial(n) * x_ms ^ n
# end

# # formulate for multiple sources
# function p2m(elements::AbstractArray{e}, leaf::Leaf, n) where e<:Element
#     val = 0.0
#     x_multipole = leaf.center
#     for i in leaf.i_start:leaf.i_end
#         val += (x_multipole - get_X(elements,i))^n * get_q(elements,i)
#     end
#     return val / vectorial(n)
# end

# # function p2m(elements::AbstractArray{e}, n, branch::Branch) where e<:Element
# #     val = 0.0
# #     for vestige in branch.branches
# #         val += p2m(elements, n, vestige)
# #     end
# #     return val
# # end

# # function p2m(elements::AbstractArray{e}, n, root::Root) where e<:Element
# #     val = 0.0
# #     for vestige in root.branches
# #         val += p2m(elements, n, vestige)
# #     end
# #     return val
# # end

# function m2l(elements::AbstractArray{e}, leaf::Leaf, p, k, x_local, kernel::Kernel{TF,dims}; verbose=false, store_series=nothing) where {TF,dims,e<:Element}
#     sum_k = sum(k)
#     x_multipole = leaf.center
#     Rho = x_local - x_multipole
#     value = 0.0
#     one_vec = ones(UInt8,dims)
#     for n in VectorIndex(ntuple(_ -> 0, dims),p-sum_k)
#         if verbose
#             val = kernel(n + k + one_vec, x_multipole, 1.0, x_local ) * p2m(elements, leaf, n)
#             # @show n val
#             if !isnothing(store_series)
#                 push!(store_series[2], val)
#             end
#         end
#         value += kernel(n + k + one_vec, x_multipole, 1.0, x_local ) * p2m(elements, leaf, n)
#     end
#     return value
# end

# function l2p(elements::AbstractArray{e}, leaf::Leaf, p, x_target, x_local, kernel::Kernel{TF,dims}; verbose=false, store_series=nothing) where {TF,dims,e<:Element}
#     Rho = x_target - x_local
#     value = 0.0
#     for k in VectorIndex(ntuple(_ -> 0, dims),p)
#         if verbose
#             val = 1/vectorial(k) * Rho^k
#             # @show k val
#             if !isnothing(store_series)
#                 push!(store_series[1], val)
#             end
#         end
#         value += 1/vectorial(k) * Rho^k * m2l(elements, leaf, p, k, x_local, kernel; verbose, store_series)
#     end
#     return value
# end

# # formulate for multiple targets
# function l2p(source_elements::AbstractArray{e1}, target_elements::AbstractArray{e2}, source_leaf::Leaf, target_leaf::Leaf, p, kernel::Kernel{TF,dims}; verbose=false, store_series=nothing) where {TF,dims,e1<:Element,e2<:Element}
#     x_local = target_leaf.center
#     n_targets = length(target_leaf)
#     Phis = zeros(n_targets)

#     for i_target in 1:n_targets
#         Rho = get_X(target_elements,target_leaf.i_start + i_target - 1) - x_local
#         value = 0.0

#         for k in VectorIndex(ntuple(_ -> 0, dims),p)
#             if verbose
#                 val = 1/vectorial(k) * Rho^k
#                 # @show k val
#                 if !isnothing(store_series)
#                     push!(store_series[1], val)
#                 end
#             end
#             value += 1/vectorial(k) * Rho^k * m2l(source_elements, source_leaf, p, k, x_local, kernel; verbose, store_series)
#         end

#         Phis[i_target] = value
#     end

#     return Phis
# end

# # putting it all together
# function l2p(elements::AbstractArray{e}, leaf::Leaf, dims, p) where e<:Element
#     vals = Vector{eltype(elements)}(undef,length(leaf))
#     for i in 1:length(leaf)
#         val = 0.0
#         for k in VectorIndex(ntuple(_ -> 0, dims),p)
#             val += 1/vectorial(k) * (get_X(elements,leaf.i_start + i - 1) - leaf.center)^k * l2l(elements, leaf, k, p)
#         end
#         vals[i] = val
#     end
# end

# function l2l(elements::AbstractArray{e}, leaf::Leaf, n, p) where e<:Element
#     val = 0.0
#     for k in VectorIndex(n,p)
#         val += 1/vectorial(k-n) * (leaf.center - leaf.parent.center)^(k-n) * m2l(elements, leaf.parent, k)
#     end
# end

# function m2l(elements::AbstractArray{e}, source_root::Root, target_root::Root, kernel::Kernel{TF,dims}) where {TF,dims,e<:Element}
#     reset_local!(target_root)

#     for target_vestige in target_root.branches
#         for source_vestige in source_root.branches
#             if !(source_vestige === target_vestige)
#                 for n in VectorIndex(ntuple(_ -> 0, dims),p-k)
#                 target_vestige.local_data[1] +=0
#                 end
#             end
#         end
#     end
# end

# function m2m(elements::AbstractArray{e}, root::Root, n) where e<:Element
#     val = 0.0
#     for vestige in root.branches
#         val += m2m(elements, vestige, n)
#     end
# end

# function m2m(elements::AbstractArray{e}, branch::Branch, n) where e<:Element
#     val = 0.0
#     x_M = branch.center
#     for vestige in branch.branches
#         x_mu = vestige.center
#         for k in VectorIndex(ntuple(_ -> 0, dims), n)
#             val += 1 / vectorial(n - k) * (x_M - x_mu)^(n-k) * m2m(elements, vestige, n)
#         end
#     end
#     return val
# end

# function m2m(elements::AbstractArray{e}, leaf::Leaf, n) where e<:Element
#     val = 0.0
#     x_M = leaf.center
#     for k in VectorIndex(ntuple(_ -> 0, dims), n)
#         val += 1 / vectorial(n - k) * p2m(elements, leaf, n)
#     end
#     return val
# end

# function reset_multipole!(root::Root)
#     for branch in root.branches
#         reset_multipole!(branch)
#     end
# end

# function reset_multipole!(branch::Branch)
#     branch.multipole_data[1] *= 0
#     for branch in branch.branches
#         reset_multipole!(branch)
#     end
# end

# function reset_multipole!(leaf::Leaf)
#     leaf.multipole_data[1] *= 0
#     return nothing
# end

# function reset_local!(root::Root)
#     for branch in root.branches
#         reset_local!(branch)
#     end
# end

# function reset_local!(branch::Branch)
#     branch.local_data[1] *= 0
#     for branch in branch.branches
#         reset_local!(branch)
#     end
# end

# function reset_local!(leaf::Leaf)
#     leaf.local_data[1] *= 0
#     return nothing
# end
