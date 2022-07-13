"Sum charges at each level."
function upward_pass!(root::Root, elements::AbstractArray{e}) where e<:Element
    # iterate over branches/leaves
    for (i,level) in enumerate(root.branches)
        for branch in level
            X_branch = get_X(branch)
            branch.coefficients .*= 0.0 # reset coefficients
            nursery = i > 1 ? root.branches[i-1][branch.children] : elements[branch.children]
            for child in nursery
                q = get_q(child)
                dx_complex = complexify(get_X(child) - X_branch)
                set_q(branch.element, get_q(branch.element) + q)
                for i_coefficient in 1:root.p_expansion
                    branch.coefficients[i_coefficient] -= q * dx_complex^i_coefficient / i_coefficient
                end
            end
        end
    end
    # THIS IS STILL MISSING A WAY OF TRANSLATING MULTIPOLE COEFFICIENTS TO PARENTS
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
    for i_level in length(root.branches):-1:1
        # sibling interactions
        for target_ci in CartesianIndices(root.branches[i_level])
            il = interaction_list(root, i_level, target_ci)
            for source_ci in il
                direct!(root.branches[i_level][source_ci].element, root.branches[i_level][target_ci].element, kernel)
            end
        end

        # translate to children
        nursery = i_level > 1 ? root.branches[i_level-1] : elements
        for branch in root.branches[i_level] # iterate over each branch after calculating sibling influence
            for child_ci in branch.children # iterate over all children
                add_V(nursery[child_ci], get_V(branch))
            end
        end
    end

    # direct calculation at leaf level
    for target_ci in CartesianIndices(root.branches[1])
        target_i = root.branches[1][target_ci].children
        for source_ci in i_nearest_neighbors(root, 1, target_ci)
            source_i = root.branches[1][source_ci].children
            direct!(elements[source_i], elements[target_i], kernel)
        end
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
