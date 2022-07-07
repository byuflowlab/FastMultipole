function l2p(p, x_target, x_local, x_multipole, x_source, q_source, kernel::Kernel{TF,dims}; verbose=false, store_series=nothing) where {TF,dims}
    Rho = x_target - x_local
    value = 0.0
    for k in VectorIndex(ntuple(_ -> 0, dims),p)
        if verbose
            val = 1/vectorial(k) * Rho^k
            # @show k val
            if !isnothing(store_series)
                push!(store_series[1], val)
            end
        end
        value += 1/vectorial(k) * Rho^k * m2l(p, k, x_local, x_multipole, x_source, q_source, kernel; verbose, store_series)
    end
    return value
end

function m2m(x_M, x_mu, n, dims, x_source, q_source)
    x_M_mu = x_M - x_mu
    val = 0.0
    for k in VectorIndex(ntuple(_ -> 0, dims), n)
        val += 1 / vectorial(n-k) * x_M_mu ^ (n-k) * p2m(n, x_mu, x_source, q_source)
    end
    return val
end

function m2m(elements::AbstractArray{e}, n, mu::l1, M::l2, kernel::Kernel{TF,dims}) where {TF,dims,e<:Element, l1<:Tree, l2<:Tree}
    x_mu = mu.center
    x_M = M.center
    x_M_mu = x_M - x_mu
    val = 0.0
    for k in VectorIndex(ntuple(_ -> 0, dims), n)
        val += 1/vectorial(n-k) * x_M_mu ^ (n-k) * p2m(n, x_mu, x_source, q_source)
    end
    return val
end

function m2l(p, k, x_local, x_multipole, x_source, q_source, kernel::Kernel{TF,dims}; verbose=false, store_series=nothing) where {TF,dims}
    sum_k = sum(k)
    Rho = x_local - x_multipole
    value = 0.0
    one_vec = ones(UInt8,dims)
    for n in VectorIndex(ntuple(_ -> 0, dims),p-sum_k)
        if verbose
            val = kernel(n + k + one_vec, x_multipole, 1.0, x_local ) * p2m(n, x_multipole, x_source, q_source)
            # @show n val
            if !isnothing(store_series)
                push!(store_series[2], val)
            end
        end
        value += kernel(n + k + one_vec, x_multipole, 1.0, x_local ) * p2m(n, x_multipole, x_source, q_source)
    end
    return value
end

# how to make this package?
# upward pass: add up p2m and m2m from the bottom up (zero first)
# horizontal pass: m2l to the root/branch/leaf where the potential is desired
# say we desire the potential at all points in a tree
# then, m2l at the root level only, from each branch to each branch
# then, l2l from each branch downward to each child
# finally, l2p from each leaf to its points

function p2m(n, x_multipole, x_source, q_source)
    x_ms = x_multipole - x_source
    return q_source / vectorial(n) * x_ms ^ n
end

# formulate for multiple sources
function p2m(elements::AbstractArray{e}, leaf::Leaf, n) where e<:Element
    val = 0.0
    x_multipole = leaf.center
    for i in leaf.i_start:leaf.i_end
        val += (x_multipole - get_X(elements,i))^n * get_q(elements,i)
    end
    return val / vectorial(n)
end

# function p2m(elements::AbstractArray{e}, n, branch::Branch) where e<:Element
#     val = 0.0
#     for vestige in branch.branches
#         val += p2m(elements, n, vestige)
#     end
#     return val
# end

# function p2m(elements::AbstractArray{e}, n, root::Root) where e<:Element
#     val = 0.0
#     for vestige in root.branches
#         val += p2m(elements, n, vestige)
#     end
#     return val
# end

function m2l(elements::AbstractArray{e}, leaf::Leaf, p, k, x_local, kernel::Kernel{TF,dims}; verbose=false, store_series=nothing) where {TF,dims,e<:Element}
    sum_k = sum(k)
    x_multipole = leaf.center
    Rho = x_local - x_multipole
    value = 0.0
    one_vec = ones(UInt8,dims)
    for n in VectorIndex(ntuple(_ -> 0, dims),p-sum_k)
        if verbose
            val = kernel(n + k + one_vec, x_multipole, 1.0, x_local ) * p2m(elements, leaf, n)
            # @show n val
            if !isnothing(store_series)
                push!(store_series[2], val)
            end
        end
        value += kernel(n + k + one_vec, x_multipole, 1.0, x_local ) * p2m(elements, leaf, n)
    end
    return value
end

function l2p(elements::AbstractArray{e}, leaf::Leaf, p, x_target, x_local, kernel::Kernel{TF,dims}; verbose=false, store_series=nothing) where {TF,dims,e<:Element}
    Rho = x_target - x_local
    value = 0.0
    for k in VectorIndex(ntuple(_ -> 0, dims),p)
        if verbose
            val = 1/vectorial(k) * Rho^k
            # @show k val
            if !isnothing(store_series)
                push!(store_series[1], val)
            end
        end
        value += 1/vectorial(k) * Rho^k * m2l(elements, leaf, p, k, x_local, kernel; verbose, store_series)
    end
    return value
end

# formulate for multiple targets
function l2p(source_elements::AbstractArray{e1}, target_elements::AbstractArray{e2}, source_leaf::Leaf, target_leaf::Leaf, p, kernel::Kernel{TF,dims}; verbose=false, store_series=nothing) where {TF,dims,e1<:Element,e2<:Element}
    x_local = target_leaf.center
    n_targets = length(target_leaf)
    Phis = zeros(n_targets)

    for i_target in 1:n_targets
        Rho = get_X(target_elements,target_leaf.i_start + i_target - 1) - x_local
        value = 0.0

        for k in VectorIndex(ntuple(_ -> 0, dims),p)
            if verbose
                val = 1/vectorial(k) * Rho^k
                # @show k val
                if !isnothing(store_series)
                    push!(store_series[1], val)
                end
            end
            value += 1/vectorial(k) * Rho^k * m2l(source_elements, source_leaf, p, k, x_local, kernel; verbose, store_series)
        end

        Phis[i_target] = value
    end

    return Phis
end

# putting it all together
function l2p(elements::AbstractArray{e}, leaf::Leaf, dims, p) where e<:Element
    vals = Vector{eltype(elements)}(undef,length(leaf))
    for i in 1:length(leaf)
        val = 0.0
        for k in VectorIndex(ntuple(_ -> 0, dims),p)
            val += 1/vectorial(k) * (get_X(elements,leaf.i_start + i - 1) - leaf.center)^k * l2l(elements, leaf, k, p)
        end
        vals[i] = val
    end
end

function l2l(elements::AbstractArray{e}, leaf::Leaf, n, p) where e<:Element
    val = 0.0
    for k in VectorIndex(n,p)
        val += 1/vectorial(k-n) * (leaf.center - leaf.parent.center)^(k-n) * m2l(elements, leaf.parent, k)
    end
end

function m2l(elements::AbstractArray{e}, source_root::Root, target_root::Root, kernel::Kernel{TF,dims}) where {TF,dims,e<:Element}
    reset_local!(target_root)

    for target_vestige in target_root.branches
        for source_vestige in source_root.branches
            if !(source_vestige === target_vestige)
                for n in VectorIndex(ntuple(_ -> 0, dims),p-k)
                target_vestige.local_data[1] +=0
                end
            end
        end
    end
end

function m2m(elements::AbstractArray{e}, root::Root, n) where e<:Element
    val = 0.0
    for vestige in root.branches
        val += m2m(elements, vestige, n)
    end
end

function m2m(elements::AbstractArray{e}, branch::Branch, n) where e<:Element
    val = 0.0
    x_M = branch.center
    for vestige in branch.branches
        x_mu = vestige.center
        for k in VectorIndex(ntuple(_ -> 0, dims), n)
            val += 1 / vectorial(n - k) * (x_M - x_mu)^(n-k) * m2m(elements, vestige, n)
        end
    end
    return val
end

function m2m(elements::AbstractArray{e}, leaf::Leaf, n) where e<:Element
    val = 0.0
    x_M = leaf.center
    for k in VectorIndex(ntuple(_ -> 0, dims), n)
        val += 1 / vectorial(n - k) * p2m(elements, leaf, n)
    end
    return val
end

function reset_multipole!(root::Root)
    for branch in root.branches
        reset_multipole!(branch)
    end
end

function reset_multipole!(branch::Branch)
    branch.multipole_data[1] *= 0
    for branch in branch.branches
        reset_multipole!(branch)
    end
end

function reset_multipole!(leaf::Leaf)
    leaf.multipole_data[1] *= 0
    return nothing
end

function reset_local!(root::Root)
    for branch in root.branches
        reset_local!(branch)
    end
end

function reset_local!(branch::Branch)
    branch.local_data[1] *= 0
    for branch in branch.branches
        reset_local!(branch)
    end
end

function reset_local!(leaf::Leaf)
    leaf.local_data[1] *= 0
    return nothing
end
