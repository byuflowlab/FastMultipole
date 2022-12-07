struct Options{TF,TI,TV<:Vector{TI}}
    expansion_order::TI
    n_per_branch::TI
    theta::TF
    targets_index::TV
    sources_index::TV
end

function Options(expansion_order, n_per_branch, theta, targets_index::Int, sources_index::Int)
    targets_index = [targets_index]
    sources_index = [sources_index]
    return Options(expansion_order, n_per_branch, theta, targets_index, sources_index)
end

Options(expansion_order, n_per_branch, theta) = Options(expansion_order, n_per_branch, theta, 1, 1)
