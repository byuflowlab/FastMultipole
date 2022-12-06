struct Options{TF,NT}
    expansion_order::Int64
    n_per_branch::Int64
    theta::TF
    targets_index::SVector{NT,Int8}
end

function Options(expansion_order, n_per_branch, theta, targets_index::Int)
    targets_index = SVector{1}(Int8(targets_index))
    return Options(expansion_order, n_per_branch, theta, targets_index)
end
