struct Options{TF,TI,TB<:Bool}
    expansion_order::TI
    n_per_branch::TI
    theta::TF
    shrinking::TB
    second_pass::TB
end

function Options(expansion_order, n_per_branch, theta)
    return Options(expansion_order,n_per_branch,theta,false,false)
end



