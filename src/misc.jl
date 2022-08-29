function vectorial(n)
    prod(factorial.(n))
end

# function ^(a::Vector{TF},b::Vector{TF}) where TF
#     prod([a[i]^b[i] for i in 1:length(a)])
# end

# function ^(a::Vector{TF1},b::Vector{TF2}) where {TF1,TF2}
#     prod([a[i]^b[i] for i in 1:length(a)])
# end

# function ^(a::Vector{TF},b::TF) where TF
#     prod([a[i]^b for i in 1:length(a)])
# end

# function ^(a::Vector{TF1},b::TF2) where {TF1,TF2}
#     prod([a[i]^b for i in 1:length(a)])
# end

function ^(a::V, i, j, k) where {V<:AbstractVector}
    a[1]^i * a[2]^j * a[3]^k
end

##
## Iterators
##
struct VectorIndex{dims,TI}
    i_startndex::NTuple{dims,TI}
    max_sum::TI
end

function iterate(vi::VectorIndex)
    return collect(vi.i_startndex), (collect(vi.i_startndex), 1)
end

function iterate(vi::VectorIndex, state)
    index, i = state
    next = increment(index, i, vi.i_startndex, vi.max_sum)
    if isnothing(next)
        return nothing
    else
        return next[1], next
    end
end

function increment(index, i, i_startndex, max_sum) where {TI,dims}
    s = sum(index)
    if s < max_sum
        index[i] += 1
        return index, i
    elseif index[i] > i_startndex[i] # max reached and ith value can be decremented
        index[i] = i_startndex[i]
        index[i+1] += 1
        return index, i
    elseif i < length(index)-1 # ith value cannot be decremented and end of index not yet reached
        i += 1
        return increment(index, i, i_startndex, max_sum)
    else # end of index reached
        return nothing
    end
end

# function increment!(index, i, index0, max_sum)
#     s = sum(index)
#     if s < max_sum
#         index[i] += 1
#         return index
#     elseif index[i] > index0[i] # max reached and ith value can be decremented
#         index[i] = index0[i]
#         index[i+1] += 1
#         return index
#     elseif i < length(index)-1 # ith value cannot be decremented and end of index not yet reached
#         i += 1
#         return increment!(index, i, index0, max_sum)
#     else # end of index reached
#         return nothing
#     end
# end

function ijk_2_index(i,j,k)
    order = i + j + k
    dimensions = 3
    index = n_terms(order-1,dimensions)
    for ii in order:-1:0
        for jj in order-ii:-1:0
            index += 1
            kk = order - ii - jj
            if ii == i && jj == j && kk == k
                return index
            end
        end
    end
    # need an explicit formula
end
