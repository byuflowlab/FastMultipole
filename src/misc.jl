function vectorial(n)
    prod(factorial.(n))
end

function ^(a::V, i, j, k) where {V<:AbstractVector}
    a[1]^i * a[2]^j * a[3]^k
end


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
