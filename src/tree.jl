struct Branch
    center
    N
    LS
    SS
    indices
end

function Branch(dims)
    Branch(zeros(dims), [0], zeros(dims), zeros(dims), Int32[])
end

struct Root
    branch_system::Vector{Vector{Branch}}
    branch_limit
    elements_per_leaf
end

function Root(elements::AbstractArray{e}, branch_limit, elements_per_leaf) where e<:Element
    # determine grid refinement
    n_cells = length(elements) / elements_per_leaf
    dims = get_dims(elements)
    n_divides = Int(floor(log2(n_cells)))
    if isodd(n_divides); n_divides += 1; end
    n_divides /= dims # number of levels in the tree to reach the desired elements per leaf (roughly)

    # bottom level cell discretization
    rect = get_rectangle(elements)
    lengths = [rect[2,i] - rect[1,i] for i in 1:dims]
    size = Int(2^n_divides)
    fences = [range(rect[1,i] - lengths[i]/100, stop=rect[2,i] + lengths[i]/100, length=1+size) for i in 1:dims]
    leaves = [Branch(dims) for i in 1:size^dims]

    for (i_element,element) in enumerate(elements)
        X = get_X(element)
        cartesian_i = CartesianIndex([findfirst(x -> x>X[i], fences[i]) - 1 for i in 1:dims]...)
        i_leaf = cartesian_2_linear(cartesian_i, size)
        push!(leaves[i_leaf].indices, i_element)
        leaves[i_leaf].N[1] += 1
    end

    # build branch system
    top_size = length(fences[1]) - 1
    branch_system = branch(leaves, top_size, dims)

    return Root(branch_system, branch_limit, elements_per_leaf)
end

function branch(leaves::Vector{Branch}, top_size, dims)
    branch_system = Vector{Vector{Branch}}(undef,1)
    branch_system[1] = leaves
    return branch(branch_system, top_size, dims)
end

function branch(branch_system::Vector{Vector{Branch}}, top_size, dims)
    top_n = length(branch_system[end])
    if top_n == 1
        return branch_system
    else
        # merge branches/leaves
        new_n = Int(top_n / 2^dims)
        new_size = Int(top_size / 2)
        new_branch = Vector{Branch}(undef,new_n)
        for (i,ci) in enumerate(CartesianIndices(Tuple(fill(new_size,dims))))
            indices = [
                cartesian_2_linear((2*ci[1]-1, 2*ci[2]-1), new_size),
                cartesian_2_linear((2*ci[1], 2*ci[2]-1), new_size),
                cartesian_2_linear((2*ci[1]-1, 2*ci[2]), new_size),
                cartesian_2_linear((2*ci[1], 2*ci[2]), new_size)
            ]
            center = S.mean([branch.center for branch in branch_system[end][indices]])
            # center = [fences[j][2*ci[j]] for j in 1:dims]
            N = [sum([branch_system[end][j].N[1] for j in indices])]
            LS = 0.0
            SS = 0.0
            new_branch[i] = Branch(center, N, LS, SS, indices)
        end
        push!(branch_system, new_branch)

        return branch(branch_system, Int(top_size/2), dims)
    end
end

function cartesian_2_linear(cartesian_i, size)
    linear_index = cartesian_i[1]
    mul = size
    for i in 2:length(cartesian_i)
        linear_index += mul * (cartesian_i[i] - 1)
        mul *= size
    end
    return linear_index
end

function get_rectangle(elements::AbstractArray{e}) where e<:Element
    dims = get_dims(elements)
    rectangle = zeros(2, dims)
    for dim in 1:dims
        xs = [get_X(element)[dim] for element in elements]
        rectangle[1,dim] = minimum(xs)
        rectangle[2,dim] = maximum(xs)
    end
    return rectangle
end

function get_center(rectangle)
    center = S.mean(rectangle, dims=1)[:]
end

function get_center(elements::AbstractArray{e}) where e<:Element
    rectangle = get_rectangle(elements)
    return get_center(rectangle)
end







# abstract type Tree{TF,TI} end

# struct Leaf{TF,TI} <: Tree{TF,TI}
#     center::Vector{TF}
#     multipole_data::Vector{TF}
#     local_data::Vector{TF}
#     i_start::TI
#     i_end::TI
# end

# struct Branch{TF,TI} <: Tree{TF,TI}
#     center::Vector{TF}
#     multipole_data::Vector{TF}
#     local_data::Vector{TF}
#     i_start::TI
#     i_end::TI
#     branches::Vector{Union{Branch{TF,TI}, Leaf{TF,TI}}}
# end

# struct Root{TF,TI} <: Tree{TF,TI}
#     center::Vector{TF}
#     branches::Vector{Union{Branch{TF,TI}, Leaf{TF,TI}}}
#     cardinality::TI
#     n_divide::TI
# end

# function Root(elements::AbstractArray{e}; n_divide=2^get_dims(elements), cardinality=10) where e<:Element
#     rectangle = get_rectangle(elements) # encompasses all elements
#     center = get_center(rectangle) # root center
#     if length(elements) <= cardinality # no partitioning necessary
#         vestiges = Vector{Union{Leaf{eltype(elements), typeof(n_divide)},Branch{eltype(elements), typeof(n_divide)}}}(undef,1)
#         vestiges[1] = Leaf(center, zeros(eltype(elements),1), zeros(eltype(elements),1), 1, length(elements))
#     else
#         vestiges = branch!(elements, 1, length(elements), center, cardinality, n_divide)
#     end
#     return Root(center, vestiges, cardinality, n_divide)
# end

# function branch!(elements::AbstractArray{e}, i_start, i_end, center, cardinality, n_divide) where e<:Element

#     # divide the data
#     i_partitions, centers = divide!(elements, i_start, i_end, center, cardinality, n_divide)

#     # preallocate vestiges
#     vestiges = Vector{Union{Branch{eltype(elements), typeof(n_divide)}, Leaf{eltype(elements), typeof(n_divide)}}}(undef,2^n_divide)

#     for i_partition in 1:2^n_divide  # loop over each partition
#         # check cardinality for the current partition
#         if i_partitions[2,i_partition] - i_partitions[1,i_partition] < cardinality
#             vestiges[i_partition] = Leaf(centers[:,i_partition], zeros(eltype(elements),1), zeros(eltype(elements),1), i_partitions[1,i_partition], i_partitions[2,i_partition])
#         else # if not satisfied, recursively branch! the partition
#             vestiges[i_partition] = Branch(centers[:,i_partition], zeros(eltype(elements),1), zeros(eltype(elements),1), i_partitions[1,i_partition], i_partitions[2,i_partition],
#                 branch!(elements, i_partitions[1,i_partition], i_partitions[2,i_partition], centers[:,i_partition], cardinality, n_divide))
#         end
#     end

#     return vestiges
# end

# "Sorts and partitions a vector of elements."
# function divide_history!(elements, center, cardinality, n_divide)
#     n_elements = length(elements)
#     if n_elements < cardinality # cardinality satisfied
#         return [Leaf(center, 1, n_elements)]
#     else # cardinality not satisfied; partition needed
#         dims = get_dims(elements)
#         i_partitions = Array{Int64,2}(undef,2,2^(n_divide+1)-1) # i_partitions[:,1] spans the entire (i_start, i_end)
#         centers = Array{eltype(elements),2}(undef,dims,2^(n_divide+1)-1) # i_partitions[:,1] spans the entire (i_start, i_end)
#         centers[:,1] .= center[:]
#         i_partitions[1,1] = 1
#         i_partitions[2,1] = n_elements
#         for i_divide in 1:n_divide

#             for i_partition in 1:2^(i_divide-1) # split each previous partition in 2; this can be parallelized

#                 # determine which dimension to partition
#                 dim = (i_divide-1) % dims + 1 # this is a cheap, simple, hopefully effective method; perhaps plug in other functions here

#                 # recall the i_start, i_end of the partition we are dividing
#                 i_dividend = 2^(i_divide-1) + i_partition - 1 # index of the centers and partitions being divided
#                 i_start, i_end = i_partitions[:,i_dividend]

#                 # sort along the chosen dimension
#                 sort!(view(elements, i_start:i_end); by=x->get_X(x)[dim], alg=QuickSort)
#                 # sort!(elements; by=x->get_X(x)[dim], alg=PartialQuickSort(i_start:i_end))

#                 # get partition index
#                 x_p = centers[dim,i_dividend]
#                 partition_index = findfirst(element -> get_X(element)[dim] > x_p, view(elements, i_start:i_end)) + i_start - 1

#                 # populate new partition
#                 i_new = 2^i_divide + 2 * (i_partition-1)
#                 i_partitions[1,i_new] = i_start
#                 i_partitions[2,i_new] = partition_index - 1
#                 i_partitions[1,i_new+1] = partition_index
#                 i_partitions[2,i_new+1] = i_end

#                 # get centers
#                 centers[:,i_new] .= get_center(view(elements, i_start:partition_index-1))[:]
#                 centers[:,i_new+1] .= get_center(view(elements, partition_index:i_end))[:]
#             end
#         end
#         return i_partitions, centers
#     end
# end

# "Sorts and partitions a vector of elements."
# function divide!(elements, i_start, i_end, center, cardinality, n_divide)
#     dims = get_dims(elements)
#     i_partitions = Array{Int64,2}(undef,2,2^(n_divide+1)-1) # i_partitions[:,1] spans the entire (i_start, i_end)
#     centers = Array{eltype(elements),2}(undef,dims,2^(n_divide+1)-1) # i_partitions[:,1] spans the entire (i_start, i_end)
#     centers[:,1] .= center[:]
#     i_partitions[1,1] = i_start
#     i_partitions[2,1] = i_end
#     for i_divide in 1:n_divide
#         for i_partition in 1:2^(i_divide-1) # split each previous partition in 2; this can be parallelized

#             # determine which dimension to partition
#             dim = (i_divide-1) % dims + 1 # this is a cheap, simple, hopefully effective method; perhaps plug in other functions here

#             # recall the i_start_local, i_end of the partition we are dividing
#             i_dividend = 2^(i_divide-1) + i_partition - 1 # index of the centers and partitions being divided
#             i_start_local, i_end_local = i_partitions[:,i_dividend]

#             # sort along the chosen dimension
#             sort!(view(elements, i_start_local:i_end_local); by=x->get_X(x)[dim], alg=QuickSort)
#             # sort!(elements; by=x->get_X(x)[dim], alg=PartialQuickSort(i_start_local:i_end_local))

#             # get partition index
#             x_p = centers[dim,i_dividend]
#             # x_debug = [get_X(e) for e in view(elements, i_start_local:i_end_local)]
#             # @show x_debug x_p
#             partition_index = findfirst(element -> get_X(element)[dim] > x_p, view(elements, i_start_local:i_end_local)) + i_start_local - 1

#             # populate new partition
#             i_new = 2^i_divide + 2 * (i_partition-1)
#             i_partitions[1,i_new] = i_start_local
#             i_partitions[2,i_new] = partition_index - 1
#             i_partitions[1,i_new+1] = partition_index
#             i_partitions[2,i_new+1] = i_end_local

#             # get centers
#             centers[:,i_new] .= get_center(view(elements, i_start_local:partition_index-1))[:]
#             centers[:,i_new+1] .= get_center(view(elements, partition_index:i_end_local))[:]
#         end
#     end
#     return i_partitions[:,2^n_divide:end], centers[:,2^n_divide:end]
# end

# function length(leaf::Leaf)
#     return leaf.i_end - leaf.i_start + 1
# end
