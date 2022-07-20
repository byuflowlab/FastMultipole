abstract type Tree{dims} end

struct Branch{dims} <: Tree{dims}
    element
    multipole_coefficients
    local_coefficients
    N
    LS
    SS
    parent
    children
end

function get_X(branch::Branch)
    get_X(branch.element)
end

function get_q(branch::Branch)
    get_q(branch.element)
end

function get_V(branch::Branch)
    get_V(branch.element)
end

function add_V(branch::Branch, summand)
    add_V(branch.element, summand)
end

function set_X(branch::Branch, new_X)
    set_X(branch.element, new_X)
end

function Branch(element_type, p_expansion, dims)
    Branch{dims}(element_type(dims), zeros(Complex{Float64},p_expansion), zeros(Complex{Float64},p_expansion+1), [0], zeros(dims), zeros(dims), CartesianIndex{dims}[], CartesianIndex{dims}[])
end

function Leaf(element_type, p_expansion, dims)
    Branch{dims}(element_type(dims), zeros(Complex{Float64},p_expansion), zeros(Complex{Float64},p_expansion+1), [0], zeros(dims), zeros(dims), CartesianIndex{dims}[], Int32[])
end

struct Root{dims} <: Tree{dims}
    branches::Vector{Array{Branch,dims}}
    branch_limit
    elements_per_leaf
    p_expansion
end

function Root(elements::AbstractArray{e}, branch_limit, elements_per_leaf, p_expansion; n_divides=nothing, rect=nothing) where e<:Element
    # determine grid refinement
    dims = get_dims(elements)
    if n_divides == nothing
        n_cells = length(elements) / elements_per_leaf
        n_divides = Int(floor(log2(n_cells)))
        if isodd(n_divides); n_divides += 1; end
        n_divides /= dims # number of levels in the tree to reach the desired elements per leaf (roughly)
    end

    # bottom level cell discretization
    if isnothing(rect)
        rect = get_rectangle(elements)
    end
    lengths = [rect[2,i] - rect[1,i] for i in 1:dims]
    size = Int(2^n_divides)
    fences = [range(rect[1,i], stop=rect[2,i], length=1+size) for i in 1:dims]
    # fences = [range(rect[1,i] - lengths[i]/100, stop=rect[2,i] + lengths[i]/100, length=1+size) for i in 1:dims]
    # initialize leaves
    element_type = typeof(elements[1])
    leaves = Array{Branch,dims}(undef,fill(size,dims)...)
    for ci in CartesianIndices(leaves)
        leaves[ci] = Leaf(element_type, p_expansion, dims)
        set_X(leaves[ci], [(fences[i][ci[i]] + fences[i][ci[i]+1])/2 for i in 1:dims])
    end

    # assign leaf membership
    for (i_element,element) in enumerate(elements)
        X = get_X(element)
        cartesian_i = CartesianIndex([findfirst(x -> x>X[i], fences[i]) - 1 for i in 1:dims]...)
        push!(leaves[cartesian_i].children, i_element)
        leaves[cartesian_i].N[1] += 1
    end

    # build branch system
    branches = [leaves]
    merge_branches!(branches, p_expansion, dims, element_type)

    return Root(branches, branch_limit, elements_per_leaf, p_expansion)
end

function merge_branches!(branches, p_expansion, dims, element_type)
    highest_level = branches[end]
    top_size = size(highest_level)[1]
    if top_size == 1 # note: if only 1 leaf exists, no branches are generated
                     #     if more than 1 leaf exists, branches are generated
        return branches
    else
        # merge branches/leaves
        new_size = Int(top_size / 2)
        next_level = Array{Branch,dims}(undef,fill(new_size,dims)...)
        for ci in CartesianIndices(Tuple(fill(new_size,dims)))
            children = CartesianIndices(Tuple([2*ci[j]-1:2*ci[j] for j in 1:dims]))
            center = S.mean([get_X(branch) for branch in highest_level[children]])
            element = element_type(dims)
            set_X(element, center)
            multipole_coefficients = zeros(Complex{eltype(center)}, p_expansion)
            local_coefficients = zeros(Complex{eltype(center)}, p_expansion+1)
            N = [sum([highest_level[j].N[1] for j in children])]
            LS = 0.0
            SS = 0.0
            next_level[ci] = Branch{dims}(element, multipole_coefficients, local_coefficients, N, LS, SS, CartesianIndex{dims}[], reshape(children,length(children)))

            # update parent info
            for child in highest_level[children]
                push!(child.parent, ci)
            end
        end
        push!(branches, next_level)

        merge_branches!(branches, p_expansion, dims, element_type)
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

function linear_2_cartesian(linear_i, size, dims)
    cartesian_index = Base._ind2sub(fill(size, dims), linear_i)
    return cartesian_index
end

function get_rectangle(elements::AbstractArray{e}; rel_offset = 1e-2) where e<:Element
    dims = get_dims(elements)
    rectangle = zeros(2, dims)
    for dim in 1:dims
        xs = [get_X(element)[dim] for element in elements]
        rectangle[1,dim] = minimum(xs)*(1-rel_offset)
        rectangle[2,dim] = maximum(xs)*(1+rel_offset)
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
