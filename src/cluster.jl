abstract type Tree{TF,TI} end

struct Leaf{TF,TI} <: Tree{TF,TI}
    center::Vector{TF}
    indices::AbstractArray{TI}
end

struct Branch{TF,TI} <: Tree{TF,TI}
    center::Vector{TF}
    branches::Vector{branch} where branch <: Tree{TF,TI}
    indices::Vector{TI}
end

struct Root{TF,TI} <: Tree{TF,TI}
    vestiges::Vector{Vector{branch}} where branch <: Tree{TF,TI}
    max_elements::TI
    max_branches::TI
end

function Root(elements::AbstractArray{e}; max_branches=2^get_dims(elements), max_elements=10) where e<:Element
    branches = Branch(elements, max_elements, max_branches)
    return Root(branches, max_elements, max_branches)
end

function Branch(elements::Vector{e}, max_elements, max_branches) where e <: Element
    rectangle = get_rectangle(elements) # encompasses all elements
    center = get_center(rectangle) # root center
    if i_end - i_start + 1 <= max_elements # form leaf
        return [Leaf(center, fill(1,length(elements)))]
    else
        nothing
    end
end

function branch!(elements::AbstractArray{e}, center, max_elements, max_branches) where e<:Element
    # divide the data
    i_partitions, centers = divide!(elements, i_start, i_end, center, max_elements, max_branches)

    # preallocate vestiges
    vestiges = Vector{Union{Branch{eltype(elements), typeof(max_branches)}, Leaf{eltype(elements), typeof(max_branches)}}}(undef,2^max_branches)

    for i_partition in 1:2^max_branches  # loop over each partition
        # check max_elements for the current partition
        if i_partitions[2,i_partition] - i_partitions[1,i_partition] < max_elements
            vestiges[i_partition] = Leaf(centers[:,i_partition], zeros(eltype(elements),1), zeros(eltype(elements),1), i_partitions[1,i_partition], i_partitions[2,i_partition])
        else # if not satisfied, recursively branch! the partition
            vestiges[i_partition] = Branch(centers[:,i_partition], zeros(eltype(elements),1), zeros(eltype(elements),1), i_partitions[1,i_partition], i_partitions[2,i_partition],
                branch!(elements, i_partitions[1,i_partition], i_partitions[2,i_partition], centers[:,i_partition], max_elements, max_branches))
        end
    end

    return vestiges
end

function Root(elements::AbstractArray{e}; max_branches=2^get_dims(elements), max_elements=10) where e<:Element
    rectangle = get_rectangle(elements) # encompasses all elements
    center = get_center(rectangle) # root center
    if length(elements) <= max_elements # no partitioning necessary
        vestiges = Vector{Union{Leaf{eltype(elements), typeof(max_branches)},Branch{eltype(elements), typeof(max_branches)}}}(undef,1)
        vestiges[1] = Leaf(center, zeros(eltype(elements),1), zeros(eltype(elements),1), 1, length(elements))
    else
        vestiges = branch!(elements, 1, length(elements), center, max_elements, max_branches)
    end
    return Root(center, vestiges, max_elements, max_branches)
end

function branch!(elements::AbstractArray{e}, i_start, i_end, center, max_elements, max_branches) where e<:Element

    # divide the data
    i_partitions, centers = divide!(elements, i_start, i_end, center, max_elements, max_branches)

    # preallocate vestiges
    vestiges = Vector{Union{Branch{eltype(elements), typeof(max_branches)}, Leaf{eltype(elements), typeof(max_branches)}}}(undef,2^max_branches)

    for i_partition in 1:2^max_branches  # loop over each partition
        # check max_elements for the current partition
        if i_partitions[2,i_partition] - i_partitions[1,i_partition] < max_elements
            vestiges[i_partition] = Leaf(centers[:,i_partition], zeros(eltype(elements),1), zeros(eltype(elements),1), i_partitions[1,i_partition], i_partitions[2,i_partition])
        else # if not satisfied, recursively branch! the partition
            vestiges[i_partition] = Branch(centers[:,i_partition], zeros(eltype(elements),1), zeros(eltype(elements),1), i_partitions[1,i_partition], i_partitions[2,i_partition],
                branch!(elements, i_partitions[1,i_partition], i_partitions[2,i_partition], centers[:,i_partition], max_elements, max_branches))
        end
    end

    return vestiges
end

"Sorts and partitions a vector of elements."
function divide_history!(elements, center, max_elements, max_branches)
    n_elements = length(elements)
    if n_elements < max_elements # max_elements satisfied
        return [Leaf(center, 1, n_elements)]
    else # max_elements not satisfied; partition needed
        dims = get_dims(elements)
        i_partitions = Array{Int64,2}(undef,2,2^(max_branches+1)-1) # i_partitions[:,1] spans the entire (i_start, i_end)
        centers = Array{eltype(elements),2}(undef,dims,2^(max_branches+1)-1) # i_partitions[:,1] spans the entire (i_start, i_end)
        centers[:,1] .= center[:]
        i_partitions[1,1] = 1
        i_partitions[2,1] = n_elements
        for i_divide in 1:max_branches

            for i_partition in 1:2^(i_divide-1) # split each previous partition in 2; this can be parallelized

                # determine which dimension to partition
                dim = (i_divide-1) % dims + 1 # this is a cheap, simple, hopefully effective method; perhaps plug in other functions here

                # recall the i_start, i_end of the partition we are dividing
                i_dividend = 2^(i_divide-1) + i_partition - 1 # index of the centers and partitions being divided
                i_start, i_end = i_partitions[:,i_dividend]

                # sort along the chosen dimension
                sort!(view(elements, i_start:i_end); by=x->get_X(x)[dim], alg=QuickSort)
                # sort!(elements; by=x->get_X(x)[dim], alg=PartialQuickSort(i_start:i_end))

                # get partition index
                x_p = centers[dim,i_dividend]
                partition_index = findfirst(element -> get_X(element)[dim] > x_p, view(elements, i_start:i_end)) + i_start - 1

                # populate new partition
                i_new = 2^i_divide + 2 * (i_partition-1)
                i_partitions[1,i_new] = i_start
                i_partitions[2,i_new] = partition_index - 1
                i_partitions[1,i_new+1] = partition_index
                i_partitions[2,i_new+1] = i_end

                # get centers
                centers[:,i_new] .= get_center(view(elements, i_start:partition_index-1))[:]
                centers[:,i_new+1] .= get_center(view(elements, partition_index:i_end))[:]
            end
        end
        return i_partitions, centers
    end
end

"Sorts and partitions a vector of elements."
function divide!(elements, i_start, i_end, center, max_elements, max_branches)
    dims = get_dims(elements)
    i_partitions = Array{Int64,2}(undef,2,2^(max_branches+1)-1) # i_partitions[:,1] spans the entire (i_start, i_end)
    centers = Array{eltype(elements),2}(undef,dims,2^(max_branches+1)-1) # i_partitions[:,1] spans the entire (i_start, i_end)
    centers[:,1] .= center[:]
    i_partitions[1,1] = i_start
    i_partitions[2,1] = i_end
    for i_divide in 1:max_branches
        for i_partition in 1:2^(i_divide-1) # split each previous partition in 2; this can be parallelized

            # determine which dimension to partition
            dim = (i_divide-1) % dims + 1 # this is a cheap, simple, hopefully effective method; perhaps plug in other functions here

            # recall the i_start_local, i_end of the partition we are dividing
            i_dividend = 2^(i_divide-1) + i_partition - 1 # index of the centers and partitions being divided
            i_start_local, i_end_local = i_partitions[:,i_dividend]

            # sort along the chosen dimension
            sort!(view(elements, i_start_local:i_end_local); by=x->get_X(x)[dim], alg=QuickSort)
            # sort!(elements; by=x->get_X(x)[dim], alg=PartialQuickSort(i_start_local:i_end_local))

            # get partition index
            x_p = centers[dim,i_dividend]
            # x_debug = [get_X(e) for e in view(elements, i_start_local:i_end_local)]
            # @show x_debug x_p
            partition_index = findfirst(element -> get_X(element)[dim] > x_p, view(elements, i_start_local:i_end_local)) + i_start_local - 1

            # populate new partition
            i_new = 2^i_divide + 2 * (i_partition-1)
            i_partitions[1,i_new] = i_start_local
            i_partitions[2,i_new] = partition_index - 1
            i_partitions[1,i_new+1] = partition_index
            i_partitions[2,i_new+1] = i_end_local

            # get centers
            centers[:,i_new] .= get_center(view(elements, i_start_local:partition_index-1))[:]
            centers[:,i_new+1] .= get_center(view(elements, partition_index:i_end_local))[:]
        end
    end
    return i_partitions[:,2^max_branches:end], centers[:,2^max_branches:end]
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

function length(leaf::Leaf)
    return leaf.i_end - leaf.i_start + 1
end
