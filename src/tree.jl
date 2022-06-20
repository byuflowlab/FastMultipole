abstract type Tree{TF,TI} end

struct Leaf{TF,TI} <: Tree{TF,TI}
    center::Vector{TF}
    i_start::TI
    i_end::TI
end

struct Branch{TF,TI} <: Tree{TF,TI}
    center::Vector{TF}
    i_start::TI
    i_end::TI
    branches::Vector{Union{Branch{TF,TI}, Leaf{TF,TI}}}
end

struct Root{TF,TI} <: Tree{TF,TI}
    center::Vector{TF}
    branches::Vector{Union{Branch{TF,TI}, Leaf{TF,TI}}}
    cardinality::TI
    n_divide::TI
end

function Root(elements::AbstractArray{e}; n_divide=2^get_dims(elements), cardinality=10) where e<:Element
    rectangle = get_rectangle(elements) # encompasses all elements
    center = get_center(rectangle) # root center
    if length(elements) <= cardinality # no partitioning necessary
        vestiges = Vector{Union{Leaf{eltype(elements), typeof(n_divide)},Branch{eltype(elements), typeof(n_divide)}}}(undef,1)
        vestiges[1] = Leaf(center, 1, length(elements))
    else
        vestiges = branch!(elements, 1, length(elements), center, cardinality, n_divide)
    end
    return Root(center, vestiges, cardinality, n_divide)
end

function branch!(elements::AbstractArray{e}, i_start, i_end, center, cardinality, n_divide) where e<:Element

    # divide the data
    i_partitions, centers = divide!(elements, i_start, i_end, center, cardinality, n_divide)

    # preallocate vestiges
    vestiges = Vector{Union{Branch{eltype(elements), typeof(n_divide)}, Leaf{eltype(elements), typeof(n_divide)}}}(undef,2^n_divide)

    for i_partition in 1:2^n_divide  # loop over each partition
        # check cardinality for the current partition
        if i_partitions[2,i_partition] - i_partitions[1,i_partition] < cardinality
            vestiges[i_partition] = Leaf(centers[:,i_partition], i_partitions[1,i_partition], i_partitions[2,i_partition])
        else # if not satisfied, recursively branch! the partition
            vestiges[i_partition] = Branch(centers[:,i_partition], i_partitions[1,i_partition], i_partitions[2,i_partition],
                branch!(elements, i_partitions[1,i_partition], i_partitions[2,i_partition], centers[:,i_partition], cardinality, n_divide))
        end
    end

    return vestiges
end

"Sorts and partitions a vector of elements."
function divide_history!(elements, center, cardinality, n_divide)
    n_elements = length(elements)
    if n_elements < cardinality # cardinality satisfied
        return [Leaf(center, 1, n_elements)]
    else # cardinality not satisfied; partition needed
        dims = get_dims(elements)
        i_partitions = Array{Int64,2}(undef,2,2^(n_divide+1)-1) # i_partitions[:,1] spans the entire (i_start, i_end)
        centers = Array{eltype(elements),2}(undef,dims,2^(n_divide+1)-1) # i_partitions[:,1] spans the entire (i_start, i_end)
        centers[:,1] .= center[:]
        i_partitions[1,1] = 1
        i_partitions[2,1] = n_elements
        for i_divide in 1:n_divide

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
function divide!(elements, i_start, i_end, center, cardinality, n_divide)
    dims = get_dims(elements)
    i_partitions = Array{Int64,2}(undef,2,2^(n_divide+1)-1) # i_partitions[:,1] spans the entire (i_start, i_end)
    centers = Array{eltype(elements),2}(undef,dims,2^(n_divide+1)-1) # i_partitions[:,1] spans the entire (i_start, i_end)
    centers[:,1] .= center[:]
    i_partitions[1,1] = i_start
    i_partitions[2,1] = i_end
    for i_divide in 1:n_divide
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
    return i_partitions[:,2^n_divide:end], centers[:,2^n_divide:end]
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
