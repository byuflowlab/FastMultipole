"""
    Tree(elements; expansion_order=2, n_per_branch=1)

Constructs an octree of the provided element objects.

# Inputs

- `elements`- a struct containing the following members:

    * `bodies::Array{Float64,2}`- 3+4+mxN array containing element positions, strengths, and m other values that must be sorted into the octree
    * `potential::Array{Float64,2}`- 4+12+36xN array of the potential, Jacobian, and Hessian that are reset every iteration (and don't require sorting)
    * `velocity::Array{Float64,2}`- 3xN array of the velocity vectors at each element, reset every iteration, and calculated in post-processing
    * `direct!::Function`- function calculates the direct influence of the body at the specified location
    * `B2M!::Function`- function converts the body's influence into a multipole expansion
"""
function Tree(systems::Tuple, expansion_order, n_per_branch)
    # initialize objects
    buffer_list = Tuple(get_buffer(system) for system in systems)
    index_list = Tuple(collect(1:length(system)) for system in systems)
    inverse_index_list = Tuple(similar(index) for index in index_list)
    buffer_index_list = Tuple(similar(index) for index in index_list)
    n_systems = length(systems)
    TF = eltype(systems[1])
    branches = Vector{MultiBranch{TF,n_systems}}(undef,1)

    # recursively build branches
    i_start = SVector{n_systems,Int32}([Int32(1) for _ in systems])
    i_end = SVector{n_systems,Int32}([Int32(length(system)) for system in systems])
    i_branch = 1
    center, radius = center_radius(systems; scale_radius = 1.00001)
    level = 0
    multi_branch!(branches, systems, buffer_list, index_list, buffer_index_list, i_start, i_end, i_branch, center, radius, level, expansion_order, n_per_branch)

    # invert index
    for (i_system,index) in enumerate(index_list)
        inverse_index = inverse_index_list[i_system]
        for i_body in eachindex(index)
            inverse_index[index[i_body]] = i_body
        end
    end

    # count leaves
    n_leaves = 0
    for branch in branches
        branch.n_branches == 0 && (n_leaves += 1)
    end

    # create leaf index and leaf count
    cumulative_count = Vector{Int64}(undef,n_leaves)
    cumulative_count[1] = 0
    leaf_index = zeros(Int,n_leaves)
    i_leaf_index = 1
    for (i_branch, branch) in enumerate(branches)
        if branch.n_branches == 0
            leaf_index[i_leaf_index] = Int32(i_branch)
            i_leaf_index < n_leaves && (cumulative_count[i_leaf_index+1] = cumulative_count[i_leaf_index] + sum(branch.n_bodies))
            i_leaf_index += 1
        end
    end

    # assemble tree
    tree = MultiTree(branches, Int16(expansion_order), Int32(n_per_branch), index_list, inverse_index_list, leaf_index, cumulative_count)

    return tree
end

function Tree(system, expansion_order, n_per_branch)
    # initialize objects
    buffer = get_buffer(system)
    index = collect(1:length(system))
    inverse_index = similar(index)
    buffer_index = similar(index)
    TF = eltype(system[1,POSITION])
    branches = Vector{SingleBranch{TF}}(undef,1)

    # recursively build branches
    i_start::Int32 = 1
    i_end::Int32 = length(system)
    i_branch = 1
    center, radius = center_radius(system; scale_radius = 1.00001)
    level = 0
    single_branch!(branches, system, buffer, index, buffer_index, i_start, i_end, i_branch, center, radius, level, expansion_order, n_per_branch)

    # invert index
    for i_body in eachindex(index)
        inverse_index[index[i_body]] = i_body
    end

    # count leaves
    n_leaves = 0
    for branch in branches
        branch.n_branches == 0 && (n_leaves += 1)
    end

    # create leaf index and leaf count
    cumulative_count = Vector{Int64}(undef,n_leaves)
    cumulative_count[1] = 0
    leaf_index = zeros(Int,n_leaves)
    i_leaf_index = 1
    for (i_branch, branch) in enumerate(branches)
        if branch.n_branches == 0
            leaf_index[i_leaf_index] = i_branch
            i_leaf_index < n_leaves && (cumulative_count[i_leaf_index+1] = cumulative_count[i_leaf_index] + sum(branch.n_bodies))
            i_leaf_index += 1
        end
    end

    # assemble tree
    tree = SingleTree(branches, Int16(expansion_order), Int32(n_per_branch), index, inverse_index, leaf_index, cumulative_count)

    return tree
end

# Base.eltype(tree::Tree{TF,<:Any}) where TF = TF
Base.eltype(tree::SingleTree{TF}) where TF = TF

function multi_branch!(branches, systems, buffer_list, index_list, buffer_index_list, i_start, i_end, i_branch, center, radius, level, expansion_order, n_per_branch)
    n_branches = Int8(0)
    n_bodies = i_end - i_start .+ Int32(1)
    n_systems = length(systems)
    multipole_expansion = initialize_expansion(expansion_order)
    local_expansion = initialize_expansion(expansion_order)
    if prod(n_bodies .<= n_per_branch) # checks for all element structs => branch is a leaf; no new branches needed
        i_child = Int32(-1)
        branch = eltype(branches)(n_branches, n_bodies, i_child, i_start, center, radius, multipole_expansion, local_expansion, ReentrantLock(), ReentrantLock())
        branches[i_branch] = branch
        return nothing
    else # not a leaf; branch children
        # count elements in each octant
        octant_attendance = zeros(Int32, 8, n_systems)
        for i_system in 1:n_systems # loop over each type
            for i_body in i_start[i_system]:i_end[i_system] # loop over each body
                x = systems[i_system][i_body, POSITION]
                i_octant = get_octant(x, center)
                octant_attendance[i_octant, i_system] += Int32(1)
            end
        end

        # determine number of children
        for i_octant in 1:8
            n_branches += true in (octant_attendance[i_octant,:] .> 0) # check for elements of any type
        end

        # create branch
        i_child = Int32(length(branches) + 1)
        branches[i_branch] = eltype(branches)(n_branches, n_bodies, i_child, i_start, center, radius, multipole_expansion, local_expansion, ReentrantLock(), ReentrantLock())

        # sort bodies
        ## write offsets
        offsets = zeros(Int32,8,n_systems)
        offsets[1,:] .= i_start
        for i=2:8
            offsets[i,:] .= offsets[i-1,:] + octant_attendance[i-1,:]
        end

        ## create counter
        counter = deepcopy(offsets)

        ## sort element indices into the buffer
        for i_system in 1:n_systems
            buffer_loop!(buffer_list[i_system], buffer_index_list[i_system], counter, systems[i_system], index_list[i_system], i_start[i_system], i_end[i_system], i_system, center)
            place_buffer!(systems[i_system], index_list[i_system], buffer_list[i_system], buffer_index_list[i_system], i_start[i_system], i_end[i_system])
            
            # # update index_list
            # index_list[i_system] .= buffer_index_list[i_system]
        end

        # recursively build new branches
        prev_branches = length(branches)
        resize!(branches, prev_branches + n_branches)

        i_tape = 1
        for i_octant in Int8(0):Int8(7) # check each octant for members
            if true in (octant_attendance[i_octant+Int8(1),:] .> 0) # if octant is populated
                child_i_start = offsets[i_octant+Int8(1),:] # index of first member for all element types
                child_i_end = child_i_start + octant_attendance[i_octant+Int8(1),:] .- Int32(1) # if attendence is 0, the end index will be less than the start
                child_radius = radius / 2#(1 << (level + 1))
                child_center = MVector{3}(center)
                for d in Int8(0):Int8(2)
                    child_center[d + Int8(1)] += child_radius * (((i_octant & Int8(1) << d) >> d) * Int8(2) - Int8(1))
                end
                multi_branch!(branches, systems, buffer_list, index_list, buffer_index_list, SVector{n_systems}(child_i_start), SVector{n_systems}(child_i_end), i_tape + prev_branches, SVector{3}(child_center), child_radius, level + 1, expansion_order, n_per_branch)
                i_tape += 1
            end
        end
    end
end

function single_branch!(branches, system, buffer, index, buffer_index, i_start, i_end, i_branch, center, radius, level, expansion_order, n_per_branch)
    n_branches::Int8 = 0
    n_bodies::Int32 = i_end - i_start + 1
    multipole_expansion = initialize_expansion(expansion_order)
    local_expansion = initialize_expansion(expansion_order)
    if n_bodies <= n_per_branch # checks for all element structs => branch is a leaf; no new branches needed
        i_child::Int32 = -1
        branch = eltype(branches)(n_branches, n_bodies, i_child, i_start, center, radius, multipole_expansion, local_expansion, ReentrantLock(), ReentrantLock())
        branches[i_branch] = branch
        return nothing
    else # not a leaf; branch children
        # count elements in each octant
        octant_attendance = zeros(Int32, 8)
        for i_body in i_start:i_end # loop over each body
            x = system[i_body, POSITION]
            i_octant = get_octant(x, center)
            octant_attendance[i_octant] += 1
        end

        # determine number of children
        for i_octant in 1:8
            n_branches += true == (octant_attendance[i_octant] > 0) # check for elements of any type
        end

        # create branch
        i_child = length(branches) + 1
        # i_child::Int32 = length(branches) + 1
        branches[i_branch] = eltype(branches)(n_branches, n_bodies, i_child, i_start, center, radius, multipole_expansion, local_expansion, ReentrantLock(), ReentrantLock())

        # sort bodies
        ## write offsets
        offsets = zeros(Int32,8)
        offsets[1] = i_start
        for i=2:8
            offsets[i] = offsets[i-1] + octant_attendance[i-1]
        end

        ## create counter
        counter = deepcopy(offsets)

        ## sort element indices into the buffer
        buffer_loop!(buffer, buffer_index, counter, system, index, i_start, i_end, center)
        place_buffer!(system, index, buffer, buffer_index, i_start, i_end)

        # recursively build new branches
        prev_branches = length(branches)
        resize!(branches, prev_branches + n_branches)

        i_tape = 1
        for i_octant in 0:7 # check each octant for members
            if octant_attendance[i_octant+Int8(1)] > 0 # if octant is populated
                child_i_start = offsets[i_octant+Int8(1)] # index of first member for all element types
                child_i_end = child_i_start + octant_attendance[i_octant+Int8(1)] - Int32(1) # if attendence is 0, the end index will be less than the start
                child_radius = radius / 2#(1 << (level + 1))
                child_center = MVector{3}(center)
                for d in 0:2
                    child_center[d + Int8(1)] += child_radius * (((i_octant & Int8(1) << d) >> d) * Int8(2) - Int8(1))
                end
                single_branch!(branches, system, buffer, index, buffer_index, child_i_start, child_i_end, i_tape + prev_branches, SVector{3}(child_center), child_radius, level + 1, expansion_order, n_per_branch)
                i_tape += 1
            end
        end
    end
end

Base.eltype(branch::SingleBranch{TF}) where TF = TF
Base.eltype(branch::MultiBranch{TF,<:Any}) where TF = TF

function buffer_loop!(buffer, buffer_index, counter, system::SortWrapper, index, i_start, i_end, center)
    for i_body in i_start:i_end
        x = system[i_body,POSITION]
        i_octant = get_octant(x, center)
        buffer_index[counter[i_octant]] = index[i_body]
        counter[i_octant] += 1
    end
end

function buffer_loop!(buffer, buffer_index, counter, system, index, i_start, i_end, center)
    for i_body in i_start:i_end
        x = system[i_body,POSITION]
        i_octant = get_octant(x, center)
        buffer[counter[i_octant]] = system[i_body]
        buffer_index[counter[i_octant]] = index[i_body]
        counter[i_octant] += 1
    end
end

function buffer_loop!(buffer, buffer_index, counter, system::SortWrapper, index, i_start, i_end, i_system, center)
    for i_body in i_start:i_end
        x = system[i_body,POSITION]
        i_octant = get_octant(x, center)
        buffer_index[counter[i_octant,i_system]] = index[i_body]
        counter[i_octant,i_system] += 1
    end
end

function buffer_loop!(buffer, buffer_index, counter, system, index, i_start, i_end, i_system, center)
    for i_body in i_start:i_end
        x = system[i_body,POSITION]
        i_octant = get_octant(x, center)
        buffer[counter[i_octant,i_system]] = system[i_body]
        buffer_index[counter[i_octant,i_system]] = index[i_body]
        counter[i_octant,i_system] += 1
    end
end

function place_buffer!(system::SortWrapper, index, buffer, buffer_index, i_start, i_end)
    system.index[i_start:i_end] .= view(buffer_index,i_start:i_end)
    index[i_start:i_end] .= view(buffer_index,i_start:i_end)
    return nothing # no need to sort bodies
end

function place_buffer!(system, index, buffer, buffer_index, i_start, i_end)
    # place sorted bodies
    for i in i_start:i_end
        system[i] = buffer[i]
    end
    index[i_start:i_end] .= view(buffer_index,i_start:i_end)
end

function get_buffer(system::SortWrapper)
    return nothing
end

function get_buffer(system)
    return deepcopy(system)
end

function get_index(system::SortWrapper)
    return system.index
end

function get_index(system)
    return collect(1:length(system))
end

"""
Undoes the sort operation performed by the tree.
"""
function unsort!(systems::Tuple, tree::MultiTree)
    for (system, inverse_index) in zip(systems, tree.inverse_index_list)
        unsort!(system, inverse_index)
    end
end

function unsort!(system, tree::SingleTree)
    unsort!(system, tree.inverse_index)
end

function unsort!(system, inverse_index)
    buffer = deepcopy(system)
    for i in 1:length(system)
        buffer[i] = system[inverse_index[i]]
    end
    for i in 1:length(system)
        system[i] = buffer[i]
    end
end

function unsort!(system::SortWrapper, inverse_index)
    system.index .= system.index[inverse_index]
end

"""
Performs the same sort operation as the tree. (Undoes `unsort!` operation.)
"""
function resort!(system, tree::SingleTree)
    buffer = deepcopy(system)
    index = tree.index
    for i in 1:length(system)
        buffer[i] = system[index[i]]
    end
    for i in 1:length(system)
        system[i] = buffer[i]
    end
end

# function resort!(systems::Tuple, tree::Tree)
#     for (i_system, system) in enumerate(systems)
#         buffer = deepcopy(system)
#         index = tree.index_list[i_system]
#         for i in 1:length(system)
#             buffer[i] = system[index[i]]
#         end
#         for i in 1:length(system)
#             system[i] = buffer[i]
#         end
#     end
# end

function center_radius(systems::Tuple; scale_radius = 1.00001)
    x_min, y_min, z_min = systems[1][1,POSITION]
    x_max, y_max, z_max = systems[1][1,POSITION]
    for system in systems
        for i in 1:length(system)
            x, y, z = system[i,POSITION]
            if x < x_min
                x_min = x
            elseif x > x_max
                x_max = x
            end
            if y < y_min
                y_min = y
            elseif y > y_max
                y_max = y
            end
            if z < z_min
                z_min = z
            elseif z > z_max
                z_max = z
            end
        end
    end
    center = SVector{3}((x_max+x_min)/2, (y_max+y_min)/2, (z_max+z_min)/2)
    # TODO: add element smoothing radius here? Note this method creates cubic cells
    radius = max(x_max-center[1], y_max-center[2], z_max-center[3]) * scale_radius # get half of the longest side length of the rectangle
    return SVector{3}(center), radius
end

function center_radius(system; scale_radius = 1.00001)
    x_min, y_min, z_min = system[1,POSITION]
    x_max, y_max, z_max = system[1,POSITION]
    for i in 1:length(system)
        x, y, z = system[i,POSITION]
        if x < x_min
            x_min = x
        elseif x > x_max
            x_max = x
        end
        if y < y_min
            y_min = y
        elseif y > y_max
            y_max = y
        end
        if z < z_min
            z_min = z
        elseif z > z_max
            z_max = z
        end
    end
    center = SVector{3}((x_max+x_min)/2, (y_max+y_min)/2, (z_max+z_min)/2)
    # TODO: add element smoothing radius here? Note this method creates cubic cells
    radius = max(x_max-center[1], y_max-center[2], z_max-center[3]) * scale_radius # get half of the longest side length of the rectangle
    return SVector{3}(center), radius
end

@inline function get_octant(x, center)
    return (UInt8((x[1] > center[1])) + UInt8(x[2] > center[2]) << 0b1 + UInt8(x[3] > center[3]) << 0b10) + 0b1
end

function n_terms(expansion_order, dimensions)
    n = 0
    for order in 0:expansion_order
        # leverage stars and bars theorem
        n += binomial(order + dimensions - 1, dimensions - 1)
    end
    return n
end

function initialize_expansion(expansion_order)
    return Tuple(zeros(Complex{Float64}, ((expansion_order+1) * (expansion_order+2)) >> 1) for _ in 1:4)
end

function reset_expansions!(tree)
    Threads.@threads for branch in tree.branches
        for dim in 1:4
            branch.multipole_expansion[dim] .*= 0
            branch.local_expansion[dim] .*= 0
        end
    end
end
