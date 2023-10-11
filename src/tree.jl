const BRANCH_TYPE = Float64
global SHRINKING_OFFSET = .000001

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
function Tree(systems, options::Options, ::Val{TF}=Val{Float64}()) where TF
    # unpack options
    expansion_order = options.expansion_order
    n_per_branch = options.n_per_branch

    # initialize objects
    buffer_list = Tuple(get_buffer(system) for system in systems)
    index_list = Tuple(collect(1:length(system)) for system in systems)
    inverse_index_list = Tuple(similar(index) for index in index_list)
    buffer_index_list = Tuple(similar(index) for index in index_list)
    n_systems = length(systems)
    branches = Vector{Branch{TF,n_systems}}(undef,1)

    # recursively build branches
    i_start = SVector{n_systems,Int32}([Int32(1) for _ in systems])
    i_end = SVector{n_systems,Int32}([Int32(length(system)) for system in systems])
    i_branch = 1
    center, radius = center_radius(systems; scale_radius = 1.00001)
    level = 0
    branch!(branches, systems, buffer_list, index_list, buffer_index_list, i_start, i_end, i_branch, center, radius, level, expansion_order, n_per_branch)

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
    tree = Tree(branches, Int16(expansion_order), Int32(n_per_branch), index_list, inverse_index_list, leaf_index, cumulative_count)

    if options.shrinking
        update_radius(systems,branches,1; options.second_pass)
    end

    return tree
end

function branch!(branches, systems, buffer_list, index_list, buffer_index_list, i_start, i_end, i_branch, center, radius, level, expansion_order, n_per_branch)
    n_branches = Int8(0)
    n_bodies = i_end - i_start .+ Int32(1)
    n_systems = length(systems)
    multipole_expansion = initialize_expansion(expansion_order)
    local_expansion = initialize_expansion(expansion_order)
    if prod(n_bodies .<= n_per_branch) # checks for all element structs => branch is a step_through_bodies; no new branches needed
        i_child = Int32(-1)
        branch = eltype(branches)(n_branches, n_bodies, i_child, i_start, center, radius, multipole_expansion, local_expansion, ReentrantLock(), ReentrantLock())
        branches[i_branch] = branch
        return nothing
    else # not a step_through_bodies; branch children
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
                branch!(branches, systems, buffer_list, index_list, buffer_index_list, SVector{n_systems}(child_i_start), SVector{n_systems}(child_i_end), i_tape + prev_branches, SVector{3}(child_center), child_radius, level + 1, expansion_order, n_per_branch)
                i_tape += 1
            end
        end
    end
end

#TODO: update radius and center
function update_radius(systems::Tuple, branches, branch_index; second_pass=true)
    branch = branches[branch_index]
    if branch.n_branches != 0
        for b = 0:branch.n_branches - 1
            update_radius(systems, branches, branch.first_branch + b; second_pass)
        end
        step_through_branches(branches, branch_index; second_pass)
    else
        step_through_bodies(systems, branches, branch_index; second_pass)
    end
    return nothing
end

function step_through_branches(branches, branch_index; second_pass=true)
    branch = branches[branch_index]

    first_child_index = branch.first_branch
    first_child = branches[first_child_index]

    rectangle = MVector{6, Float64}(
        first_child.center[1]+first_child.radius, first_child.center[1]-first_child.radius,
        first_child.center[2]+first_child.radius, first_child.center[2]-first_child.radius,
        first_child.center[3]+first_child.radius, first_child.center[3]-first_child.radius
        )

    for index = 0:branch.n_branches - 1
        child = branches[first_child_index + index]
        rectangle[1] = max(child.center[1]+child.radius, rectangle[1])
        rectangle[3] = max(child.center[2]+child.radius, rectangle[3])
        rectangle[5] = max(child.center[3]+child.radius, rectangle[5])
        rectangle[2] = min(child.center[1]-child.radius, rectangle[2])
        rectangle[4] = min(child.center[2]-child.radius, rectangle[4])
        rectangle[6] = min(child.center[3]-child.radius, rectangle[6])
    end

    new_center = typeof(branch.center)(
        (rectangle[1] + rectangle[2]) / 2,
        (rectangle[3] + rectangle[4]) / 2,
        (rectangle[5] + rectangle[6]) / 2
    )


    # if branch.n_branches == 1
    #     branch.center[1] = (rectangle[1] + rectangle[2]) / 2
    #     branch.center[2] = (rectangle[3] + rectangle[4]) / 2
    #     branch.center[3] = (rectangle[5] + rectangle[6]) / 2
    # else
    #     branch.center[1] = (rectangle[1] + rectangle[2]) / 2
    #     branch.center[2] = (rectangle[3] + rectangle[4]) / 2
    #     branch.center[3] = (rectangle[5] + rectangle[6]) / 2
    # end

    new_radius = zero(branch.radius)
    if second_pass
        for index = 0:branch.n_branches - 1
            child = branches[first_child_index + index]
            distance = sqrt((child.center[1] - new_center[1])^2 + (child.center[2] - new_center[2])^2 + (child.center[3] - new_center[3])^2) + child.radius
            new_radius = max(new_radius, distance)
        end
    else 
        new_radius = sqrt(
            (rectangle[1] - new_center[1])^2 + 
            (rectangle[3] - new_center[2])^2 + 
            (rectangle[5] - new_center[3])^2
            )
    end

    (; n_branches, n_bodies, first_branch, first_body, center, multipole_expansion, local_expansion, lock, child_lock) = branch
    new_branch = Branch(n_branches, n_bodies, first_branch, first_body, new_center, new_radius, multipole_expansion, local_expansion, lock, child_lock)
    branches[branch_index] = new_branch
    return nothing
end

function step_through_bodies(systems::Tuple, branches, branch_index; second_pass=true)
    branch = branches[branch_index]
    
    first_index = branch.first_body[1]
    rectangle = MVector{6, Float64}(
        systems[1][first_index, POSITION][1], systems[1][first_index, POSITION][1],
        systems[1][first_index, POSITION][2], systems[1][first_index, POSITION][2],
        systems[1][first_index, POSITION][3], systems[1][first_index, POSITION][3]
        )
    for (i_system, system) in enumerate(systems)
        for i_body in branch.first_body[i_system]:(branch.first_body[i_system]+branch.n_bodies[i_system]-1)
            body = system[i_body,POSITION]
            body_r = system[i_body,RADIUS]
            rectangle[1] = max(body[1]+body_r, rectangle[1])
            rectangle[3] = max(body[2]+body_r, rectangle[3])
            rectangle[5] = max(body[3]+body_r, rectangle[5])
            rectangle[2] = min(body[1]-body_r, rectangle[2])
            rectangle[4] = min(body[2]-body_r, rectangle[4])
            rectangle[6] = min(body[3]-body_r, rectangle[6])

        end
    end

    new_center = (branch.n_bodies == 1) ? typeof(branch.center)(
        (rectangle[1] + rectangle[2]) / 2 + SHRINKING_OFFSET,
        (rectangle[3] + rectangle[4]) / 2 + SHRINKING_OFFSET,
        (rectangle[5] + rectangle[6]) / 2 + SHRINKING_OFFSET
        ) : typeof(branch.center)(
            (rectangle[1] + rectangle[2]) / 2,
            (rectangle[3] + rectangle[4]) / 2,
            (rectangle[5] + rectangle[6]) / 2
            )

    (; n_branches, n_bodies, first_branch, first_body, center, radius, multipole_expansion, local_expansion, lock, child_lock) = branch

    if second_pass
        new_radius = zero(branch.radius)
        for (i_system, system) in enumerate(systems)
            for i_body in branch.first_body[i_system]:(branch.first_body[i_system]+branch.n_bodies[i_system]-1)
                body = system[i_body,POSITION]
                body_r = system[i_body,RADIUS]
                distance = sqrt((body[1] - new_center[1])^2 + (body[2]- new_center[2])^2 + (body[3] - new_center[3])^2) + body_r
                new_radius = max(new_radius, distance)
            end
        end
    else
        (; n_branches, n_bodies, first_branch, first_body, center, multipole_expansion, local_expansion, lock, child_lock) = branch
        new_radius = sqrt((rectangle[1] - new_center[1])^2 + (rectangle[3] - new_center[2])^2 + (rectangle[5] - new_center[3])^2)
    end

    new_branch = Branch(n_branches, n_bodies, first_branch, first_body, new_center, new_radius, multipole_expansion, local_expansion, lock, child_lock)
    branches[branch_index] = new_branch
    return nothing

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
function unsort!(systems::Tuple, tree::Tree)
    for (i_system, system) in enumerate(systems)
        unsort!(system, tree.inverse_index_list[i_system])
    end
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
function resort!(systems::Tuple, tree::Tree)
    for (i_system, system) in enumerate(systems)
        buffer = deepcopy(system)
        index = tree.index_list[i_system]
        for i in 1:length(system)
            buffer[i] = system[index[i]]
        end
        for i in 1:length(system)
            system[i] = buffer[i]
        end
    end
end

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
