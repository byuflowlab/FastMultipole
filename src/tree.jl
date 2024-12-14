#####
##### tree constructor
#####
function Tree(system; expansion_order=7, leaf_size=100, n_divisions=20, shrink_recenter=false, accumulate_charge=false, allocation_safety_factor=1.0, estimate_cost=false, read_cost_file=false, write_cost_file=false)
    # ensure `system` isn't empty; otherwise return an empty tree
    if get_n_bodies(system) > 0
        # initialize variables
        octant_container = get_octant_container(system) # preallocate octant counter; records the total number of bodies in the first octant to start out, but can be reused to count bodies per octant later on
        cumulative_octant_census = get_octant_container(system) # for creating a cumsum of octant populations
        bodies_index = get_bodies_index(system)

        # root branch
        center, box = center_box(system)
        bx, by, bz = box

        # initial octree generation uses cubic cells
        bmax = max(max(bx,by),bz)
        radius = sqrt(bmax*bmax*3.0)
        box = SVector{3}(bmax, bmax, bmax)

        # prepare to divide
        i_first_branch = 2
        buffer = get_buffer(system)
        sort_index = get_sort_index(system)
        sort_index_buffer = get_sort_index_buffer(system)

        # grow root branch
        root_branch, n_children, i_leaf = Branch(system, sort_index, octant_container, buffer, sort_index_buffer, i_first_branch, bodies_index, center, center, radius, radius, box, box, 0, 1, leaf_size, expansion_order) # even though no sorting needed for creating this branch, it will be needed later on; so `Branch` not ony_min creates the root_branch, but also sorts itself into octants and returns the number of children it will have so we can plan array size
        branches = [root_branch] # this first branch will already have its child branches encoded
        # estimated_n_branches = estimate_n_branches(system, leaf_size, allocation_safety_factor)
        # sizehint!(branches, estimated_n_branches)
        parents_index = 1:1

        # grow branches
        levels_index = [parents_index] # store branches at each level
        for i_divide in 1:n_divisions
            if n_children > 0
                parents_index, n_children, i_leaf = child_branches!(branches, system, sort_index, buffer, sort_index_buffer, i_leaf, leaf_size, parents_index, cumulative_octant_census, octant_container, n_children, expansion_order)
                push!(levels_index, parents_index)
            end
        end

        # check depth
        if n_children > 0
            n_children_prewhile = n_children
            while n_children > 0
                parents_index, n_children, i_leaf = child_branches!(branches, system, sort_index, buffer, sort_index_buffer, i_leaf, leaf_size, parents_index, cumulative_octant_census, octant_container, n_children, expansion_order)
                push!(levels_index, parents_index)
            end
            if WARNING_FLAG_LEAF_SIZE[]
                @warn "leaf_size not reached in for loop, so while loop used to build octree; to improve performance, increase `n_divisions` > $(length(levels_index))"
                WARNING_FLAG_LEAF_SIZE[] = false
            end
        end

        # invert index
        invert_index!(sort_index_buffer, sort_index)
        inverse_sort_index = sort_index_buffer # reuse the buffer as the inverse index

        # accumulate charge (for error prediction and dynamic expansion order)
        if accumulate_charge
            for i_level in length(levels_index):-1:1
                accumulate_charge!(branches, system, levels_index[i_level])
            end
        end

        # shrink and recenter branches to account for bodies of nonzero radius
        shrink_recenter && shrink_recenter!(branches, levels_index, system)

        # store leaves
        leaf_index = Vector{Int}(undef,0)
        sizehint!(leaf_index, length(levels_index[end]))
        for (i_branch,branch) in enumerate(branches)
            if branch.n_branches == 0
                push!(leaf_index, i_branch)
            end
        end

        # # cost parameters
        # if estimate_cost
        #     params, errors, nearfield_params, nearfield_errors = estimate_tau(system; read_cost_file=read_cost_file, write_cost_file=write_cost_file)
        #     cost_parameters = CostParameters(params..., nearfield_params)
        # else
        #     cost_parameters = CostParameters(system)
        # end
        # cost_parameters = Threads.nthreads() > 1 ? direct_cost_estimate(system, leaf_size) : dummy_direct_cost_estimate(system, leaf_size)

        # assemble tree
        tree = Tree(branches, levels_index, leaf_index, sort_index, inverse_sort_index, buffer, Val(expansion_order), leaf_size)#, cost_parameters)

    else
        tree = EmptyTree(system)
    end

    return tree
end

Tree(branches::Vector{<:SingleBranch}, levels_index, leaf_index, sort_index, inverse_sort_index, buffer, expansion_order, leaf_size) =
    SingleTree(branches, levels_index, leaf_index, sort_index, inverse_sort_index, buffer, expansion_order, leaf_size)#, cost_parameters)

Tree(branches::Vector{<:MultiBranch}, levels_index, leaf_index, sort_index, inverse_sort_index, buffer, expansion_order, leaf_size) =
    MultiTree(branches, levels_index, leaf_index, sort_index, inverse_sort_index, buffer, expansion_order, leaf_size)#, cost_parameters)

"""
    EmptyTree(system)

Returns an empty tree. Used if `system` is empty.

# Arguments

* `system`: the system from which a tree is to be created

# Returns

* `tree`: if `typeof(system)<:Tuple`, a `::MultiTree` is returned; otherwise, a `::SingleTree` is returned

"""
function EmptyTree(system)
    branches = Vector{SingleBranch{Float64}}(undef,0)
    levels_index = Vector{UnitRange{Int64}}(undef,0)
    leaf_index = Int[]
    sort_index = Int[]
    inverse_sort_index = Int[]
    buffer = nothing
    expansion_order = Val(-1)
    leaf_size = -1
    return Tree(branches, levels_index, leaf_index, sort_index, inverse_sort_index, buffer, expansion_order, leaf_size)
end

function EmptyTree(system::Tuple)
    N = length(system)
    branches = Vector{MultiBranch{Float64,N}}(undef,0)
    levels_index = Vector{UnitRange{Int64}}(undef,0)
    leaf_index = Int[]
    sort_index_list = Tuple(Int[] for _ in 1:N)
    inverse_sort_index = Tuple(Int[] for _ in 1:N)
    buffer = nothing
    expansion_order = Val(-1)
    leaf_size = -1
    return Tree(branches, levels_index, leaf_index, sort_index_list, inverse_sort_index_list, buffer, expansion_order, leaf_size)
end

@inline function get_n_bodies(systems::Tuple)
    n_bodies = 0
    for system in systems
        n_bodies += get_n_bodies(system)
    end
    return n_bodies
end

function estimate_n_branches(system, leaf_size, allocation_safety_factor)
    n_bodies = get_n_bodies(system)
    estimated_n_divisions = Int(ceil(log(8,n_bodies/leaf_size)))
    estimated_n_branches = div(8^(estimated_n_divisions+1) - 1,7)
    return Int(ceil(estimated_n_branches * allocation_safety_factor))
end

function child_branches!(branches, system, sort_index, buffer, sort_index_buffer, i_leaf, leaf_size, parents_index, cumulative_octant_census, octant_container, n_children, expansion_order)
    i_first_branch = parents_index[end] + n_children + 1
    for i_parent in parents_index
        parent_branch = branches[i_parent]
        if parent_branch.n_branches > 0
            # radius of the child branches
            child_radius = parent_branch.target_radius * 0.5
            child_box = parent_branch.target_box * 0.5

            # count bodies per octant
            census!(cumulative_octant_census, system, parent_branch.bodies_index, parent_branch.target_center) # doesn't need to sort them here; just count them; the alternative is to save census data for EVERY CHILD BRANCH EACH GENERATION; then I save myself some effort at the expense of more memory allocation, as the octant_census would already be available; then again, the allocation might cost more than I save (which is what my intuition suggests)
            update_octant_accumulator!(cumulative_octant_census)

            # number of child branches
            if get_population(cumulative_octant_census) > leaf_size
                for i_octant in 1:8
                    if get_population(cumulative_octant_census, i_octant) > 0
                        bodies_index = get_bodies_index(cumulative_octant_census, parent_branch.bodies_index, i_octant)
                        child_center = get_child_center(parent_branch.target_center, parent_branch.target_box, i_octant)
                        child_branch, n_grandchildren, i_leaf = Branch(system, sort_index, octant_container, buffer, sort_index_buffer, i_first_branch, bodies_index, child_center, child_center, child_radius, child_radius, child_box, child_box, i_parent, i_leaf, leaf_size, expansion_order)
                        i_first_branch += n_grandchildren
                        push!(branches, child_branch)
                    end
                end
            end
        end
    end
    n_children = i_first_branch - length(branches) - 1 # the grandchildren of branches[parents_index]
    parents_index = parents_index[end]+1:length(branches) # the parents of the next generation
    return parents_index, n_children, i_leaf
end

function Branch(system, sort_index, octant_container, buffer, sort_index_buffer, i_first_branch, bodies_index, source_center, target_center, source_radius, target_radius, source_box, target_box, i_parent, i_leaf, leaf_size, expansion_order)
    # count bodies in each octant
    census!(octant_container, system, bodies_index, target_center)

    # cumsum
    update_octant_accumulator!(octant_container)
    # number of child branches
    n_children = get_n_children(octant_container, leaf_size)

    if n_children > 0
        # get beginning index of sorted bodies
        octant_beginning_index!(octant_container, bodies_index)

        # sort bodies into octants
        sort_bodies!(system, sort_index, octant_container, buffer, sort_index_buffer, bodies_index, target_center)
    end

    # get child branch information
    branch_index = i_first_branch : i_first_branch + n_children - 1
    n_branches = length(branch_index)

    if n_branches == 0
        i_leaf_index = i_leaf
        i_leaf += 1
    else
        i_leaf_index = -1
    end

    return Branch(bodies_index, n_branches, branch_index, i_parent, i_leaf_index, source_center, target_center, source_radius, target_radius, source_box, target_box, expansion_order), n_children, i_leaf
end

function Branch(bodies_index::UnitRange, n_branches, branch_index, i_parent, i_leaf_index, source_center, target_center, source_radius, target_radius, source_box, target_box, expansion_order)
    return SingleBranch(bodies_index, n_branches, branch_index, i_parent, i_leaf_index, source_center, target_center, source_radius, target_radius, source_box, target_box, initialize_expansion(expansion_order, typeof(source_radius)), initialize_expansion(expansion_order, typeof(source_radius)), initialize_harmonics(expansion_order, typeof(source_radius)), ReentrantLock())
end

function Branch(bodies_index, n_branches, branch_index, i_parent, i_leaf_index, source_center, target_center, source_radius, target_radius, source_box, target_box, expansion_order)
    return MultiBranch(bodies_index, n_branches, branch_index, i_parent, i_leaf_index, source_center, target_center, source_radius, target_radius, source_box, target_box, initialize_expansion(expansion_order, typeof(source_radius)), initialize_expansion(expansion_order, typeof(source_radius)), initialize_expansion(expansion_order, typeof(source_radius)), ReentrantLock())
end

@inline get_body_positions(system, bodies_index::UnitRange) = (system[i,POSITION] for i in bodies_index)

@inline function get_octant(position, center)
    octant = 1
    position[1] > center[1] && (octant += 1)
    position[2] > center[2] && (octant += 2)
    position[3] > center[3] && (octant += 4)
    return octant
end

@inline get_population(cumulative_octant_census::AbstractVector, i_octant) = i_octant == 1 ? cumulative_octant_census[1] : cumulative_octant_census[i_octant] - cumulative_octant_census[i_octant-1]

@inline get_population(cumulative_octant_census::AbstractMatrix, i_octant) = i_octant == 1 ? sum(cumulative_octant_census[:,1]) : sum(cumulative_octant_census[:,i_octant]) - sum(cumulative_octant_census[:,i_octant-1])

@inline get_population(cumulative_octant_census::AbstractVector) = cumulative_octant_census[end]

@inline get_population(cumulative_octant_census::AbstractMatrix) = sum(cumulative_octant_census[:,end])

@inline function get_child_center(parent_center, parent_target_box::SVector, i_octant)
    delta = parent_target_box[1] * 0.5
    i_octant -= 1
    dx = iseven(i_octant) ? -delta : delta
    delta = parent_target_box[2] * 0.5
    i_octant >>= 1
    dy = iseven(i_octant) ? -delta : delta
    delta = parent_target_box[3] * 0.5
    i_octant >>= 1
    dz = iseven(i_octant) ? -delta : delta
    child_center = SVector{3,eltype(parent_center)}(parent_center[1] + dx, parent_center[2] + dy, parent_center[3] + dz)
end

function get_octant_container(systems::Tuple)
    census = MMatrix{length(systems),8,Int64}(undef)
    for (i,system) in enumerate(systems)
        census[i,1] = get_n_bodies(system)
    end
    return census
end

function get_octant_container(system)
    census = MVector{8,Int64}(undef)
    census[1] = get_n_bodies(system)
    return census
end

function census!(octant_container::AbstractVector, system, bodies_index, center)
    octant_container .= zero(eltype(octant_container))
    for position in get_body_positions(system, bodies_index)
        i_octant = get_octant(position, center)
        octant_container[i_octant] += 1
    end
end

function census!(octant_container::AbstractMatrix, systems, bodies_indices, center)
    octant_container .= zero(eltype(octant_container))
    for (i_system, (system, bodies_index)) in enumerate(zip(systems, bodies_indices))
        census!(view(octant_container,i_system,:), system, bodies_index, center)
    end
end

@inline update_octant_accumulator!(octant_population::AbstractMatrix) = cumsum!(octant_population, octant_population, dims=2)

@inline update_octant_accumulator!(octant_population::AbstractVector) = cumsum!(octant_population, octant_population)

@inline function octant_beginning_index!(cumulative_octant_census::AbstractVector, bodies_index::UnitRange)
    if length(bodies_index) > 0
        for i_octant in 8:-1:2
            cumulative_octant_census[i_octant] = cumulative_octant_census[i_octant-1] + bodies_index[1]
        end
        cumulative_octant_census[1] = bodies_index[1]
    end
    return cumulative_octant_census
end

@inline function octant_beginning_index!(cumulative_octant_census::AbstractMatrix, bodies_indices::AbstractVector)
    for (i_system,bodies_index) in enumerate(bodies_indices)
        octant_beginning_index!(view(cumulative_octant_census,i_system,:), bodies_index)
    end
    return cumulative_octant_census
end

@inline function get_bodies_index(cumulative_octant_census::AbstractVector, parent_bodies_index::UnitRange, i_octant)
    if length(parent_bodies_index) > 0
        first_offset = i_octant == 1 ? 0 : cumulative_octant_census[i_octant-1]
        bodies_index = parent_bodies_index[1] + first_offset : parent_bodies_index[1] + cumulative_octant_census[i_octant] - 1
    else
        bodies_index = 1:0
    end
    return bodies_index
end

@inline function get_bodies_index(cumulative_octant_census::AbstractMatrix, parent_bodies_indices::AbstractVector, i_octant)
    n_systems = size(cumulative_octant_census,1)
    bodies_index = SVector{n_systems,UnitRange{Int64}}([get_bodies_index(view(cumulative_octant_census,i_system,:), parent_bodies_indices[i_system], i_octant) for i_system in 1:n_systems])
    # TODO why does this require square brackets? ah, because StaticArrays tries to collect the unit range ONLY IF THERE IS ONLY ONE;
    return bodies_index
end

@inline get_bodies_index(system) = 1:get_n_bodies(system)

@inline function get_bodies_index(systems::Tuple)
    n_systems = length(systems)
    return SVector{n_systems,UnitRange{Int64}}([get_bodies_index(system) for system in systems])
    # TODO why does this require square brackets? ah, because StaticArrays tries to collect the unit range ONLY IF THERE IS ONLY ONE;
end

#####
##### create buffers; overload for SortWrappers---in case we cannot sort bodies in place,
#####                 wrap the body in a SortWrapper object with an index which can be sorted
@inline function get_buffer(system::SortWrapper)
    return nothing
end

@inline function get_buffer(system)
    # return Vector{Int64}(undef,length(system))
    # [buffer_element(system) for _ in 1:length(system)]
    buffer = [buffer_element(system)]
    resize!(buffer,get_n_bodies(system))
    return buffer
end

@inline function get_buffer(systems::Tuple)
    return Tuple(get_buffer(system) for system in systems)
end

@inline function get_sort_index(system::SortWrapper)
    return system.index
end

@inline function get_sort_index(system)
    return collect(1:get_n_bodies(system))
end

@inline function get_sort_index(systems::Tuple)
    return Tuple(get_sort_index(system) for system in systems)
end

@inline function get_sort_index_buffer(system) # need not be overloaded for SortWrapper as it will be the same
    return Vector{Int64}(undef,get_n_bodies(system))
end

@inline function get_sort_index_buffer(systems::Tuple)
    return Tuple(get_sort_index_buffer(system) for system in systems)
end

#####
##### determine the number of descendants
#####
@inline function get_n_children(cumulative_octant_census, leaf_size)
    n_children = 0
    if get_population(cumulative_octant_census) > leaf_size
        for i_octant in 1:8
            get_population(cumulative_octant_census,i_octant) > 0 && (n_children += 1)
        end
    end
    return n_children
end

#####
##### sort bodies into the octree
#####
function sort_bodies!(system::SortWrapper, sort_index, octant_indices::AbstractVector, buffer, sort_index_buffer, bodies_index::UnitRange, center)
    # sort indices
    for i_body in bodies_index
        i_octant = get_octant(system[i_body,POSITION], center)
        sort_index_buffer[octant_indices[i_octant]] = sort_index[i_body]
        octant_indices[i_octant] += 1
    end
    # place buffers
    sort_index[bodies_index] .= view(sort_index_buffer,bodies_index)
end

function sort_bodies!(system, sort_index, octant_indices::AbstractVector, buffer, sort_index_buffer, bodies_index::UnitRange, center)
    # sort indices
    for i_body in bodies_index
        i_octant = get_octant(system[i_body, POSITION], center)
        buffer[octant_indices[i_octant]] = system[i_body, BODY]
        sort_index_buffer[octant_indices[i_octant]] = sort_index[i_body]
        octant_indices[i_octant] += 1
    end
    # place buffers
    for i_body in bodies_index
        system[i_body, BODY] = buffer[i_body]
    end
    sort_index[bodies_index] .= view(sort_index_buffer,bodies_index)
end

function sort_bodies!(systems, sort_indices, octant_indices::AbstractMatrix, buffers, sort_index_buffers, bodies_indices::AbstractVector, center)
    for (i_system, (system, sort_index, buffer, sort_index_buffer, bodies_index)) in enumerate(zip(systems, sort_indices, buffers, sort_index_buffers, bodies_indices))
        sort_bodies!(system, sort_index, view(octant_indices,i_system,:), buffer, sort_index_buffer, bodies_index, center)
    end
end

#####
##### invert the sort permutation for undoing the sort operation
#####
function invert_index!(inverse_sort_index::AbstractVector{Int64}, sort_index::AbstractVector{Int64})
    for i_body in eachindex(sort_index)
        inverse_sort_index[sort_index[i_body]] = i_body
    end
end

function invert_index!(inverse_sort_indices::Tuple, sort_indices::Tuple)
    for (inverse_sort_index, sort_index) in zip(inverse_sort_indices, sort_indices)
        invert_index!(inverse_sort_index, sort_index)
    end
end

#####
##### undo/redo the sort operation used to create the octree
#####
"""
Undoes the sort operation performed by the tree.
"""
function unsort!(systems::Tuple, tree::MultiTree)
    for (system, buffer, inverse_sort_index) in zip(systems, tree.buffers, tree.inverse_sort_index_list)
        unsort!(system, buffer, inverse_sort_index)
    end
end

function unsort!(system, tree::SingleTree)
    unsort!(system, tree.buffer, tree.inverse_sort_index)
end

function unsort!(system, buffer, tree::Tree)
    unsort!(system, buffer, tree.inverse_sort_index)
end

@inline function unsort!(system, buffer, inverse_sort_index)
    for i_body in 1:get_n_bodies(system)
        buffer[i_body] = system[inverse_sort_index[i_body], BODY]
    end
    for i_body in 1:get_n_bodies(system)
        system[i_body, BODY] = buffer[i_body]
    end
end

@inline function unsort!(system::SortWrapper, buffer, inverse_sort_index)
    system.index .= view(system.index,inverse_sort_index)
end

"""
Performs the same sort operation as the tree. (Undoes `unsort!` operation.)
"""
function resort!(systems::Tuple, tree::MultiTree)
    for (system, buffer, inverse_sort_index) in zip(systems, tree.buffers, tree.inverse_sort_index_list)
        resort!(system, buffer, inverse_sort_index)
    end
end

function resort!(system, tree::SingleTree)
    resort!(system, tree.buffer, tree.sort_index)
end

function resort!(system, buffer, tree::Tree)
    resort!(system, buffer, tree.sort_index)
end

function resort!(system, buffer, sort_index)
    for i_body in 1:get_n_bodies(system)
        buffer[i_body] = system[sort_index[i_body], BODY]
    end
    for i_body in 1:get_n_bodies(system)
        system[i_body, BODY] = buffer[i_body]
    end
end

@inline function unsorted_index_2_sorted_index(i_unsorted, tree::SingleTree)
    return tree.inverse_sort_index[i_unsorted]
end

@inline function unsorted_index_2_sorted_index(i_unsorted, i_system, tree::MultiTree)
    return tree.inverse_sort_index_list[i_system][i_unsorted]
end

@inline function sorted_index_2_unsorted_index(i_sorted, tree::SingleTree)
    return tree.sort_index[i_sorted]
end

@inline function sorted_index_2_unsorted_index(i_unsorted, i_system, tree::MultiTree)
    return tree.sort_index_list[i_system][i_unsorted]
end

#--- find the center of a (group of) system(s) of bodies of zero radius ---#

@inline function max_xyz(x_min, x_max, y_min, y_max, z_min, z_max, x, y, z)
    x_min = min(x_min, x)
    x_max = max(x_max, x)
    y_min = min(y_min, y)
    y_max = max(y_max, y)
    z_min = min(z_min, z)
    z_max = max(z_max, z)

    return x_min, x_max, y_min, y_max, z_min, z_max
end

@inline function max_xyz(x_min, x_max, y_min, y_max, z_min, z_max, system, bodies_index)
    for i in bodies_index
        x, y, z = system[i,Position()]
        x_min, x_max, y_min, y_max, z_min, z_max = max_xyz(x_min, x_max, y_min, y_max, z_min, z_max, x, y, z)
    end

    return x_min, x_max, y_min, y_max, z_min, z_max
end

@inline function max_xyz(x_min, x_max, y_min, y_max, z_min, z_max, branches::Vector{<:Branch}, child_index)

    # loop over child branches
    for i_child in child_index

        # extract child branch
        child_branch = branches[i_child]
        cx, cy, cz = child_branch.target_center
        dx, dy, dz = child_branch.target_box

        # get bounding box
        x_min = min(x_min, cx-dx)
        x_max = max(x_max, cx+dx)
        y_min = min(y_min, cy-dy)
        y_max = max(y_max, cy+dy)
        z_min = min(z_min, cz-dz)
        z_max = max(z_max, cz+dz)
    end

    return x_min, x_max, y_min, y_max, z_min, z_max
end

@inline function get_center_box(x_min, x_max, y_min, y_max, z_min, z_max)
    center = SVector{3}((x_max+x_min)*0.5, (y_max+y_min)*0.5, (z_max+z_min)*0.5)
    bounding_box = SVector{3}(x_max-center[1], y_max-center[2], z_max-center[3])
    dx, dy, dz = bounding_box
    return center, bounding_box
end

function center_box(systems)
    bodies_indices = get_bodies_index(systems)
    return center_box(systems, bodies_indices)
end

function center_box(systems::Tuple, bodies_indices)
    x_min, y_min, z_min = first_body_position(systems, bodies_indices)
    x_max, y_max, z_max = x_min, y_min, z_min
    for (system, bodies_index) in zip(systems, bodies_indices)
        x_min, x_max, y_min, y_max, z_min, z_max = max_xyz(x_min, x_max, y_min, y_max, z_min, z_max, system, bodies_index)
    end

    return get_center_box(x_min, x_max, y_min, y_max, z_min, z_max)
end

function center_box(system, bodies_index)
    x_min, y_min, z_min = first_body_position(system, bodies_index)
    x_max, y_max, z_max = x_min, y_min, z_min
    x_min, x_max, y_min, y_max, z_min, z_max = max_xyz(x_min, x_max, y_min, y_max, z_min, z_max, system, bodies_index)
    return get_center_box(x_min, x_max, y_min, y_max, z_min, z_max)
end

@inline function center_box(branches::Vector{<:Branch}, child_index)

    # smallest bounding rectangle
    first_branch = branches[child_index[1]]
    x_min, y_min, z_min = first_branch.target_center
    x_max = x_min
    y_max = y_min
    z_max = z_min
    x_min, x_max, y_min, y_max, z_min, z_max = max_xyz(x_min, x_max, y_min, y_max, z_min, z_max, branches, child_index)

    # find center
    center = SVector{3}((x_min+x_max)*0.5, (y_min+y_max)*0.5, (z_min+z_max)*0.5)
    bounding_box = SVector{3}(x_max - center[1], y_max - center[2], z_max - center[3])

    return center, bounding_box
end

@inline function source_max_xyz(x_max, y_max, z_max, x, y, z)
    x_max = max(x_max, x)
    y_max = max(y_max, y)
    z_max = max(z_max, z)

    return x_max, y_max, z_max
end

@inline function source_min_xyz(x_min, y_min, z_min, x, y, z)
    x_min = min(x_min, x)
    y_min = min(y_min, y)
    z_min = min(z_min, z)

    return x_min, y_min, z_min
end

@inline function source_max_xyz(x_min, x_max, y_min, y_max, z_min, z_max, system, bodies_index)
    for i in bodies_index
        x, y, z = system[i,Position()]
        r = system[i,Radius()]
        x_min, y_min, z_min = source_min_xyz(x_min, y_min, z_min, x-r, y-r, z-r)
        x_max, y_max, z_max = source_max_xyz(x_max, y_max, z_max, x+r, y+r, z+r)
    end

    return x_min, x_max, y_min, y_max, z_min, z_max
end

@inline function source_max_xyz(x_min, x_max, y_min, y_max, z_min, z_max, branches::Vector{<:Branch}, child_index)

    # loop over child branches
    for i_child in child_index

        # extract child branch
        child_branch = branches[i_child]
        cx, cy, cz = child_branch.source_center
        dx, dy, dz = child_branch.source_box

        # get bounding box
        x_min = min(x_min, cx-dx)
        x_max = max(x_max, cx+dx)
        y_min = min(y_min, cy-dy)
        y_max = max(y_max, cy+dy)
        z_min = min(z_min, cz-dz)
        z_max = max(z_max, cz+dz)
    end

    return x_min, x_max, y_min, y_max, z_min, z_max
end

function source_center_box(systems::Tuple, bodies_indices)
    x_min, y_min, z_min = first_body_position(systems, bodies_indices)
    x_max, y_max, z_max = x_min, y_min, z_min
    for (system, bodies_index) in zip(systems, bodies_indices)
        x_min, x_max, y_min, y_max, z_min, z_max = source_max_xyz(x_min, x_max, y_min, y_max, z_min, z_max, system, bodies_index)
    end

    return get_center_box(x_min, x_max, y_min, y_max, z_min, z_max)
end

function source_center_box(system, bodies_index)
    x_min, y_min, z_min = first_body_position(system, bodies_index)
    x_max, y_max, z_max = x_min, y_min, z_min
    x_min, x_max, y_min, y_max, z_min, z_max = source_max_xyz(x_min, x_max, y_min, y_max, z_min, z_max, system, bodies_index)
    return get_center_box(x_min, x_max, y_min, y_max, z_min, z_max)
end

@inline function source_center_box(branches::Vector{<:Branch}, child_index)

    # smallest bounding rectangle
    first_branch = branches[child_index[1]]
    x_min, y_min, z_min = first_branch.source_center
    x_max = x_min
    y_max = y_min
    z_max = z_min
    x_min, x_max, y_min, y_max, z_min, z_max = source_max_xyz(x_min, x_max, y_min, y_max, z_min, z_max, branches, child_index)

    # find center
    source_center = SVector{3}((x_min+x_max)*0.5, (y_min+y_max)*0.5, (z_min+z_max)*0.5)
    bounding_box = SVector{3}(x_max - source_center[1], y_max - source_center[2], z_max - source_center[3])

    return source_center, bounding_box
end

#------- shrinking method for bodies of non-zero radius -------#

@inline function first_body_position(system, bodies_index)
    return system[bodies_index[1],POSITION]
end

@inline function first_body_position(systems::Tuple, bodies_indices)
    for (system, bodies_index) in zip(systems, bodies_indices)
        length(bodies_index) > 0 && (return system[bodies_index[1],POSITION])
    end
end

@inline get_n_bodies(bodies_index::UnitRange) = length(bodies_index)

@inline function get_n_bodies(bodies_indices::AbstractVector{<:UnitRange})
    n_bodies = 0
    for bodies_index in bodies_indices
        n_bodies += get_n_bodies(bodies_index)
    end
    return n_bodies
end

@inline get_n_bodies(branch::Branch) = get_n_bodies(branch.bodies_index)

@inline get_n_bodies_vec(branch::SingleBranch) = length(branch.bodies_index)

@inline get_n_bodies_vec(branch::MultiBranch) = SVector{Int}(length(bodies_index) for bodies_index in branch.bodies_index)


@inline function shrink_radius(target_radius, source_radius, target_center, source_center, systems::Tuple, bodies_indices)
    # loop over systems
    for (system, bodies_index) in zip(systems, bodies_indices)
        target_radius, source_radius = shrink_radius(target_radius, source_radius, target_center, source_center, system, bodies_index)
    end
    return target_radius, source_radius
end

@inline function shrink_radius(target_radius, source_radius, target_center, source_center, system, bodies_index)
    # extract target/source centers
    cxt, cyt, czt = target_center
    cxs, cys, czs = source_center

    # loop over all bodies
    for i_body in bodies_index

        # extract body position and size
        x, y, z = system[i_body,Position()]
        body_radius = system[i_body,Radius()]

        # max target radius
        dx = x - cxt
        dy = y - cyt
        dz = z - czt
        distance = sqrt(dx*dx + dy*dy + dz*dz)
        target_radius = max(target_radius, distance)

        # max source radius
        dx = x - cxs
        dy = y - cys
        dz = z - czs
        distance = sqrt(dx*dx + dy*dy + dz*dz)
        source_radius = max(source_radius, distance + body_radius)

    end

    return target_radius, source_radius
end

@inline function shrink_radius(target_center, source_center, system, bodies_index)
    # initialize values
    target_radius = zero(eltype(target_center))
    source_radius = zero(eltype(target_center))

    # get radii
    target_radius, source_radius = shrink_radius(target_radius, source_radius, target_center, source_center, system, bodies_index)

    return target_radius, source_radius
end

@inline function shrink_radius(target_center, source_center, branches::Vector{<:Branch}, child_index)
    # initialize values
    target_radius = zero(eltype(target_center))
    source_radius = zero(eltype(target_center))
    first_child = branches[child_index[1]]

    # get radii
    target_radius, source_radius = shrink_radius(target_radius, source_radius, target_center, source_center, branches, child_index)

    return target_radius, source_radius
end

@inline function shrink_radius(target_radius, source_radius, target_center, source_center, branches::Vector{<:Branch}, child_index)
    # extract target/source centers
    cxt, cyt, czt = target_center
    cxs, cys, czs = source_center

    # loop over all bodies
    for i_child in child_index

        # extract body position and size
        child_branch = branches[i_child]
        xt, yt, zt = child_branch.target_center
        xs, ys, zs = child_branch.source_center
        child_target_radius = child_branch.target_radius
        child_source_radius = child_branch.source_radius

        # max target radius
        dx = xt - cxt
        dy = yt - cyt
        dz = zt - czt
        distance = sqrt(dx*dx + dy*dy + dz*dz)
        target_radius = max(target_radius, distance + child_target_radius)

        # max source radius
        dx = xs - cxs
        dy = ys - cys
        dz = zs - czs
        distance = sqrt(dx*dx + dy*dy + dz*dz)
        source_radius = max(source_radius, distance + child_source_radius)

    end

    return target_radius, source_radius
end

function shrink_leaf!(branches, i_branch, system)
    # unpack
    branch = branches[i_branch]
    bodies_index = branch.bodies_index

    # recenter and target box
    new_target_center, new_target_box = center_box(system, bodies_index)
    new_source_center, new_source_box = source_center_box(system, bodies_index)

    # shrink radii and create source box
    new_target_radius, new_source_radius = shrink_radius(new_target_center, new_source_center, system, bodies_index)

    # replace branch
    replace_branch!(branches, i_branch, new_source_center, new_target_center, new_source_radius, new_target_radius, new_source_box, new_target_box)

end

"""
Computes the smallest bounding box to completely bound all child boxes.

Shrunk radii are merely the distance from the center to the corner of the box.
"""
function shrink_branch!(branches, i_branch, child_index)

    # recenter and target box
    new_target_center, new_target_box = center_box(branches, child_index)
    new_source_center, new_source_box = source_center_box(branches, child_index)

    # shrink radii and create source box
    new_target_radius, new_source_radius = shrink_radius(new_target_center, new_source_center, branches, child_index)

    # replace branch
    replace_branch!(branches, i_branch, new_source_center, new_target_center, new_source_radius, new_target_radius, new_source_box, new_target_box)
end

function shrink_recenter!(branches, levels_index, system)
    for i_level in length(levels_index):-1:1 # start at the bottom level
        level_index = levels_index[i_level]
        for i_branch in level_index
            branch = branches[i_branch]
            if branch.n_branches == 0 # leaf
                shrink_leaf!(branches, i_branch, system)
            else
                shrink_branch!(branches, i_branch, branch.branch_index)
            end
        end
    end
end

#--- accumulate charge ---#

"""
    accumulate_charge!(branches, systems, branch_index)

Computes the sum of the absolute value of body source strengths, as well as the sum of the L2 norm of the dipole strengths, and store them in each branch, inheriting values from child branches if available.

"""
function accumulate_charge!(branches::AbstractVector{<:Branch{TF}}, systems, branch_index) where TF

    # loop over branches
    for i_branch in branch_index
        branch = branches[i_branch]

        # inherit from child branches
        if branch.n_branches > 0
            charge = zero(TF)
            dipole = zero(TF)
            for i_child in branch.branch_index
                child = branches[i_child]
                charge += child.charge
                dipole += child.dipole
            end

            # replace branch
            replace_branch!(branches, i_branch, charge, dipole)

        # compute based on bodies directly
        else
            accumulate_charge_bodies!(branches, i_branch, systems)
        end
    end
end

"""
    accumulate_charge!(branch, systems)

Computes the sum of the absolute value of body source strengths, as well as the sum of the L2 norm of the dipole strengths, based on source bodies directly, and store them in each branch.

"""
function accumulate_charge_bodies!(branches::Vector{SingleBranch{TF}}, i_branch, system) where TF
    # reset multipole expansion, just in case
    branch = branches[i_branch]
    branch.multipole_expansion .= zero(TF)

    # only need through p=1 for this
    expansion_order = Val(1)

    # accumulate charge and dipole
    charge = zero(TF)
    dipole = zero(TF)

    # loop over bodies
    for i_body in branch.bodies_index

        # reset required expansions
        branch.multipole_expansion[1:2,1:2,1:3] .= zero(TF)

        # multipole coefficients
        body_to_multipole!(system, branch, i_body:i_body, branch.harmonics, expansion_order)

        # get effective charge
        charge += abs(branch.multipole_expansion[1,1,1])

        # get effective dipole moment
        qz = branch.multipole_expansion[1,1,2]
        qy = branch.multipole_expansion[1,1,3] * 2.0
        qx = branch.multipole_expansion[2,1,3] * 2.0
        dipole += sqrt(qx*qx + qy*qy + qz*qz)
    end

    # reset multipole expansion again
    branch.multipole_expansion .= zero(TF)

    # replace branch
    replace_branch!(branches, i_branch, charge, dipole)

end

function accumulate_charge_bodies!(branches::Vector{MultiBranch{TF,N}}, i_branch, systems) where {TF,N}
    # reset multipole expansion, just in case
    branch = branches[i_branch]
    branch.multipole_expansion .= zero(TF)

    # only need through p=1 for this
    expansion_order = Val(1)

    # accumulate charge and dipole
    charge = zero(TF)
    dipole = zero(TF)

    # loop over systems
    for i_system in 1:N
        system = systems[i_system]

        # loop over bodies
        for i_body in branch.bodies_index[i_system]

            # reset required expansions
            branch.multipole_expansion[1:2,1:2,1:3] .= zero(TF)

            # multipole coefficients
            body_to_multipole!(system, branch, i_body:i_body, branch.harmonics, expansion_order)

            # get effective charge
            charge += abs(branch.multipole_expansion[1,1,1])

            # get effective dipole moment
            qz = branch.multipole_expansion[1,1,2]
            qy = branch.multipole_expansion[1,1,3] * 2.0
            qx = branch.multipole_expansion[2,1,3] * 2.0
            dipole += sqrt(qx*qx + qy*qy + qz*qz)
        end
    end

    # reset multipole expansion again
    branch.multipole_expansion .= zero(TF)

    # replace branch
    replace_branch!(branches, i_branch, charge, dipole)

end

#--- helper function ---#

@inline function replace_branch!(branches::Vector{TB}, i_branch, new_source_center, new_target_center, new_source_radius, new_target_radius, new_source_box, new_target_box) where TB
    # (; bodies_index, n_branches, branch_index, i_parent, center, radius, multipole_expansion, local_expansion, harmonics, lock) = branch[]
    branch = branches[i_branch]
    bodies_index = branch.bodies_index
    n_branches = branch.n_branches
    branch_index = branch.branch_index
    i_parent = branch.i_parent
    i_leaf = branch.i_leaf
    charge = branch.charge
    dipole = branch.dipole
    error = branch.error
    multipole_expansion = branch.multipole_expansion
    local_expansion = branch.local_expansion
    harmonics = branch.harmonics
    lock = branch.lock
    branches[i_branch] = TB(bodies_index, n_branches, branch_index, i_parent, i_leaf, new_source_center, new_target_center, new_source_radius, new_target_radius, new_source_box, new_target_box, charge, dipole, error, multipole_expansion, local_expansion, harmonics, lock)
end

@inline function replace_branch!(branches::Vector{TB}, i_branch, new_charge, new_dipole) where TB
    # (; bodies_index, n_branches, branch_index, i_parent, center, radius, multipole_expansion, local_expansion, harmonics, lock) = branch[]
    branch = branches[i_branch]
    bodies_index = branch.bodies_index
    n_branches = branch.n_branches
    branch_index = branch.branch_index
    i_parent = branch.i_parent
    i_leaf = branch.i_leaf
    source_center = branch.source_center
    target_center = branch.target_center
    source_radius = branch.source_radius
    target_radius = branch.target_radius
    source_box = branch.source_box
    target_box = branch.target_box
    error = branch.error
    multipole_expansion = branch.multipole_expansion
    local_expansion = branch.local_expansion
    harmonics = branch.harmonics
    lock = branch.lock
    branches[i_branch] = TB(bodies_index, n_branches, branch_index, i_parent, i_leaf, source_center, target_center, source_radius, target_radius, source_box, target_box, new_charge, new_dipole, error, multipole_expansion, local_expansion, harmonics, lock)
end

function initialize_expansion(expansion_order, type=Float64)
    return zeros(type, 2, 2, ((expansion_order+1) * (expansion_order+2)) >> 1)
end

function initialize_velocity_n_m(expansion_order, type=Float64)
    p = expansion_order
    n_harmonics = harmonic_index(p,p)
    return zeros(type, 2, 3, n_harmonics)
end

function initialize_harmonics(expansion_order, type=Float64)
    p = expansion_order # +1 to allow for evaluate_multipole!()
    n_harmonics = harmonic_index(p,p)
    return zeros(type, 2, 2, n_harmonics)
end

function reset_expansions!(tree)
    T = eltype(tree.branches[1].multipole_expansion)
    for branch in tree.branches
        branch.multipole_expansion .= zero(T)
        branch.local_expansion .= zero(T)
    end
end
