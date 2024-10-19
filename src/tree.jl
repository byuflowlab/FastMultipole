# const BRANCH_TYPE = Float64
# global SHRINKING_OFFSET = .000001

const WARNING_FLAG_LEAF_SIZE = Array{Bool,0}(undef)
WARNING_FLAG_LEAF_SIZE[] = true

#####
##### tree constructor
#####
function Tree(system; expansion_order=7, leaf_size=100, n_divisions=20, shrink=false, recenter=false, allocation_safety_factor=1.0, estimate_cost=false, read_cost_file=false, write_cost_file=false)
    # ensure `system` isn't empty; otherwise return an empty tree
    if get_n_bodies(system) > 0
        # initialize variables
        octant_container = get_octant_container(system) # preallocate octant counter; records the total number of bodies in the first octant to start out, but can be reused to count bodies per octant later on
        cumulative_octant_census = get_octant_container(system) # for creating a cumsum of octant populations
        bodies_index = get_bodies_index(system)
        center, radius, target_box = center_radius_box(system)
        bx, by, bz = target_box
        source_box = SVector{6}(bx, bx, by, by, bz, bz)
        i_first_branch = 2
        buffer = get_buffer(system)
        sort_index = get_sort_index(system)
        sort_index_buffer = get_sort_index_buffer(system)

        # grow root branch
        root_branch, n_children, i_leaf = Branch(system, sort_index, octant_container, buffer, sort_index_buffer, i_first_branch, bodies_index, center, radius, radius, source_box, target_box, 0, 1, leaf_size, expansion_order) # even though no sorting needed for creating this branch, it will be needed later on; so `Branch` not ony_min creates the root_branch, but also sorts itself into octants and returns the number of children it will have so we can plan array size
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

        # # check depth
        # if n_children > 0
        #     error("leaf_size not reached; n_children = $n_children")
        # end

        # invert index
        invert_index!(sort_index_buffer, sort_index)
        inverse_sort_index = sort_index_buffer # reuse the buffer as the inverse index

        # shrink and recenter branches to account for bodies of nonzero radius
        if shrink
            shrink_recenter!(branches, levels_index, system, recenter)
        end

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
            child_source_box = parent_branch.source_box * 0.5
            child_target_box = parent_branch.target_box * 0.5

            # count bodies per octant
            census!(cumulative_octant_census, system, parent_branch.bodies_index, parent_branch.center) # doesn't need to sort them here; just count them; the alternative is to save census data for EVERY CHILD BRANCH EACH GENERATION; then I save myself some effort at the expense of more memory allocation, as the octant_census would already be available; then again, the allocation might cost more than I save (which is what my intuition suggests)
            update_octant_accumulator!(cumulative_octant_census)

            # number of child branches
            if get_population(cumulative_octant_census) > leaf_size
                for i_octant in 1:8
                    if get_population(cumulative_octant_census, i_octant) > 0
                        bodies_index = get_bodies_index(cumulative_octant_census, parent_branch.bodies_index, i_octant)
                        child_center = get_child_center(parent_branch.center, parent_branch.target_box, i_octant)
                        child_branch, n_grandchildren, i_leaf = Branch(system, sort_index, octant_container, buffer, sort_index_buffer, i_first_branch, bodies_index, child_center, child_radius, child_radius, child_source_box, child_target_box, i_parent, i_leaf, leaf_size, expansion_order)
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

function Branch(system, sort_index, octant_container, buffer, sort_index_buffer, i_first_branch, bodies_index, center, source_radius, target_radius, source_box, target_box, i_parent, i_leaf, leaf_size, expansion_order)
    # count bodies in each octant
    census!(octant_container, system, bodies_index, center)

    # cumsum
    update_octant_accumulator!(octant_container)
    # number of child branches
    n_children = get_n_children(octant_container, leaf_size)

    if n_children > 0
        # get beginning index of sorted bodies
        octant_beginning_index!(octant_container, bodies_index)

        # sort bodies into octants
        sort_bodies!(system, sort_index, octant_container, buffer, sort_index_buffer, bodies_index, center)
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

    return Branch(bodies_index, n_branches, branch_index, i_parent, i_leaf_index, center, source_radius, target_radius, source_box, target_box, expansion_order), n_children, i_leaf
end

function Branch(bodies_index::UnitRange, n_branches, branch_index, i_parent, i_leaf_index, center, source_radius, target_radius, source_box, target_box, expansion_order)
    return SingleBranch(bodies_index, n_branches, branch_index, i_parent, i_leaf_index, center, source_radius, target_radius, source_box, target_box, initialize_expansion(expansion_order, typeof(source_radius)), initialize_expansion(expansion_order, typeof(source_radius)), initialize_harmonics(expansion_order, typeof(source_radius)), ReentrantLock())
end

function Branch(bodies_index, n_branches, branch_index, i_parent, i_leaf_index, center, source_radius, target_radius, source_box, target_box, expansion_order)
    return MultiBranch(bodies_index, n_branches, branch_index, i_parent, i_leaf_index, center, source_radius, target_radius, source_box, target_box, initialize_expansion(expansion_order, typeof(source_radius)), initialize_expansion(expansion_order, typeof(source_radius)), initialize_expansion(expansion_order, typeof(source_radius)), ReentrantLock())
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
    delta = parent_target_box[1]
    i_octant -= 1
    dx = iseven(i_octant) ? -delta : delta
    delta = parent_target_box[2]
    i_octant >>= 1
    dy = iseven(i_octant) ? -delta : delta
    delta = parent_target_box[3]
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

#--- find the center and radius of a (group of) system(s) of bodies of zero radius ---#

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

@inline get_radius(dx,dy,dz) = sqrt(dx*dx + dy*dy + dz*dz) # assume spherical cells

@inline function get_center_radius_box(x_min, x_max, y_min, y_max, z_min, z_max)
    center = SVector{3}((x_max+x_min)*0.5, (y_max+y_min)*0.5, (z_max+z_min)*0.5)
    bounding_box = SVector{3}(x_max-center[1], y_max-center[2], z_max-center[3])
    radius = get_radius(bounding_box[1], bounding_box[2], bounding_box[3])
    return center, radius, bounding_box
end

function center_radius_box(systems)
    bodies_indices = get_bodies_index(systems)
    return center_radius_box(systems, bodies_indices)
end

function center_radius_box(systems::Tuple, bodies_indices)
    x_min, y_min, z_min = first_body_position(systems, bodies_indices)
    x_max, y_max, z_max = x_min, y_min, z_min
    for (system, bodies_index) in zip(systems, bodies_indices)
        x_min, x_max, y_min, y_max, z_min, z_max = max_xyz(x_min, x_max, y_min, y_max, z_min, z_max, system, bodies_index)
    end

    return get_center_radius_box(x_min, x_max, y_min, y_max, z_min, z_max)
end

function center_radius_box(system, bodies_index)
    x_min, y_min, z_min = first_body_position(system, bodies_index)
    x_max, y_max, z_max = x_min, y_min, z_min
    x_min, x_max, y_min, y_max, z_min, z_max = max_xyz(x_min, x_max, y_min, y_max, z_min, z_max, system, bodies_index)
    return get_center_radius_box(x_min, x_max, y_min, y_max, z_min, z_max)
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


@inline function shrink_radius_box(radius, dx_min, dx_max, dy_min, dy_max, dz_min, dz_max, center, systems::Tuple, bodies_indices)
    # loop over systems
    for (system, bodies_index) in zip(systems, bodies_indices)
        radius, dx_min, dx_max, dy_min, dy_max, dz_min, dz_max = shrink_radius_box(radius, dx_max, dx_min, dy_max, dy_min, dz_max, dz_min, center, system, bodies_index)
    end
    return radius, dx_min, dx_max, dy_min, dy_max, dz_min, dz_max
end

@inline function shrink_radius_box(radius, dx_max, dx_min, dy_max, dy_min, dz_max, dz_min, center, system, bodies_index)
    # extract center
    cx, cy, cz = center

    # loop over all bodies
    for i_body in bodies_index

        # extract body position and size
        x, y, z = system[i_body,Position()]
        body_radius = system[i_body,Radius()]

        # max radius
        dx = x - cx
        dy = y - cy
        dz = z - cz
        distance = sqrt(dx*dx + dy*dy + dz*dz)
        radius = max(radius, distance + body_radius)

        # bounding box
        dx_min = min(dx_min, dx-body_radius)
        dx_max = max(dx_max, dx+body_radius)
        dy_min = min(dy_min, dy-body_radius)
        dy_max = max(dy_max, dy+body_radius)
        dz_min = min(dz_min, dz-body_radius)
        dz_max = max(dz_max, dz+body_radius)
    end

    return radius, dx_min, dx_max, dy_min, dy_max, dz_min, dz_max
end

@inline function shrink_radius_box(center, system, bodies_index)
    # initialize values
    radius = zero(eltype(center))
    x, y, z = first_body_position(system, bodies_index)
    dx_max, dy_max, dz_max = x - center[1], y - center[2], z - center[3]
    dx_min, dy_min, dz_min = dx_max, dy_max, dz_max

    # max/min values
    radius, dx_min, dx_max, dy_min, dy_max, dz_min, dz_max = shrink_radius_box(radius, dx_min, dx_max, dy_min, dy_max, dz_min, dz_max, center, system, bodies_index)

    # form bounding box
    bounding_box = SVector{6}(-dx_min, dx_max, -dy_min, dy_max, -dz_min, dz_max)
    for i in 1:6
        @assert bounding_box[i] >= 0.0
    end

    return radius, bounding_box
end

@inline function replace_branch!(branches::Vector{TB}, i_branch, new_center, new_source_radius, new_target_radius, new_source_box, new_target_box) where TB
    # (; bodies_index, n_branches, branch_index, i_parent, center, radius, multipole_expansion, local_expansion, harmonics, lock) = branch[]
    branch = branches[i_branch]
    bodies_index = branch.bodies_index
    n_branches = branch.n_branches
    branch_index = branch.branch_index
    i_parent = branch.i_parent
    i_leaf = branch.i_leaf
    multipole_expansion = branch.multipole_expansion
    local_expansion = branch.local_expansion
    harmonics = branch.harmonics
    lock = branch.lock
    branches[i_branch] = TB(bodies_index, n_branches, branch_index, i_parent, i_leaf, new_center, new_source_radius, new_target_radius, new_source_box, new_target_box, multipole_expansion, local_expansion, harmonics, lock)
end

function shrink_leaf!(branches, i_branch, system, recenter)
    # unpack
    branch = branches[i_branch]
    bodies_index = branch.bodies_index

    # recenter and target box
    new_center, new_target_radius, new_target_box = center_radius_box(system, bodies_index)

    if !recenter
        new_center = branch.center
        new_target_box = branch.target_box
        bx, by, bz = new_target_box
        new_target_radius = sqrt(bx*bx+by*by+bz*bz)
    end

    # shrink radius and source box
    new_source_radius, new_source_box = shrink_radius_box(new_center, system, bodies_index)

    # replace branch
    replace_branch!(branches, i_branch, new_center, new_source_radius, new_target_radius, new_source_box, new_target_box)
end

@inline function smallest_target_box(x_min, x_max, y_min, y_max, z_min, z_max, branch, branches, child_index)

    # loop over child branches
    for i_child in child_index

        # extract child branch
        child_branch = branches[i_child]
        cx, cy, cz = child_branch.center
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

@inline function center_box_target_branch(branch::Branch, branches, child_index)

    # smallest bounding rectangle
    x_min, y_min, z_min = branch.center
    x_max = x_min
    y_max = y_min
    z_max = z_min
    x_min, x_max, y_min, y_max, z_min, z_max = smallest_target_box(x_min, x_max, y_min, y_max, z_min, z_max, branch, branches, child_index)

    # find center
    center = SVector{3}((x_min+x_max)*0.5, (y_min+y_max)*0.5, (z_min+z_max)*0.5)
    bounding_box = SVector{3}(x_max - center[1], y_max - center[2], z_max - center[3])

    return center, bounding_box
end

@inline function shrink_radius(radius, center, branch::Branch, branches, branch_index)
    for i_branch in branch_index
        child_branch = branches[i_branch]
        x, y, z = child_branch.center
        body_radius = child_branch.radius
        radius = max(radius, get_distance(x, y, z, center) + body_radius)
    end
    return radius
end

"""
returns coordinates on extrema of the box
"""
@inline function smallest_source_box(x_min, x_max, y_min, y_max, z_min, z_max, branch, branches, child_index)

    # loop over child branches
    for i_child in child_index

        # extract child branch
        child_branch = branches[i_child]
        cx, cy, cz = child_branch.center
        dx_min, dx_max, dy_min, dy_max, dz_min, dz_max = child_branch.source_box

        # get bounding box
        x_min = min(x_min, cx-dx_min)
        x_max = max(x_max, cx+dx_max)
        y_min = min(y_min, cy-dy_min)
        y_max = max(y_max, cy+dy_max)
        z_min = min(z_min, cz-dz_min)
        z_max = max(z_max, cz+dz_max)
    end

    return x_min, x_max, y_min, y_max, z_min, z_max
end

function center_box_source_branch(branch, branches, child_index)

    # smallest bounding rectangle
    cx, cy, cz = branch.center
    x_min, x_max, y_min, y_max, z_min, z_max = smallest_source_box(cx, cx, cy, cy, cz, cz, branch, branches, child_index)

    # should be entirely positive values
    bounding_box = SVector{6}(cx-x_min, x_max-cx, cy-y_min, y_max-cy, cz-z_min, z_max-cz)

    return bounding_box
end

"""
Computes the smallest bounding box to completely bound all child boxes.

Shrunk radii are merely the distance from the center to the corner of the box.
"""
function shrink_branch!(branches, i_branch, child_index, recenter)

    # recenter about targets
    branch = branches[i_branch]
    new_center, new_target_box = center_box_target_branch(branch, branches, child_index)
    if !recenter
        new_center = branch.center
        cx, cy, cz = new_center
        x_min, x_max, y_min, y_max, z_min, z_max = smallest_target_box(cx, cx, cy, cy, cz, cz, branch, branches, child_index)
        new_target_box = SVector{3}(max(x_max-cx,cx-x_min), max(y_max-cy,cy-y_min), max(z_max-cz,cz-z_min))
    end
    bx, by, bz = new_target_box
    new_target_radius = sqrt(bx*bx + by*by + bz*bz)

    # shrink source box
    new_source_box = center_box_source_branch(branch, branches, child_index)
    bx_max, bx_min, by_max, by_min, bz_max, bz_min = new_source_box
    bx = max(bx_max, bx_min)
    by = max(by_max, by_min)
    bz = max(bz_max, bz_min)
    new_source_radius = sqrt(bx*bx + by*by + bz*bz)

    # replace branch
    replace_branch!(branches, i_branch, new_center, new_source_radius, new_target_radius, new_source_box, new_target_box)
end

function shrink_recenter!(branches, levels_index, system, recenter)
    for i_level in length(levels_index):-1:1 # start at the bottom level
        level_index = levels_index[i_level]
        for i_branch in level_index
            branch = branches[i_branch]
            if branch.n_branches == 0 # leaf
                shrink_leaf!(branches, i_branch, system, recenter)
            else
                shrink_branch!(branches, i_branch, branch.branch_index, recenter)
            end
        end
    end
end

#--- helper function ---#

function initialize_expansion(expansion_order, type=Float64)
    return zeros(type, 2, 2, ((expansion_order+1) * (expansion_order+2)) >> 1)
end

function initialize_harmonics(expansion_order, type=Float64)
    p = expansion_order
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
