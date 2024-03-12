# const BRANCH_TYPE = Float64
# global SHRINKING_OFFSET = .000001
# away_from_center! is a bandaid fix

const WARNING_FLAG_N_PER_BRANCH = Array{Bool,0}(undef)
WARNING_FLAG_N_PER_BRANCH[] = true

#####
##### tree constructor
#####
function Tree(system; expansion_order=7, n_per_branch=100, ndivisions=7, scale_radius=1.00001, shrink_recenter=false, allocation_safety_factor=1.0, estimate_cost=false, read_cost_file=true, write_cost_file=false)
    # initialize variables
    octant_container = get_octant_container(system) # preallocate octant counter; records the total number of bodies in the first octant to start out, but can be reused to count bodies per octant later on
    cumulative_octant_census = get_octant_container(system) # for creating a cumsum of octant populations
    bodies_index = get_bodies_index(system)
    center, radius = center_radius(system; scale_radius=scale_radius)
    i_first_branch = 2
    buffer = get_buffer(system)
    sort_index = get_sort_index(system)
    sort_index_buffer = get_sort_index_buffer(system)
    
    # grow root branch
    root_branch, n_children = Branch(system, sort_index, octant_container, buffer, sort_index_buffer, i_first_branch, bodies_index, center, radius, 0, n_per_branch, expansion_order) # even though no sorting needed for creating this branch, it will be needed later on; so `Branch` not ony_min creates the root_branch, but also sorts itself into octants and returns the number of children it will have so we can plan array size
    branches = [root_branch] # this first branch will already have its child branches encoded
    # estimated_n_branches = estimate_n_branches(system, n_per_branch, allocation_safety_factor)
    # sizehint!(branches, estimated_n_branches)
    parents_index = 1:1

    # grow branches
    levels_index = [parents_index] # store branches at each level
    for i_divide in 1:ndivisions
        if n_children > 0
            parents_index, n_children = child_branches!(branches, system, sort_index, buffer, sort_index_buffer, n_per_branch, parents_index, cumulative_octant_census, octant_container, n_children, expansion_order)
            push!(levels_index, parents_index)
        end
    end
    
    # check depth
    if n_children > 0
        n_children_prewhile = n_children
        while n_children > 0
            parents_index, n_children = child_branches!(branches, system, sort_index, buffer, sort_index_buffer, n_per_branch, parents_index, cumulative_octant_census, octant_container, n_children, expansion_order)
            push!(levels_index, parents_index)
        end
        if WARNING_FLAG_N_PER_BRANCH[]
            @warn "n_per_branch not reached in for loop, so while loop used to build octree; to improve performance, increase `ndivisions` > $(length(levels_index))"
            WARNING_FLAG_N_PER_BRANCH[] = false
        end
    end

    # # check depth
    # if n_children > 0
    #     error("n_per_branch not reached; n_children = $n_children")
    # end

    # invert index
    invert_index!(sort_index_buffer, sort_index)
    inverse_sort_index = sort_index_buffer # reuse the buffer as the inverse index

    # shrink and recenter branches to account for bodies of nonzero radius
    if shrink_recenter
        shrink_recenter!(branches, levels_index, system)
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
    # cost_parameters = Threads.nthreads() > 1 ? direct_cost_estimate(system, n_per_branch) : dummy_direct_cost_estimate(system, n_per_branch)

    # assemble tree
    tree = Tree(branches, levels_index, leaf_index, sort_index, inverse_sort_index, buffer, Val(expansion_order), n_per_branch)#, cost_parameters)

    return tree
end

Tree(branches::Vector{<:SingleBranch}, levels_index, leaf_index, sort_index, inverse_sort_index, buffer, expansion_order, n_per_branch) = 
    SingleTree(branches, levels_index, leaf_index, sort_index, inverse_sort_index, buffer, expansion_order, n_per_branch)#, cost_parameters)

Tree(branches::Vector{<:MultiBranch}, levels_index, leaf_index, sort_index, inverse_sort_index, buffer, expansion_order, n_per_branch) = 
    MultiTree(branches, levels_index, leaf_index, sort_index, inverse_sort_index, buffer, expansion_order, n_per_branch)#, cost_parameters)

@inline function get_n_bodies(systems::Tuple)
    n_bodies = 0
    for system in systems
        n_bodies += get_n_bodies(system)
    end
    return n_bodies
end

function estimate_n_branches(system, n_per_branch, allocation_safety_factor)
    n_bodies = get_n_bodies(system)
    estimated_n_divisions = Int(ceil(log(8,n_bodies/n_per_branch)))
    estimated_n_branches = div(8^(estimated_n_divisions+1) - 1,7)
    return Int(ceil(estimated_n_branches * allocation_safety_factor))
end

function child_branches!(branches, system, sort_index, buffer, sort_index_buffer, n_per_branch, parents_index, cumulative_octant_census, octant_container, n_children, expansion_order)
    i_first_branch = parents_index[end] + n_children + 1
    for i_parent in parents_index
        parent_branch = branches[i_parent]
        if parent_branch.n_branches > 0
            # radius of the child branches
            child_radius = parent_branch.radius / 2.0

            # count bodies per octant
            census!(cumulative_octant_census, system, parent_branch.bodies_index, parent_branch.center) # doesn't need to sort them here; just count them; the alternative is to save census data for EVERY CHILD BRANCH EACH GENERATION; then I save myself some effort at the expense of more memory allocation, as the octant_census would already be available; then again, the allocation might cost more than I save (which is what my intuition suggests)
            update_octant_accumulator!(cumulative_octant_census)

            # @show cumulative_octant_census parents_index get_population(cumulative_octant_census)
            
            # number of child branches
            if get_population(cumulative_octant_census) > n_per_branch
                for i_octant in 1:8
                    if get_population(cumulative_octant_census, i_octant) > 0  
                        bodies_index = get_bodies_index(cumulative_octant_census, parent_branch.bodies_index, i_octant)
                        child_center = get_child_center(parent_branch.center, parent_branch.radius, i_octant)
                        child_branch, n_grandchildren = Branch(system, sort_index, octant_container, buffer, sort_index_buffer, i_first_branch, bodies_index, child_center, child_radius, i_parent, n_per_branch, expansion_order)
                        i_first_branch += n_grandchildren
                        push!(branches, child_branch)
                    end
                end
            end
        end
    end
    n_children = i_first_branch - length(branches) - 1 # the grandchildren of branches[parents_index]
    parents_index = parents_index[end]+1:length(branches) # the parents of the next generation
    return parents_index, n_children
end

function Branch(system, sort_index, octant_container, buffer, sort_index_buffer, i_first_branch, bodies_index, center, radius, i_parent, n_per_branch, expansion_order)
    # count bodies in each octant
    census!(octant_container, system, bodies_index, center)
    
    # cumsum
    update_octant_accumulator!(octant_container)
    # number of child branches
    n_children = get_n_children(octant_container, n_per_branch)
    
    if n_children > 0
        # get beginning index of sorted bodies
        octant_beginning_index!(octant_container, bodies_index)

        # sort bodies into octants
        sort_bodies!(system, sort_index, octant_container, buffer, sort_index_buffer, bodies_index, center)
    end

    # get child branch information
    branch_index = i_first_branch : i_first_branch + n_children - 1
    n_branches = length(branch_index)

    return Branch(bodies_index, n_branches, branch_index, i_parent, center, radius, expansion_order), n_children
end

function Branch(bodies_index::UnitRange, n_branches, branch_index, i_parent, center, radius, expansion_order)
    return SingleBranch(bodies_index, n_branches, branch_index, i_parent, center, radius, initialize_expansion(expansion_order, typeof(radius)), initialize_expansion(expansion_order, typeof(radius)), initialize_harmonics(expansion_order, typeof(radius)), initialize_ML(expansion_order, typeof(radius)), ReentrantLock())
end

function Branch(bodies_index, n_branches, branch_index, i_parent, center, radius, expansion_order)
    return MultiBranch(bodies_index, n_branches, branch_index, i_parent, center, radius, initialize_expansion(expansion_order, typeof(radius)), initialize_expansion(expansion_order, typeof(radius)), initialize_harmonics(expansion_order, typeof(radius)), initialize_ML(expansion_order, typeof(radius)), ReentrantLock())
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

@inline function get_child_center(parent_center, parent_radius, i_octant)
    delta = parent_radius / 2.0
    i_octant -= 1
    dx = iseven(i_octant) ? -delta : delta
    i_octant >>= 1
    dy = iseven(i_octant) ? -delta : delta
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
@inline function get_n_children(cumulative_octant_census, n_per_branch)
    n_children = 0
    if get_population(cumulative_octant_census) > n_per_branch
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

# """
# Performs the same sort operation as the tree. (Undoes `unsort!` operation.)
# """
# function resort!(system, tree::SingleTree)
#     buffer = deepcopy(system)
#     index = tree.index
#     for i in 1:length(system)
#         buffer[i] = system[index[i]]
#     end
#     for i in 1:length(system)
#         system[i] = buffer[i]
#     end
# end

# function resort!(systems::Tuple, tree::MultiTree)
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

function unsorted_index_2_sorted_index(i_unsorted, tree::SingleTree)
    return tree.inverse_sort_index[i_unsorted]
end

function unsorted_index_2_sorted_index(i_unsorted, i_system, tree::MultiTree)
    return tree.inverse_sort_index_list[i_system][i_unsorted]
end

function sorted_index_2_unsorted_index(i_sorted, tree::SingleTree)
    return tree.sort_index[i_sorted]
end

function sorted_index_2_unsorted_index(i_unsorted, i_system, tree::MultiTree)
    return tree.sort_index_list[i_system][i_unsorted]
end

#####
##### find the center and radius of a (group of) system(s) of bodies of zero radius
#####
@inline function max_xyz(x_min, y_min, z_min, x_max, y_max, z_max, x, y, z)
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
    
    return x_min, y_min, z_min, x_max, y_max, z_max
end

@inline function max_xyz(x_min, y_min, z_min, x_max, y_max, z_max, system)
    for i in 1:get_n_bodies(system)
        x, y, z = system[i,POSITION]
        x_min, y_min, z_min, x_max, y_max, z_max = max_xyz(x_min, y_min, z_min, x_max, y_max, z_max, x, y, z)
    end

    return x_min, y_min, z_min, x_max, y_max, z_max
end

@inline get_radius(dx,dy,dz,scale_radius) = max(dx,dy,dz) * scale_radius # assume cubic cells
# @inline get_radius(dx,dy,dz,scale_radius) = sqrt(dx*dx + dy*dy + dz*dz) * scale_radius # assume spherical cells

@inline function get_center_radius(x_min, y_min, z_min, x_max, y_max, z_max, scale_radius)
    center = SVector{3}((x_max+x_min)/2, (y_max+y_min)/2, (z_max+z_min)/2)
    radius = get_radius(x_max-center[1], y_max-center[2], z_max-center[3], scale_radius)
    return center, radius
end

function center_radius(systems::Tuple; scale_radius = 1.00001)
    x_min, y_min, z_min = systems[1][1,POSITION]
    x_max, y_max, z_max = systems[1][1,POSITION]
    for system in systems
        x_min, y_min, z_min, x_max, y_max, z_max = max_xyz(x_min, y_min, z_min, x_max, y_max, z_max, system)
    end
    return get_center_radius(x_min, y_min, z_min, x_max, y_max, z_max, scale_radius)
end

function center_radius(system; scale_radius = 1.00001)
    x_min, y_min, z_min = system[1,POSITION]
    x_max, y_max, z_max = system[1,POSITION]
    x_min, y_min, z_min, x_max, y_max, z_max = max_xyz(x_min, y_min, z_min, x_max, y_max, z_max, system)
    return get_center_radius(x_min, y_min, z_min, x_max, y_max, z_max, scale_radius)
end

#####
##### shrinking method for bodies of non-zero radius
#####
@inline function max_xyz_nonzero_radius(x_min, y_min, z_min, x_max, y_max, z_max, x, y, z, radius)
    x_min = min(x_min, x-radius)
    y_min = min(y_min, y-radius)
    z_min = min(z_min, z-radius)
    x_max = max(x_max, x+radius)
    y_max = max(y_max, y+radius)
    z_max = max(z_max, z+radius)
    
    return x_min, y_min, z_min, x_max, y_max, z_max
end

@inline function max_xyz_nonzero_radius(x_min, y_min, z_min, x_max, y_max, z_max, systems::Tuple, bodies_indices)
    for (system, bodies_index) in zip(systems, bodies_indices)
        x_min, y_min, z_min, x_max, y_max, z_max = max_xyz_nonzero_radius(x_min, y_min, z_min, x_max, y_max, z_max, system, bodies_index)
    end
    return x_min, y_min, z_min, x_max, y_max, z_max
end

@inline function max_xyz_nonzero_radius(x_min, y_min, z_min, x_max, y_max, z_max, system, bodies_index::UnitRange)
    for i_body in bodies_index
        x, y, z = system[i_body,POSITION]
        radius = system[i_body,RADIUS]
        x_min, y_min, z_min, x_max, y_max, z_max = max_xyz_nonzero_radius(x_min, y_min, z_min, x_max, y_max, z_max, x, y, z, radius)
    end
    return x_min, y_min, z_min, x_max, y_max, z_max
end

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


@inline function center_nonzero_radius(system, bodies_index)
    x_min, y_min, z_min = first_body_position(system, bodies_index)
    x_max = x_min
    y_max = y_min
    z_max = z_min
    x_min, y_min, z_min, x_max, y_max, z_max = max_xyz_nonzero_radius(x_min, y_min, z_min, x_max, y_max, z_max, system, bodies_index)
    # if get_n_bodies(bodies_index) == 1 # singularity issues in the local expansion if we center the expansion on point where we want to evaluate it
    #     # TODO with a smarter local expansion evaluation, we could get rid of this provision
    #     center = SVector{3}((x_min+x_max)/2.0 + SHRINKING_OFFSET, (y_min+y_max)/2.0, (z_min+z_max)/2.0)
    # else
        center = SVector{3}((x_min+x_max)/2.0, (y_min+y_max)/2.0, (z_min+z_max)/2.0)
    # end
    delta = away_from_center!(center, system, bodies_index)
    
    return center + delta
end

@inline function get_distance_leaf(x, y, z, center)
    dx = x - center[1]
    dy = y - center[2]
    dz = z - center[3]
    return dx, dy, dz #sqrt(dx*dx + dy*dy + dz*dz)
end

@inline function get_distance(x, y, z, center)
    dx = x - center[1]
    dy = y - center[2]
    dz = z - center[3]
    return sqrt(dx*dx + dy*dy + dz*dz)
end

@inline function shrink_radius(radius, center, system, bodies_index)
    for i_body in bodies_index
        x, y, z = system[i_body,POSITION]
        body_radius = system[i_body,RADIUS]
        distance_2_body_center = get_distance(x, y, z, center)
        # if distance_2_body_center < 1e-7
        #     system[i_body,POSITION] .+= 1e-6
        #     distance_2_body_center = get_distance(x, y, z, center)
        # end
        radius = max(radius, distance_2_body_center + body_radius)
    end

    return radius
end

@inline function away_from_center!(center, systems::Tuple, bodies_indices)
    delta = @SVector zeros(3)
    for (system,bodies_index) in zip(systems, bodies_indices)
        delta += away_from_center!(center, system, bodies_index)
    end

    return delta
end

@inline function away_from_center!(center, system, bodies_index)
    for i_body in bodies_index
        x, y, z = system[i_body,POSITION]
        distance_2_body_center = get_distance(x, y, z, center)
        if distance_2_body_center < 1e-7
            return SVector{3}(1e-6,1e-6,1e-6)
        end
    end

    return @SVector zeros(3)
end

@inline function shrink_radius(radius, center, systems::Tuple, bodies_indices)
    for (system, bodies_index) in zip(systems, bodies_indices)
        radius = shrink_radius(radius, center, system, bodies_index)
    end
    return radius
end

@inline function replace_branch!(branch::SubArray{TB,0,<:Any,<:Any,<:Any}, new_center, new_radius) where TB
    # (; bodies_index, n_branches, branch_index, i_parent, center, radius, multipole_expansion, local_expansion, harmonics, ML, lock) = branch[]
    bodies_index = branch[].bodies_index
    n_branches = branch[].n_branches
    branch_index = branch[].branch_index
    i_parent = branch[].i_parent
    center = branch[].center
    radius = branch[].radius
    multipole_expansion = branch[].multipole_expansion
    local_expansion = branch[].local_expansion
    harmonics = branch[].harmonics
    ML = branch[].ML
    lock = branch[].lock
    branch[] = TB(bodies_index, n_branches, branch_index, i_parent, new_center, new_radius, multipole_expansion, local_expansion, harmonics, ML, lock)
end

function shrink_leaf!(branch, system)
    # unpack
    bodies_index = branch[].bodies_index
    
    # recenter # turn this off for now- only shrink/expand radius
    # new_center = center_nonzero_radius(system, bodies_index)
    new_center = branch[].center
    
    # shrink radius 
    new_radius = zero(branch[].radius)
    new_radius = shrink_radius(new_radius, new_center, system, bodies_index)

    # replace branch
    replace_branch!(branch, new_center, new_radius)
end

@inline function max_xyz_nonzero_radius(x_min, y_min, z_min, x_max, y_max, z_max, branch::Branch, child_branches)
    for child_branch in child_branches
        x, y, z = child_branch.center
        radius = child_branch.radius
        x_min, y_min, z_min, x_max, y_max, z_max = max_xyz_nonzero_radius(x_min, y_min, z_min, x_max, y_max, z_max, x, y, z, radius)
    end
    return x_min, y_min, z_min, x_max, y_max, z_max
end

@inline function center_nonzero_radius(branch::Branch, child_branches)
    # bounding rectangle
    x_min, y_min, z_min = branch.center
    x_max = x_min
    y_max = y_min
    z_max = z_min
    x_min, y_min, z_min, x_max, y_max, z_max = max_xyz_nonzero_radius(x_min, y_min, z_min, x_max, y_max, z_max, branch, child_branches)
    
    # find center
    center = SVector{3}((x_min+x_max)/2.0, (y_min+y_max)/2.0, (z_min+z_max)/2.0)
    
    return center
end

@inline function shrink_radius(radius, center, branch::Branch, child_branches)
    for child_branch in child_branches
        x, y, z = child_branch.center
        body_radius = child_branch.radius
        radius = max(radius, get_distance(x, y, z, center) + body_radius)
    end
    return radius
end

function shrink_branch!(branch, child_branches)
    # recenter # turn this off for now
    # new_center = center_nonzero_radius(branch[], child_branches)
    new_center = branch[].center

    # shrink radius
    new_radius = zero(branch[].radius)
    new_radius = shrink_radius(new_radius, new_center, branch[], child_branches)

    # replace branch
    replace_branch!(branch, new_center, new_radius)
end

function shrink_recenter!(branches, levels_index, system)
    for level_index in reverse(levels_index) # start at the bottom level
        for i_branch in level_index
            branch = view(branches,i_branch)
            if branch[].n_branches == 0 # leaf
                shrink_leaf!(branch, system)
            else
                children = view(branches, branch[].branch_index)
                shrink_branch!(branch, children)
            end
        end
    end
end

#####
##### helper function
#####
function initialize_expansion(expansion_order, type=Float64)
    return zeros(type, 2, 4, ((expansion_order+1) * (expansion_order+2)) >> 1)
end

function initialize_harmonics(expansion_order, type=Float64)
    root_n = expansion_order<<1 + 1
    return zeros(type, 2, 2 * root_n * root_n)
end

function initialize_ML(expansion_order, type=Float64)
    # return MArray{Tuple{2,4}, type}(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)
    return zeros(type,2,4)
end

function reset_expansions!(tree)
    T = eltype(tree.branches[1].multipole_expansion)
    for branch in tree.branches
        branch.multipole_expansion .= zero(T)
        branch.local_expansion .= zero(T)
    end
end
