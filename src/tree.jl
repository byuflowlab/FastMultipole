const BRANCH_TYPE = Float64
global SHRINKING_OFFSET = .000001

#####
##### tree constructor
#####
function Tree(system; expansion_order=7, n_per_branch=100, ndivisions=5, scale_radius=1.00001, shrink_recenter=false, allocation_safety_factor=1.0)
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
    root_branch, n_children = Branch(system, sort_index, octant_container, buffer, sort_index_buffer, i_first_branch, bodies_index, center, radius, n_per_branch, expansion_order) # even though no sorting needed for creating this branch, it will be needed later on; so `Branch` not ony_min creates the root_branch, but also sorts itself into octants and returns the number of children it will have so we can plan array size
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
        println("n_per_branch not reached; n_children = $n_children")
    end

    # invert index
    invert_index!(sort_index_buffer, sort_index)
    inverse_sort_index = sort_index_buffer # reuse the buffer as the inverse index

    # shrink and recenter branches to account for bodies of nonzero radius
    if shrink_recenter
        shrink_recenter!(branches, levels_index, system)
    end

    # # count leaves
    # n_leaves = 0
    # for branch in branches
    #     branch.n_branches == 0 && (n_leaves += 1)
    # end

    # # create leaf index and leaf count
    # cumulative_count = Vector{Int64}(undef,n_leaves)
    # cumulative_count[1] = 0
    # leaf_index = zeros(Int,n_leaves)
    # i_leaf_index = 1
    # for (i_branch, branch) in enumerate(branches)
    #     if branch.n_branches == 0
    #         leaf_index[i_leaf_index] = Int32(i_branch)
    #         i_leaf_index < n_leaves && (cumulative_count[i_leaf_index+1] = cumulative_count[i_leaf_index] + sum(branch.n_bodies))
    #         i_leaf_index += 1
    #     end
    # end

    # assemble tree
    tree = Tree(branches, levels_index, sort_index, inverse_sort_index, buffer, expansion_order, n_per_branch)

    return tree
end

Tree(branches::Vector{<:SingleBranch}, levels_index, sort_index, inverse_sort_index, buffer, expansion_order, n_per_branch) = 
    SingleTree(branches, levels_index, sort_index, inverse_sort_index, buffer, expansion_order, n_per_branch)

Tree(branches::Vector{<:MultiBranch}, levels_index, sort_index, inverse_sort_index, buffer, expansion_order, n_per_branch) = 
    MultiTree(branches, levels_index, sort_index, inverse_sort_index, buffer, expansion_order, n_per_branch)

@inline total_n_bodies(system) = length(system)

@inline function total_n_bodies(systems::Tuple)
    n_bodies = 0
    for system in systems
        n_bodies += total_n_bodies(system)
    end
    return n_bodies
end

function estimate_n_branches(system, n_per_branch, allocation_safety_factor)
    n_bodies = total_n_bodies(system)
    estimated_n_divisions = Int(ceil(log(8,n_bodies/n_per_branch)))
    estimated_n_branches = div(8^(estimated_n_divisions+1) - 1,7)
    return Int(ceil(estimated_n_branches * allocation_safety_factor))
end

function child_branches!(branches, system, sort_index, buffer, sort_index_buffer, n_per_branch, parents_index, cumulative_octant_census, octant_container, n_children, expansion_order)
    i_first_branch = parents_index[end] + n_children + 1
    for parent_branch in branches[parents_index]
        if parent_branch.n_branches > 0
            # radius of the child branches
            child_radius = parent_branch.radius / 2.0

            # count bodies per octant
            census!(cumulative_octant_census, system, parent_branch.bodies_index, parent_branch.center) # doesn't need to sort them here; just count them; the alternative is to save census data for EVERY CHILD BRANCH EACH GENERATION; then I save myself some effort at the expense of more memory allocation, as the octant_census would already be available; then again, the allocation might cost more than I save (which is what my intuition suggests)
            update_octant_accumulator!(cumulative_octant_census)
            
            # number of child branches
            if get_population(cumulative_octant_census) > n_per_branch
                for i_octant in 1:8
                    if get_population(cumulative_octant_census, i_octant) > 0  
                        bodies_index = get_bodies_index(cumulative_octant_census, parent_branch.bodies_index, i_octant)
                        child_center = get_child_center(parent_branch.center, parent_branch.radius, i_octant)
                        child_branch, n_grandchildren = Branch(system, sort_index, octant_container, buffer, sort_index_buffer, i_first_branch, bodies_index, child_center, child_radius, n_per_branch, expansion_order)
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

function Branch(system, sort_index, octant_container, buffer, sort_index_buffer, i_first_branch, bodies_index, center, radius, n_per_branch, expansion_order)
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

    return Branch(bodies_index, n_branches, branch_index, center, radius, expansion_order), n_children
end

function Branch(bodies_index::UnitRange, n_branches, branch_index, center, radius, expansion_order)
    return SingleBranch(bodies_index, n_branches, branch_index, center, radius, initialize_expansion(expansion_order), initialize_expansion(expansion_order), ReentrantLock())
end

function Branch(bodies_index, n_branches, branch_index, center, radius, expansion_order)
    return MultiBranch(bodies_index, n_branches, branch_index, center, radius, initialize_expansion(expansion_order), initialize_expansion(expansion_order), ReentrantLock())
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
        census[i,1] = length(system)
    end
    return census
end

function get_octant_container(system)
    census = MVector{8,Int64}(undef)
    census[1] = length(system)
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
    for i_octant in 8:-1:2
        cumulative_octant_census[i_octant] = cumulative_octant_census[i_octant-1] + bodies_index[1]
    end
    cumulative_octant_census[1] = bodies_index[1]
    return cumulative_octant_census
end

@inline function octant_beginning_index!(cumulative_octant_census::AbstractMatrix, bodies_indices::AbstractVector)
    for (i_system,bodies_index) in enumerate(bodies_indices)
        octant_beginning_index!(view(cumulative_octant_census,i_system,:), bodies_index)
    end
    return cumulative_octant_census
end

@inline function get_bodies_index(cumulative_octant_census::AbstractVector, parent_bodies_index::UnitRange, i_octant)
    first_offset = i_octant == 1 ? 0 : cumulative_octant_census[i_octant-1]
    bodies_index = parent_bodies_index[1] + first_offset : parent_bodies_index[1] + cumulative_octant_census[i_octant] - 1
    return bodies_index
end

@inline function get_bodies_index(cumulative_octant_census::AbstractMatrix, parent_bodies_indices::AbstractVector, i_octant)
    n_systems = size(cumulative_octant_census,1)
    bodies_index = SVector{n_systems,UnitRange{Int64}}([get_bodies_index(view(cumulative_octant_census,i_system,:), parent_bodies_indices[i_system], i_octant) for i_system in 1:n_systems])
    # TODO why does this require square brackets? ah, because StaticArrays tries to collect the unit range ONLY IF THERE IS ONLY ONE;
    return bodies_index
end

@inline get_bodies_index(system) = 1:length(system)

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
    return [buffer_element(system) for _ in 1:length(system)]
end

@inline function get_buffer(systems::Tuple)
    return Tuple(get_buffer(system) for system in systems)
end

@inline function get_sort_index(system::SortWrapper)
    return system.index
end

@inline function get_sort_index(system)
    return collect(1:length(system))
end

@inline function get_sort_index(systems::Tuple)
    return Tuple(get_sort_index(system) for system in systems)
end

@inline function get_sort_index_buffer(system) # need not be overloaded for SortWrapper as it will be the same
    return Vector{Int64}(undef,length(system))
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
function sort_bodies!(system, sort_index, octant_indices::AbstractVector, buffer, sort_index_buffer, bodies_index::UnitRange, center)
    # sort indices
    for i_body in bodies_index
        i_octant = get_octant(system[i_body,POSITION], center)
        buffer[octant_indices[i_octant]] = system[i_body]
        sort_index_buffer[octant_indices[i_octant]] = sort_index[i_body]
        octant_indices[i_octant] += 1
    end
    # place buffers
    for i_body in bodies_index
        system[i_body] = buffer[i_body]
    end
    sort_index[bodies_index] .= sort_index_buffer[bodies_index]
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
    for i_body in 1:length(system)
        buffer[i_body] = system[inverse_sort_index[i_body]]
    end
    for i_body in 1:length(system)
        system[i_body] = buffer[i_body]
    end
end

@inline function unsort!(system::SortWrapper, buffer, inverse_sort_index)
    system.index .= system.index[inverse_sort_index]
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

# function get_sorted_body(vanilla_system, tree::SingleTree, i_sorted)
#     return vanilla_system[tree.index[i_sorted]]
# end

# function get_sorted_body(vanilla_system, tree::MultiTree, i_system, i_sorted)
#     return vanilla_system[tree.index[i_system][i_sorted]]
# end

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
    for i in 1:length(system)
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
    return systems[1][bodies_indices[1][1],POSITION]
end

@inline get_n_bodies(bodies_index::UnitRange) = length(bodies_index)

@inline function get_n_bodies(bodies_indices)
    n_bodies = 0
    for bodies_index in bodies_indices
        n_bodies += get_n_bodies(bodies_index)
    end
    return n_bodies
end

@inline function center_nonzero_radius(system, bodies_index)
    x_min, y_min, z_min = first_body_position(system, bodies_index)
    x_max = x_min
    y_max = y_min
    z_max = z_min
    x_min, y_min, z_min, x_max, y_max, z_max = max_xyz_nonzero_radius(x_min, y_min, z_min, x_max, y_max, z_max, system, bodies_index)
    if get_n_bodies(bodies_index) == 1 # singularity issues in the local expansion if we center the expansion on point where we want to evaluate it
        # TODO with a smarter local expansion evaluation, we could get rid of this provision
        center = SVector{3}((x_min+x_max)/2.0 + SHRINKING_OFFSET, (y_min+y_max)/2.0, (z_min+z_max)/2.0)
    else
        center = SVector{3}((x_min+x_max)/2.0, (y_min+y_max)/2.0, (z_min+z_max)/2.0)
    end
    
    return center
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
        radius = max(radius, get_distance(x, y, z, center) + body_radius)
    end
    return radius
end

@inline function shrink_radius(radius, center, systems::Tuple, bodies_indices)
    for (system, bodies_index) in zip(systems, bodies_indices)
        radius = shrink_radius(radius, center, system, bodies_index)
    end
    return radius
end

@inline function replace_branch!(branch::SubArray{TB,0,<:Any,<:Any,<:Any}, new_center, new_radius) where TB
    (; bodies_index, n_branches, branch_index, center, radius, multipole_expansion, local_expansion, lock) = branch[]
    branch[] = TB(bodies_index, n_branches, branch_index, new_center, new_radius, multipole_expansion, local_expansion, lock)
end

function shrink_leaf!(branch, system)
    # unpack
    bodies_index = branch[].bodies_index
    
    # recenter
    new_center = center_nonzero_radius(system, bodies_index)
    
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
    # recenter
    new_center = center_nonzero_radius(branch[], child_branches)

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
function initialize_expansion(expansion_order)
    return zeros(Complex{Float64}, 4, ((expansion_order+1) * (expansion_order+2)) >> 1)
end

function reset_expansions!(tree)
    T = eltype(tree.branches[1].multipole_expansion)
    for branch in tree.branches
        branch.multipole_expansion .= zero(T)
        branch.local_expansion .= zero(T)
    end
end
