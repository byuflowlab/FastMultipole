#------- tree constructor -------#

function Tree(systems::Tuple, target::Bool, TF=get_type(systems); buffers=allocate_buffers(systems, target, TF), small_buffers = allocate_small_buffers(systems, TF), expansion_order=7, leaf_size=default_leaf_size(systems), n_divisions=20, shrink_recenter=false, allocation_safety_factor=1.0, estimate_cost=false, read_cost_file=false, write_cost_file=false, interaction_list_method=SelfTuning())

    # ensure `systems` isn't empty; otherwise return an empty tree
    if get_n_bodies(systems) > 0

        # determine float type
        TF = Float32
        for system in systems
            TF = promote_type(TF, eltype(system))
        end

        # initialize variables
        octant_container = get_octant_container(systems) # preallocate octant counter; records the total number of bodies in the first octant to start out, but can be reused to count bodies per octant later on
        cumulative_octant_census = get_octant_container(systems) # for creating a cumsum of octant populations
        bodies_index = get_bodies_index(systems)

        # root branch
        center, box = center_box(systems, TF)
        bx, by, bz = box

        # initial octree generation uses cubic cells
        bmax = max(max(bx,by),bz)
        radius = sqrt(bmax*bmax*3.0)
        box = SVector{3}(bmax, bmax, bmax)

        # prepare to divide
        i_first_branch = 2
        sort_index = get_sort_index(systems)
        sort_index_buffer = get_sort_index_buffer(systems)

        # check buffer size
        for (system, buffer) in zip(systems, buffers)
            @assert get_n_bodies(system) == size(buffer, 2) "Buffer doesn't match system size"
            if target
                @assert size(buffer, 1) == 16
            else
                @assert size(buffer, 1) == data_per_body(system)
            end
        end

        # zero buffers
        for buffer in buffers
            buffer .= zero(TF)
        end

        # update buffers with system positions
        target_to_buffer!(buffers, systems)

        # grow root branch
        root_branch, n_children, i_leaf = branch!(buffers, small_buffers, sort_index, octant_container, sort_index_buffer, i_first_branch, bodies_index, center, center, radius, radius, box, box, 0, 1, leaf_size, interaction_list_method) # even though no sorting needed for creating this branch, it will be needed later on; so `branch!` not ony_min creates the root_branch, but also sorts itself into octants and returns the number of children it will have so we can plan array size
        branches = [root_branch] # this first branch will already have its child branches encoded
        # estimated_n_branches = estimate_n_branches(systems, leaf_size, allocation_safety_factor)
        # sizehint!(branches, estimated_n_branches)
        parents_index = 1:1

        # grow branches
        levels_index = [parents_index] # store branches at each level
        for i_divide in 1:n_divisions
            if n_children > 0
                parents_index, n_children, i_leaf = child_branches!(branches, buffers, sort_index, small_buffers, sort_index_buffer, i_leaf, leaf_size, parents_index, cumulative_octant_census, octant_container, n_children, expansion_order, interaction_list_method)
                push!(levels_index, parents_index)
            end
        end

        # check depth
        if n_children > 0
            n_children_prewhile = n_children
            while n_children > 0
                parents_index, n_children, i_leaf = child_branches!(branches, buffers, sort_index, small_buffers, sort_index_buffer, i_leaf, leaf_size, parents_index, cumulative_octant_census, octant_container, n_children, expansion_order, interaction_list_method)
                push!(levels_index, parents_index)
            end
            if WARNING_FLAG_LEAF_SIZE[]
                @warn "leaf_size not reached in for loop, so while loop used to build octree; to improve performance, increase `n_divisions` > $(length(levels_index))"
                WARNING_FLAG_LEAF_SIZE[] = false
            end
        end

        # update buffers with full system data
        if !target
            # old_buffers = deepcopy(buffers)
            system_to_buffer!(buffers, systems, sort_index)
        end

        # invert index
        invert_index!(sort_index_buffer, sort_index)
        inverse_sort_index = sort_index_buffer # reuse the buffer as the inverse index

        # shrink and recenter branches to account for bodies of nonzero radius
        shrink_recenter && shrink_recenter!(branches, levels_index, buffers)

        # store leaves
        leaf_index = Vector{Int}(undef,0)
        sizehint!(leaf_index, length(levels_index[end]))
        for (i_branch,branch) in enumerate(branches)
            if branch.n_branches == 0
                push!(leaf_index, i_branch)
            end
        end

        # allocate expansions
        expansions = initialize_expansions(expansion_order, length(branches), TF)

        # # cost parameters
        # if estimate_cost
        #     params, errors, nearfield_params, nearfield_errors = estimate_tau(systems; read_cost_file=read_cost_file, write_cost_file=write_cost_file)
        #     cost_parameters = CostParameters(params..., nearfield_params)
        # else
        #     cost_parameters = CostParameters(systems)
        # end
        # cost_parameters = Threads.nthreads() > 1 ? direct_cost_estimate(systems, leaf_size) : dummy_direct_cost_estimate(systems, leaf_size)

        # assemble tree
        tree = Tree(branches, expansions, levels_index, leaf_index, sort_index, inverse_sort_index, buffers, small_buffers, expansion_order, leaf_size)#, cost_parameters)

    else
        tree = EmptyTree(systems)
    end

    return tree
end

function Tree(system, target::Bool, TF=eltype(system); optargs...)
    return Tree((system,), target, TF; optargs...)
end

"""
    EmptyTree(system)

Returns an empty tree. Used if `system` is empty.

# Arguments

* `system`: the system from which a tree is to be created

# Returns

* `tree`: if `typeof(system)<:Tuple`, a `::MultiTree` is returned; otherwise, a `::SingleTree` is returned

"""
function EmptyTree(system)
    return EmptyTree((system,))
end

function EmptyTree(system::Tuple)
    TF = get_type(system)
    N = length(system)
    branches = Vector{Branch{TF,N}}(undef,0)
    expansions = Array{TF,4}(undef,0,0,0,0)
    levels_index = Vector{UnitRange{Int64}}(undef,0)
    leaf_index = Int[]
    sort_index_list = Tuple(Int[] for _ in 1:N)
    inverse_sort_index_list = Tuple(Int[] for _ in 1:N)
    buffer = (Matrix{TF}(undef,0,0),)
    small_buffer = [Matrix{TF}(undef,0,0)]
    expansion_order = -1
    leaf_size = SVector{N}(-1 for _ in 1:N)
    return Tree(branches, expansions, levels_index, leaf_index, sort_index_list, inverse_sort_index_list, buffer, small_buffer, expansion_order, leaf_size)
end

#------- construct tree by level rather than leaf size -------#

function getTF(systems::Tuple)
    TF = Float32
    for system in systems
        TF = promote_type(TF, eltype(system))
    end
    return TF
end

"""
Doesn't stop subdividing until ALL child branches have satisfied the leaf size.
"""
function TreeByLevel(systems::Tuple, target::Bool, TF=get_type(systems); centerbox=center_box(systems, getTF(systems)), buffers=allocate_buffers(systems, target, TF), small_buffers = allocate_small_buffers(systems, TF), expansion_order=7, n_levels=5)

    # ensure `systems` isn't empty; otherwise return an empty tree
    if get_n_bodies(systems) > 0

        # determine float type
        TF = Float32
        for system in systems
            TF = promote_type(TF, eltype(system))
        end

        # initialize variables
        octant_container = get_octant_container(systems) # preallocate octant counter; records the total number of bodies in the first octant to start out, but can be reused to count bodies per octant later on
        cumulative_octant_census = get_octant_container(systems) # for creating a cumsum of octant populations
        bodies_index = get_bodies_index(systems)

        # root branch
        center, box = centerbox
        bx, by, bz = box

        # initial octree generation uses cubic cells
        bmax = max(max(bx,by),bz)
        radius = sqrt(bmax*bmax*3.0)
        box = SVector{3}(bmax, bmax, bmax)

        # prepare to divide
        i_first_branch = 2
        sort_index = get_sort_index(systems)
        sort_index_buffer = get_sort_index_buffer(systems)

        # check buffer size
        for (system, buffer) in zip(systems, buffers)
            @assert get_n_bodies(system) == size(buffer, 2) "Buffer doesn't match system size"
            if target
                @assert size(buffer, 1) == 16
            else
                @assert size(buffer, 1) == data_per_body(system)
            end
        end

        # zero buffers
        for buffer in buffers
            buffer .= zero(TF)
        end

        # update buffers with system positions
        target_to_buffer!(buffers, systems)

        # grow root branch
        leaf_size = @SVector zeros(Int, length(systems))
        interaction_list_method = Barba()
        root_branch, n_children, i_leaf = branch!(buffers, small_buffers, sort_index, octant_container, sort_index_buffer, i_first_branch, bodies_index, center, center, radius, radius, box, box, 0, 1, leaf_size, interaction_list_method) # even though no sorting needed for creating this branch, it will be needed later on; so `branch!` not ony_min creates the root_branch, but also sorts itself into octants and returns the number of children it will have so we can plan array size
        branches = [root_branch] # this first branch will already have its child branches encoded
        parents_index = 1:1

        # grow branches
        levels_index = [parents_index] # store branches at each level
        for i_level in 2:n_levels
            last_level = i_level == n_levels
            parents_index, n_children, i_leaf = child_branches_level!(branches, buffers, sort_index, small_buffers, sort_index_buffer, i_leaf, leaf_size, parents_index, cumulative_octant_census, octant_container, n_children, expansion_order, interaction_list_method, last_level)
            push!(levels_index, parents_index)
        end

        # update buffers with full system data
        if !target
            # old_buffers = deepcopy(buffers)
            system_to_buffer!(buffers, systems, sort_index)
        end

        # invert index
        invert_index!(sort_index_buffer, sort_index)
        inverse_sort_index = sort_index_buffer # reuse the buffer as the inverse index

        # shrink and recenter branches to account for bodies of nonzero radius
        shrink_recenter = false
        shrink_recenter && shrink_recenter!(branches, levels_index, buffers)

        # store leaves
        leaf_index = collect(levels_index[end])

        # allocate expansions
        expansions = initialize_expansions(expansion_order, length(branches), TF)

        # # cost parameters
        # if estimate_cost
        #     params, errors, nearfield_params, nearfield_errors = estimate_tau(systems; read_cost_file=read_cost_file, write_cost_file=write_cost_file)
        #     cost_parameters = CostParameters(params..., nearfield_params)
        # else
        #     cost_parameters = CostParameters(systems)
        # end
        # cost_parameters = Threads.nthreads() > 1 ? direct_cost_estimate(systems, leaf_size) : dummy_direct_cost_estimate(systems, leaf_size)

        # assemble tree
        tree = Tree(branches, expansions, levels_index, leaf_index, sort_index, inverse_sort_index, buffers, small_buffers, expansion_order, leaf_size)#, cost_parameters)

    else
        tree = EmptyTree(systems)
    end

    return tree
end

#--- buffers ---#

function allocate_target_buffer(TF, system)
    buffer = zeros(TF, 16, get_n_bodies(system))
    return buffer
end

function allocate_source_buffer(TF, system)
    buffer = zeros(TF, data_per_body(system), get_n_bodies(system))
    return buffer
end

function get_type(systems::Tuple)
    TF = eltype(first(systems))
    for system in systems
        TF = promote_type(TF, eltype(system))
    end
    return TF
end

function get_type(target_systems::Tuple, source_systems::Tuple)
    TF = promote_type(get_type(source_systems), get_type(target_systems))
    return TF
end

"""
    allocate_buffers(systems::Tuple, target::Bool)

Allocates buffers for the given systems. If `target` is `true`, it allocates space for position, scalar potential, gradient, and hessian matrix. Otherwise, it allocates enough memory for the user-defined `source_system_to_buffer!`.
"""
function allocate_buffers(systems::Tuple, target::Bool, TF)
    # create buffers
    if target
        buffers = Tuple(allocate_target_buffer(TF, system) for system in systems)
    else
        buffers = Tuple(allocate_source_buffer(TF, system) for system in systems)
    end

    return buffers
end

"""
    allocate_small_buffers(systems::Tuple; target=false)

Allocates small buffers for the given systems. These buffers are used for temporary storage of body positions for octree sorting.
"""
function allocate_small_buffers(systems::Tuple, TF)
    # create buffers
    small_buffers = Vector{Matrix{TF}}(undef, length(systems))
    for i in eachindex(systems)
        small_buffers[i] = zeros(3, get_n_bodies(systems[i]))
    end

    return small_buffers
end

#--- auxilliary functions ---#

@inline default_leaf_size(systems::Tuple) = SVector{length(systems)}(20 for _ in eachindex(systems))

@inline full_leaf_size(systems::Tuple) = SVector{length(systems)}(get_n_bodies(system) for system in systems)

# function estimate_n_branches(system, leaf_size, allocation_safety_factor)
#     n_bodies = get_n_bodies(system)
#     estimated_n_divisions = Int(ceil(log(8,n_bodies/leaf_size)))
#     estimated_n_branches = div(8^(estimated_n_divisions+1) - 1,7)
#     return Int(ceil(estimated_n_branches * allocation_safety_factor))
# end

function child_branches!(branches, buffers, sort_index, small_buffers, sort_index_buffer, i_leaf, leaf_size, parents_index, cumulative_octant_census, octant_container, n_children, expansion_order, interaction_list_method)
    i_first_branch = parents_index[end] + n_children + 1
    for i_parent in parents_index
        parent_branch = branches[i_parent]
        if parent_branch.n_branches > 0
            # radius of the child branches
            child_radius = parent_branch.target_radius * 0.5
            child_box = parent_branch.target_box * 0.5

            # count bodies per octant
            census!(cumulative_octant_census, buffers, parent_branch.bodies_index, parent_branch.target_center) # doesn't need to sort them here; just count them; the alternative is to save census data for EVERY CHILD BRANCH EACH GENERATION; then I save myself some effort at the expense of more memory allocation, as the octant_census would already be available; then again, the allocation might cost more than I save (which is what my intuition suggests)
            update_octant_accumulator!(cumulative_octant_census)

            # number of child branches
            if exceeds(cumulative_octant_census, leaf_size, interaction_list_method)
                for i_octant in 1:8
                    if get_population(cumulative_octant_census, i_octant) > 0
                        bodies_index = get_bodies_index(cumulative_octant_census, parent_branch.bodies_index, i_octant)
                        child_center = get_child_center(parent_branch.target_center, parent_branch.target_box, i_octant)
                        child_branch, n_grandchildren, i_leaf = branch!(buffers, small_buffers, sort_index, octant_container, sort_index_buffer, i_first_branch, bodies_index, child_center, child_center, child_radius, child_radius, child_box, child_box, i_parent, i_leaf, leaf_size, interaction_list_method)
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

function child_branches_level!(branches, buffers, sort_index, small_buffers, sort_index_buffer, i_leaf, leaf_size, parents_index, cumulative_octant_census, octant_container, n_children, expansion_order, interaction_list_method, last_level)
    i_first_branch = parents_index[end] + n_children + 1
    for i_parent in parents_index
        parent_branch = branches[i_parent]
        if parent_branch.n_branches > 0
            # radius of the child branches
            child_radius = parent_branch.target_radius * 0.5
            child_box = parent_branch.target_box * 0.5

            # count bodies per octant
            census!(cumulative_octant_census, buffers, parent_branch.bodies_index, parent_branch.target_center) # doesn't need to sort them here; just count them; the alternative is to save census data for EVERY CHILD BRANCH EACH GENERATION; then I save myself some effort at the expense of more memory allocation, as the octant_census would already be available; then again, the allocation might cost more than I save (which is what my intuition suggests)
            update_octant_accumulator!(cumulative_octant_census)

            # number of child branches
            if exceeds(cumulative_octant_census, leaf_size, interaction_list_method)
                for i_octant in 1:8
                    if get_population(cumulative_octant_census, i_octant) > 0
                        bodies_index = get_bodies_index(cumulative_octant_census, parent_branch.bodies_index, i_octant)
                        child_center = get_child_center(parent_branch.target_center, parent_branch.target_box, i_octant)
                        child_branch, n_grandchildren, i_leaf = branch!(buffers, small_buffers, sort_index, octant_container, sort_index_buffer, i_first_branch, bodies_index, child_center, child_center, child_radius, child_radius, child_box, child_box, i_parent, i_leaf, leaf_size, interaction_list_method)

                        # remove sub branches if this is the last level
                        if last_level
                            n_bodies = child_branch.n_bodies
                            bodies_index = child_branch.bodies_index
                            n_branches = child_branch.n_branches
                            branch_index = child_branch.branch_index
                            i_parent = child_branch.i_parent
                            i_leaf = child_branch.i_leaf
                            source_center = child_branch.source_center
                            target_center = child_branch.target_center
                            source_radius = child_branch.source_radius
                            target_radius = child_branch.target_radius
                            source_box = child_branch.source_box
                            target_box = child_branch.target_box
                            max_influence = child_branch.max_influence

                            branch_index = 1:0
                            child_branch = typeof(child_branch)(n_bodies, bodies_index, n_branches, 
                                branch_index, i_parent, i_leaf, source_center, target_center, 
                                source_radius, target_radius, source_box, target_box, lock, max_influence)
                        end

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

function branch!(buffer, small_buffer, sort_index, octant_container, sort_index_buffer, i_first_branch, bodies_index, source_center, target_center, source_radius, target_radius, source_box, target_box, i_parent, i_leaf, leaf_size, interaction_list_method)
    # count bodies in each octant
    census!(octant_container, buffer, bodies_index, target_center)

    # cumsum
    update_octant_accumulator!(octant_container)

    # number of child branches
    n_children = get_n_children(octant_container, leaf_size, interaction_list_method)

    if n_children > 0
        # get beginning index of sorted bodies
        octant_beginning_index!(octant_container, bodies_index)

        # sort bodies into octants
        sort_bodies!(buffer, small_buffer, sort_index, octant_container, sort_index_buffer, bodies_index, target_center)
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

    return Branch(bodies_index, n_branches, branch_index, i_parent, i_leaf_index, source_center, target_center, source_radius, target_radius, source_box, target_box), n_children, i_leaf
end

function index_in(bodies_indices, checks)
    for (bodies_index, check) in zip(bodies_indices, checks)
        check && length(bodies_index) > 0 && (return true)
    end
    return false
end

# @inline get_body_positions(system, bodies_index::UnitRange) = (system[i,Position()] for i in bodies_index)

@inline function get_octant(position, center)
    octant = 1
    position[1] > center[1] && (octant += 1)
    position[2] > center[2] && (octant += 2)
    position[3] > center[3] && (octant += 4)
    return octant
end

# @inline get_population(cumulative_octant_census::AbstractVector, i_octant) = i_octant == 1 ? cumulative_octant_census[1] : cumulative_octant_census[i_octant] - cumulative_octant_census[i_octant-1]

@inline get_population(cumulative_octant_census::AbstractMatrix, i_octant) = i_octant == 1 ? sum(cumulative_octant_census[:,1]) : sum(cumulative_octant_census[:,i_octant]) - sum(cumulative_octant_census[:,i_octant-1])

# @inline get_population(cumulative_octant_census::AbstractVector) = cumulative_octant_census[end]

@inline get_population(cumulative_octant_census::AbstractMatrix) = sum(cumulative_octant_census[:,end])

# @inline exceeds(cumulative_octant_census::AbstractVector, leaf_size) = cumulative_octant_census[end] > leaf_size

@inline function exceeds(cumulative_octant_census::AbstractMatrix, leaf_size::SVector, ::SelfTuning)
    fraction = 0.0
    n_bodies = 0
    for i_element in 1:size(cumulative_octant_census, 1)
        # fraction += cumulative_octant_census[i_element,end] / leaf_size[i_element]
        n = cumulative_octant_census[i_element,end]
        n_bodies += n
        fraction += n / (leaf_size[i_element] * leaf_size[i_element])
        # cumulative_octant_census[i_element,end] > leaf_size[i_element] && (return true)
    end
    fraction *= minimum(leaf_size)

    return fraction > 0.5 && n_bodies > 1
    # return fraction > nextfloat(1.0)
end

@inline function exceeds(branch::Branch, leaf_size, ::SelfTuning)
    fraction = 0.0
    n_bodies = 0
    for i_element in 1:length(leaf_size)
        n = length(leaf_size.bodies_index[i_element])
        n_bodies += n
        fraction += n / (leaf_size[i_element] * leaf_size[i_element])
    end
    fraction *= minimum(leaf_size)

    return fraction > 0.5 && n_bodies > 1
end

@inline function exceeds(cumulative_octant_census::AbstractMatrix, leaf_size::SVector, ::SelfTuningTreeStop)
    fraction = 0.0
    n_bodies = 0
    for i_element in 1:size(cumulative_octant_census, 1)
        # fraction += cumulative_octant_census[i_element,end] / leaf_size[i_element]
        n = cumulative_octant_census[i_element,end]
        n_bodies += n
        fraction += n / (leaf_size[i_element] * leaf_size[i_element])
        # cumulative_octant_census[i_element,end] > leaf_size[i_element] && (return true)
    end
    fraction *= minimum(leaf_size)

    return fraction > 2.8284271247461903 && n_bodies > 1 # sqrt(8) == 2.8284271247461903; leads to the expectation value
                                                         # of n_source * n_target => leaf_size^2 (break-even point)
    # return fraction > nextfloat(1.0)
end

@inline function exceeds(branch::Branch, leaf_size, ::SelfTuningTreeStop)
    fraction = 0.0
    n_bodies = 0
    for i_element in 1:length(leaf_size)
        n = length(leaf_size.bodies_index[i_element])
        n_bodies += n
        fraction += n / (leaf_size[i_element] * leaf_size[i_element])
    end
    fraction *= minimum(leaf_size)

    return fraction > 2.8284271247461903 && n_bodies > 1
end

@inline function exceeds(cumulative_octant_census::AbstractMatrix, leaf_size, ::Barba)
    not_over = true
    for i_element in 1:size(cumulative_octant_census, 1)
        n = cumulative_octant_census[i_element,end]
        not_over = not_over && n <= leaf_size[i_element]
    end

    return !not_over
end

@inline function exceeds(branch::Branch, leaf_size, ::Barba)
    not_over = true
    for i_element in 1:length(leaf_size)
        n = length(leaf_size.bodies_index[i_element])
        not_over = not_over && n <= leaf_size[i_element]
    end

    return !not_over
end

# @inline function exceeds(cumulative_octant_census::AbstractMatrix, leaf_size)
#     n_bodies = 0
#     for i_element in 1:size(cumulative_octant_census, 1)
#         n = cumulative_octant_census[i_element,end]
#         n_bodies += n
#     end
#
#     return n_bodies > maximum(leaf_size)
# end

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

function get_octant_container(systems)
    census = MMatrix{length(systems),8,Int64}(undef)
    for (i,system) in enumerate(systems)
        census[i,1] = get_n_bodies(system)
    end
    return census
end

# function get_octant_container(system)
#     census = MVector{8,Int64}(undef)
#     census[1] = get_n_bodies(system)
#     return census
# end

function census!(octant_container::AbstractVector, system::Matrix, bodies_index, center)
    octant_container .= zero(eltype(octant_container))
    # for position in get_body_positions(system, bodies_index)
    for i_body in bodies_index
        position = get_position(system, i_body)
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

# @inline update_octant_accumulator!(octant_population::AbstractVector) = cumsum!(octant_population, octant_population)

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
    return bodies_index
end

@inline get_bodies_index(system) = 1:get_n_bodies(system)

@inline function get_bodies_index(systems::Tuple)
    n_systems = length(systems)
    return SVector{n_systems,UnitRange{Int64}}([get_bodies_index(system) for system in systems])
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

#--- determine the number of descendants ---#

@inline function get_n_children(cumulative_octant_census, leaf_size, interaction_list_method)
    n_children = 0

    if exceeds(cumulative_octant_census, leaf_size, interaction_list_method)
        for i_octant in 1:8
            get_population(cumulative_octant_census,i_octant) > 0 && (n_children += 1)
        end
    end
    return n_children
end

#####
##### sort bodies into the octree
#####
# function sort_bodies!(system::SortWrapper, sort_index, octant_indices::AbstractVector, buffer, sort_index_buffer, bodies_index::UnitRange, center)
#     # sort indices
#     for i_body in bodies_index
#         i_octant = get_octant(system[i_body,Position()], center)
#         sort_index_buffer[octant_indices[i_octant]] = sort_index[i_body]
#         octant_indices[i_octant] += 1
#     end
#     # place buffers
#     sort_index[bodies_index] .= view(sort_index_buffer,bodies_index)
# end

function sort_bodies!(buffer::Matrix, small_buffer::Matrix, sort_index, octant_indices::AbstractVector, sort_index_buffer, bodies_index::UnitRange, center)
    # sort indices
    for i_body in bodies_index
        # identify octant
        i_octant = get_octant(get_position(buffer, i_body), center)
        this_i = octant_indices[i_octant]

        # update small buffer
        small_buffer[1:3,this_i] .= view(buffer, 1:3, i_body)
        # tmp = system[i_body, Body()]
        # buffer[this_i] = tmp

        # update sort index
        sort_index_buffer[octant_indices[i_octant]] = sort_index[i_body]

        # increment octant census
        octant_indices[i_octant] += 1
    end

    # place buffers
    for i_body in bodies_index
        buffer[1:3, i_body] .= view(small_buffer, 1:3, i_body)
    end

    sort_index[bodies_index] .= view(sort_index_buffer, bodies_index)
end

function sort_bodies!(buffers, small_buffers, sort_indices, octant_indices::AbstractMatrix, sort_index_buffers, bodies_indices::AbstractVector, center)
    for (i_system, (sort_index, buffer, small_buffer, sort_index_buffer, bodies_index)) in enumerate(zip(sort_indices, buffers, small_buffers, sort_index_buffers, bodies_indices))
        sort_bodies!(buffer, small_buffer, sort_index, view(octant_indices,i_system,:), sort_index_buffer, bodies_index, center)
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
# """
# Undoes the sort operation performed by the tree.
# """
# function unsort!(systems::Tuple, tree::Tree)
#     for (system, buffer, inverse_sort_index) in zip(systems, tree.buffers, tree.inverse_sort_index_list)
#         unsort!(system, buffer, inverse_sort_index)
#     end
# end

# function unsort!(system, tree::SingleTree)
#     unsort!(system, tree.buffer, tree.inverse_sort_index)
# end

# @inline function unsort!(systems, buffer, tree::Tree)
#     unsort!(systems, buffer, tree.inverse_sort_index)
# end
#
# @inline function unsort!(system, buffer, inverse_sort_index)
#     for i_body in 1:get_n_bodies(system)
#         buffer[i_body] = system[inverse_sort_index[i_body], Body()]
#     end
#     for i_body in 1:get_n_bodies(system)
#         system[i_body, Body()] = buffer[i_body]
#     end
# end

# @inline function unsort!(system::SortWrapper, buffer, inverse_sort_index)
#     system.index .= view(system.index,inverse_sort_index)
# end

# """
# Performs the same sort operation as the tree. (Undoes `unsort!` operation.)
# """
# function resort!(systems::Tuple, tree::Tree)
#     for (system, buffer, sort_index) in zip(systems, tree.buffers, tree.sort_index_list)
#         resort!(system, buffer, sort_index)
#     end
# end

# function resort!(system::SortWrapper, buffer, sort_index)
#     system.index .= sort_index
# end

# function resort!(system, tree::SingleTree)
#     resort!(system, tree.buffer, tree.sort_index)
# end

# function resort!(system, buffer, tree::Tree)
#     resort!(system, buffer, tree.sort_index)
# end

# function resort!(system, buffer, sort_index)
#     for i_body in 1:get_n_bodies(system)
#         buffer[i_body] = system[sort_index[i_body], Body()]
#     end
#     for i_body in 1:get_n_bodies(system)
#         system[i_body, Body()] = buffer[i_body]
#     end
# end

# function resort!(system::SortWrapper, buffer, sort_index)
#     system.index .= sort_index
# end

# @inline function unsorted_index_2_sorted_index(i_unsorted, tree::SingleTree)
#     return tree.inverse_sort_index[i_unsorted]
# end

@inline function unsorted_index_2_sorted_index(i_unsorted, i_system, tree::Tree)
    return tree.inverse_sort_index_list[i_system][i_unsorted]
end

# @inline function sorted_index_2_unsorted_index(i_sorted, tree::SingleTree)
#     return tree.sort_index[i_sorted]
# end

@inline function sorted_index_2_unsorted_index(i_unsorted, i_system, tree::Tree)
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
        x, y, z = get_position(system, i)
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

function center_box(systems, TF)
    bodies_indices = get_bodies_index(systems)
    return center_box(systems, bodies_indices, TF)
end

function center_box(systems::Tuple, bodies_indices, TF)
    x_min, y_min, z_min = first_body_position(systems, bodies_indices, TF)
    x_max, y_max, z_max = x_min, y_min, z_min
    for (system, bodies_index) in zip(systems, bodies_indices)
        x_min, x_max, y_min, y_max, z_min, z_max = max_xyz(x_min, x_max, y_min, y_max, z_min, z_max, system, bodies_index)
    end

    return get_center_box(x_min, x_max, y_min, y_max, z_min, z_max)
end

# function center_box(system, bodies_index)
#     x_min, y_min, z_min = first_body_position(system, bodies_index)
#     x_max, y_max, z_max = x_min, y_min, z_min
#     x_min, x_max, y_min, y_max, z_min, z_max = max_xyz(x_min, x_max, y_min, y_max, z_min, z_max, system, bodies_index)
#     return get_center_box(x_min, x_max, y_min, y_max, z_min, z_max)
# end

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
        x, y, z = get_position(system, i)
        r = get_radius(system, i)
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

function source_center_box(systems::Tuple, bodies_indices, TF)
    x_min, y_min, z_min = first_body_position(systems, bodies_indices, TF)
    x_max, y_max, z_max = x_min, y_min, z_min
    for (system, bodies_index) in zip(systems, bodies_indices)
        x_min, x_max, y_min, y_max, z_min, z_max = source_max_xyz(x_min, x_max, y_min, y_max, z_min, z_max, system, bodies_index)
    end

    return get_center_box(x_min, x_max, y_min, y_max, z_min, z_max)
end

# function source_center_box(system, bodies_index)
#     x_min, y_min, z_min = first_body_position(system, bodies_index)
#     x_max, y_max, z_max = x_min, y_min, z_min
#     x_min, x_max, y_min, y_max, z_min, z_max = source_max_xyz(x_min, x_max, y_min, y_max, z_min, z_max, system, bodies_index)
#     return get_center_box(x_min, x_max, y_min, y_max, z_min, z_max)
# end

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

# @inline function first_body_position(system, bodies_index)
#     return system[bodies_index[1],POSITION]
# end

@inline function first_body_position(systems::Tuple, bodies_indices, TF)
    for (system, bodies_index) in zip(systems, bodies_indices)
        length(bodies_index) > 0 && (return get_position(system, bodies_index[1]))
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

# @inline get_n_bodies_vec(branch::SingleBranch) = length(branch.bodies_index)

@inline get_n_bodies_vec(branch::Branch) = SVector{Int}(length(bodies_index) for bodies_index in branch.bodies_index)


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
        x, y, z = get_position(system, i_body)
        body_radius = get_radius(system, i_body)

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

function shrink_leaf!(branches::Vector{Branch{TF,N}}, i_branch, system) where {TF,N}
    # unpack
    branch = branches[i_branch]
    bodies_index = branch.bodies_index

    # recenter and target box
    new_target_center, new_target_box = center_box(system, bodies_index, TF)
    new_source_center, new_source_box = source_center_box(system, bodies_index, TF)

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

# """
#     accumulate_charge!(branches, systems, branch_index)
#
# Computes the sum of the absolute value of body source strengths, as well as the sum of the L2 norm of the dipole strengths, and store them in each branch, inheriting values from child branches if available.
#
# """
# function accumulate_charge!(branches::AbstractVector{<:Branch{TF}}, systems, branch_index) where TF
#
#     # loop over branches
#     for i_branch in branch_index
#         branch = branches[i_branch]
#
#         # inherit from child branches
#         if branch.n_branches > 0
#             charge = zero(TF)
#             dipole = zero(TF)
#             for i_child in branch.branch_index
#                 child = branches[i_child]
#                 charge += child.charge
#                 dipole += child.dipole
#             end
#
#             # replace branch
#             replace_branch!(branches, i_branch, charge, dipole)
#
#         # compute based on bodies directly
#         else
#             accumulate_charge_bodies!(branches, i_branch, systems)
#         end
#     end
# end
#
# """
#     accumulate_charge!(branch, systems)
#
# Computes the sum of the absolute value of body source strengths, as well as the sum of the L2 norm of the dipole strengths, based on source bodies directly, and store them in each branch.
#
# """
# function accumulate_charge_bodies!(branches::Vector{SingleBranch{TF}}, i_branch, system) where TF
#     # reset multipole expansion, just in case
#     branch = branches[i_branch]
#     branch.multipole_expansion .= zero(TF)
#
#     # only need through p=1 for this
#     expansion_order = 1
#
#     # accumulate charge and dipole
#     charge = zero(TF)
#     dipole = zero(TF)
#
#     # loop over bodies
#     for i_body in branch.bodies_index
#
#         # reset required expansions
#         branch.multipole_expansion[1:2,1:2,1:3] .= zero(TF)
#
#         # multipole coefficients
#         body_to_multipole!(system, branch, i_body:i_body, branch.harmonics, expansion_order)
#
#         # get effective charge
#         charge += abs(branch.multipole_expansion[1,1,1])
#
#         # get effective dipole moment
#         qz = branch.multipole_expansion[1,1,2]
#         qy = branch.multipole_expansion[1,1,3] * 2.0
#         qx = branch.multipole_expansion[2,1,3] * 2.0
#         dipole += sqrt(qx*qx + qy*qy + qz*qz)
#     end
#
#     # reset multipole expansion again
#     branch.multipole_expansion .= zero(TF)
#
#     # replace branch
#     replace_branch!(branches, i_branch, charge, dipole)
#
# end
#
# function accumulate_charge_bodies!(branches::Vector{MultiBranch{TF,N}}, i_branch, systems) where {TF,N}
#     # reset multipole expansion, just in case
#     branch = branches[i_branch]
#     branch.multipole_expansion .= zero(TF)
#
#     # only need through p=1 for this
#     expansion_order = 1
#
#     # accumulate charge and dipole
#     charge = zero(TF)
#     dipole = zero(TF)
#
#     # loop over systems
#     for i_system in 1:N
#         system = systems[i_system]
#
#         # loop over bodies
#         for i_body in branch.bodies_index[i_system]
#
#             # reset required expansions
#             branch.multipole_expansion[1:2,1:2,1:3] .= zero(TF)
#
#             # multipole coefficients
#             body_to_multipole!(system, branch, i_body:i_body, branch.harmonics, expansion_order)
#
#             # get effective charge
#             charge += abs(branch.multipole_expansion[1,1,1])
#
#             # get effective dipole moment
#             qz = branch.multipole_expansion[1,1,2]
#             qy = branch.multipole_expansion[1,1,3] * 2.0
#             qx = branch.multipole_expansion[2,1,3] * 2.0
#             dipole += sqrt(qx*qx + qy*qy + qz*qz)
#         end
#     end
#
#     # reset multipole expansion again
#     branch.multipole_expansion .= zero(TF)
#
#     # replace branch
#     replace_branch!(branches, i_branch, charge, dipole)
#
# end

#--- helper function ---#

@inline function replace_branch!(branches::Vector{TB}, i_branch, new_source_center, new_target_center, new_source_radius, new_target_radius, new_source_box, new_target_box) where TB
    branch = branches[i_branch]
    n_bodies = branch.n_bodies
    bodies_index = branch.bodies_index
    n_branches = branch.n_branches
    branch_index = branch.branch_index
    i_parent = branch.i_parent
    i_leaf = branch.i_leaf
    max_influence = branch.max_influence
    branches[i_branch] = TB(n_bodies, bodies_index, n_branches, branch_index, i_parent, i_leaf, new_source_center, new_target_center, new_source_radius, new_target_radius, new_source_box, new_target_box, max_influence)
end

function initialize_expansion(expansion_order, type=Float64)
    # incrememnt expansion order to make room for error predictions
    # expansion_order += 1

    return zeros(type, 2, 2, ((expansion_order+1) * (expansion_order+2)) >> 1)
end

function initialize_expansions(expansion_order, n_branches, type=Float64)
    # incrememnt expansion order to make room for error predictions
    # expansion_order += 1

    return zeros(type, 2, 2, ((expansion_order+1) * (expansion_order+2)) >> 1, n_branches)
end

function initialize_gradient_n_m(expansion_order, type=Float64)
    # incrememnt expansion order to make room for error predictions
    # expansion_order += 1

    p = expansion_order
    n_harmonics = harmonic_index(p,p)
    return zeros(type, 2, 3, n_harmonics)
end

function initialize_harmonics(expansion_order, type=Float64)
    # incrememnt expansion order to make room for error predictions
    # expansion_order += 1

    p = expansion_order+2
    n_harmonics = harmonic_index(p,p)
    return zeros(type, 2, 2, n_harmonics)
end

function reset_expansions!(tree)
    tree.expansions .= zero(eltype(tree.expansions))
end

#--- debugging functions ---#

function update_family_tree!(family_tree, tree, i_branch)
	branch = tree.branches[i_branch]
	if branch.i_parent > 0
		push!(family_tree, branch.i_parent)
		update_family_tree!(family_tree, tree, branch.i_parent)
	end
end

function get_interaction_list(tree, m2l_list, i_target)
	# get family tree
	family_tree = [i_target]
	update_family_tree!(family_tree, tree, i_target)

	# check m2l list
	interaction_list = Int[]
	for (i, j) in m2l_list
		if i in family_tree
			push!(interaction_list, j)
		end
	end

	return interaction_list
end
