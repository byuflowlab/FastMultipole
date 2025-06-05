#------- matrix storage -------#

function Matrices(sizes::Vector{Tuple{Int,Int}}, TF=Float64)
    # preallocate matrix storage
    n_matrix = sum(m * n for (m, n) in sizes)
    n_rhs = sum(m for (m,_) in sizes)
    data = Vector{TF}(undef, n_matrix)
    rhs = Vector{TF}(undef, n_rhs)
    
    # offsets
    matrix_offsets = Vector{Int}(undef, length(sizes))
    rhs_offsets = Vector{Int}(undef, length(sizes))
    matrix_offset = 1
    rhs_offset = 1
    for i in 1:length(sizes)
        matrix_offsets[i] = matrix_offset
        rhs_offsets[i] = rhs_offset

        m, n = sizes[i]
        matrix_offset += m * n
        rhs_offset += m
    end

    return Matrices{TF}(data, rhs, sizes, matrix_offsets, rhs_offsets)
end

function EmptyMatrices(TF=Float64)
    # preallocate empty matrix storage
    return Matrices{TF}(zeros(TF, 0), zeros(TF, 0), Tuple{Int,Int}[], Int[], Int[])
end

function get_matrix_vector(ms::Matrices, k::Int)
    m, n = ms.sizes[k]
    matrix_offset = ms.matrix_offsets[k]
    rhs_offset = ms.rhs_offsets[k]
    mat = @view ms.data[matrix_offset:matrix_offset + m*n - 1]
    return reshape(mat, m, n), view(ms.rhs, rhs_offset:rhs_offset + m - 1)
end

function set_unit_strength!(source_buffers::Tuple, source_systems::Tuple)
    for i_source_system in eachindex(source_systems)
        source_system = source_systems[i_source_system]
        source_buffer = source_buffers[i_source_system]
        
        set_unit_strength!(source_buffer, source_system)
    end
end

function set_unit_strength!(source_buffer::AbstractMatrix{TF}, source_system) where TF
    # set the source strength to unit
    n_bodies = get_n_bodies(source_system)
    unit_strength = one(TF)
    dim = strength_dims(source_system)
    for j in 1:n_bodies
        value_to_strength!(source_buffer, source_system, j, one(TF))
    end
end

function save_strengths(source_buffers::Tuple, source_systems::Tuple)
    # save the strengths of the source systems
    TF = promote_type(eltype.(source_buffers)...)
    old_strengths = Tuple(Matrix{TF}(undef, strength_dims(source_systems[k]), get_n_bodies(source_systems[k])) for k in eachindex(source_systems))
    for k in eachindex(source_systems)
        source_system = source_systems[k]
        source_buffer = source_buffers[k]
        old_strengths[k] .= view(source_buffer, 4:4 + strength_dims(source_system) - 1, :)
    end

    return old_strengths
end

function restore_strengths!(source_buffers::Tuple, source_systems::Tuple, old_strengths)
    # restore the strengths of the source systems
    for k in eachindex(source_systems)
        source_system = source_systems[k]
        source_buffer = source_buffers[k]
        source_buffer[4:4 + strength_dims(source_system) - 1, :] .= old_strengths[k]
    end
end

function nonself_influence_matrices(target_buffers::Tuple, source_buffers::Tuple, source_systems::Tuple, target_tree::Tree{TF,<:Any}, source_tree::Tree, direct_list, derivatives_switches) where TF
    
    #--- sort by source ---#

    source_branches = source_tree.branches
    target_branches = target_tree.branches
    sorted_list = sort_by_source(direct_list, source_branches)

    #--- pre-allocate influence matrices ---#
    
    if length(sorted_list) > 0
        
        # preallocate sizes
        sizes = Vector{Tuple{Int,Int}}(undef, length(source_tree.leaf_index))

        # number of non-empty matrices
        this_source = 0
        n_matrices = 0
        for (i_target,j_source) in sorted_list
            if j_source != this_source
                n_matrices += 1
                this_source = j_source
            end
        end

        # generate matrix map
        matrix_map = Vector{Int}(undef, n_matrices)
        i_leaf = 1
        i_matrix = 1
        this_source = 0
        for (_, j_source) in sorted_list

            # found a non-empty matrix
            if j_source != this_source
                # upper bound of number of leaves we've skipped
                n_missing_ub = j_source - source_tree.leaf_index[i_leaf]

                # increment i_leaf until we line up with j_source
                for i_plus in 1:n_missing_ub
                    # i_leaf is empty
                    sizes[i_leaf] = (0,0)

                    # move on to the next leaf
                    i_leaf += 1

                    # check if we've reached j_source
                    if source_tree.leaf_index[i_leaf] == j_source
                        break
                    end
                end

                # update matrix_map
                matrix_map[i_matrix] = i_leaf

                # recurse
                i_leaf += 1
                i_matrix += 1
                this_source = j_source
            end
        end

        # fill in the rest with zeros
        for i in i_leaf:length(source_tree.leaf_index)
            sizes[i] = (0,0)
        end
        
        # populate sizes
        this_leaf = source_tree.leaf_index[1]
        this_source = sorted_list[1][2]
        i_matrix = 1
        n = 0
        for (i_target, j_source) in sorted_list
            
            # just finished a matrix
            if j_source != this_source
                
                # save this size
                m = get_n_bodies(source_branches[this_source].bodies_index)
                sizes[matrix_map[i_matrix]] = (n,m)

                # increment indices for the next
                i_matrix += 1
                this_source = j_source
                n = 0
            end
            
            # accumulate n
            n += get_n_bodies(target_branches[i_target].bodies_index)
        end

        # add the final matrix
        this_source = sorted_list[end][2]
        m = get_n_bodies(source_branches[this_source].bodies_index)
        sizes[matrix_map[i_matrix]] = (n,m)

        # construct influence matrices
        matrices = Matrices(sizes, TF)

        #--- populate influence matrices ---#

        # store strengths for later
        old_strengths = save_strengths(source_buffers, source_systems)

        # set sources to unit strength
        set_unit_strength!(source_buffers, source_systems)

        # loop over direct list
        this_source = 0
        i_matrix = 0
        i_target_start = 1
        i_source_start = 1
        matrix, influence = get_matrix_vector(matrices, 1)

        for (i_target,j_source) in sorted_list
                
            # start on the next matrix
            if j_source != this_source
                i_matrix += 1
                this_source = j_source
                i_source_start = 1
                i_target_start = 1
                matrix, influence = get_matrix_vector(matrices, matrix_map[i_matrix])
            end

            # loop over source systems
            this_i_source_start = i_source_start

            for i_source_system in eachindex(source_systems)
                
                # get view of matrix corresponding to this source system
                source_index = source_branches[j_source].bodies_index[i_source_system]
                n_sources = length(source_index)

                # unpack source buffer
                source_system = source_systems[i_source_system]
                source_buffer = source_buffers[i_source_system]
                this_source_buffer = view(source_buffer, :, source_index)

                # loop over target systems
                for i_target_system in eachindex(target_buffers)
                    
                    # get view of matrix corresponding to this source and target system
                    target_index = target_branches[i_target].bodies_index[i_target_system]
                    n_targets = length(target_index)
                    this_matrix = view(matrix, i_target_start:i_target_start + n_targets - 1, this_i_source_start:this_i_source_start + n_sources - 1)
                    this_influence = view(influence, i_target_start:i_target_start + n_targets - 1)

                    # unpack target buffer
                    target_buffer = target_buffers[i_target_system]
                    this_target_buffer = view(target_buffer, :, target_index)
                    derivatives_switch = derivatives_switches[i_target_system]

                    # loop over source bodies
                    for (isb,i_source_body) in enumerate(source_index)
                    
                        # reset targets
                        reset!(target_buffer, target_index)

                        # direct influence on these targets
                        direct!(target_buffer, target_index, derivatives_switch, source_system, source_buffer, i_source_body:i_source_body)

                        # compute influences
                        influence!(this_influence, this_target_buffer, source_system, this_source_buffer)

                        # update matrix
                        this_matrix[:,isb] .= this_influence

                    end

                    # update target starting index
                    i_target_start += n_targets
                end

                # update source starting index
                this_i_source_start += length(source_index)
            end
        end

        # restore old strengths
        restore_strengths!(source_buffers, source_systems, old_strengths)

    else
        matrices = EmptyMatrices(TF)
    end

    return matrices, sorted_list
end

function index_by_source(sorted_list::Vector{SVector{2,Int}}, leaf_index::Vector{Int})
    
    # preallocate index map
    index_map = Vector{UnitRange{Int}}(undef, length(leaf_index))

    # loop over direct list
    this_source = 0
    i_start = 1
    i_leaf = 1
    this_source = sorted_list[1][2] # start with the first source
    for (i_list, (i_target, j_source)) in enumerate(sorted_list)

        # finished a non-empty matrix
        if j_source != this_source

            # upper bound of number of leaves we've skipped
            n_missing_ub = this_source - leaf_index[i_leaf]

            # increment i_leaf until we line up with j_source
            for i_plus in 1:n_missing_ub

                # no interactions for this leaf
                index_map[i_leaf] = i_start:i_start - 1

                # move on to the next leaf
                i_leaf += 1

                # check if we've reached j_source
                if leaf_index[i_leaf] == this_source
                    break
                end
            end

            # store leaf index range
            index_map[i_leaf] = i_start:i_list-1

            # recurse
            i_leaf += 1
            this_source = j_source
            i_start = i_list
        end
    end

    #--- last index ---#

    # upper bound of number of leaves we've skipped
    n_missing_ub = this_source - leaf_index[i_leaf]

    # increment i_leaf until we line up with j_source
    for i_plus in 1:n_missing_ub

        # no interactions for this leaf
        index_map[i_leaf] = i_start:i_start - 1

        # move on to the next leaf
        i_leaf += 1

        # check if we've reached j_source
        if leaf_index[i_leaf] == this_source
            break
        end
    end

    # store leaf index range
    index_map[i_leaf] = i_start:length(sorted_list)
    i_start = length(sorted_list) + 1

    #--- all remaining leaves do not act as sources ---#

    for i in i_leaf + 1:length(leaf_index)
        index_map[i] = i_start:i_start - 1
    end

    return index_map
end

function self_influence_matrices(target_buffers, source_buffers, source_systems, target_tree::Tree{TF,<:Any}, source_tree, derivatives_switches) where TF
    
    #--- pre-allocate influence matrices ---#

    # get sizes
    sizes = Vector{Tuple{Int,Int}}(undef,length(source_tree.leaf_index))
    for (i, i_branch) in enumerate(source_tree.leaf_index)
        n_sources = get_n_bodies(source_tree.branches[i_branch].bodies_index)
        sizes[i] = (n_sources, n_sources)
    end
    
    # construct influence matrices
    matrices = Matrices(sizes, TF)

    #--- populate influence matrices ---#

    # store strengths for later
    old_strengths = save_strengths(source_buffers, source_systems)

    # set sources to unit strength
    set_unit_strength!(source_buffers, source_systems)

    # populate influence matrices
    for (i_matrix, i_leaf) in enumerate(source_tree.leaf_index)
        # get branch
        branch = source_tree.branches[i_leaf]

        # get view of matrix corresponding to this branch
        matrix, influence = get_matrix_vector(matrices, i_matrix)

        # loop over source systems
        i_source_start = 0
        for i_source_system in eachindex(source_systems)

            # unpack source system and buffer
            source_system = source_systems[i_source_system]
            source_buffer = source_buffers[i_source_system]
            source_bodies_index = branch.bodies_index[i_source_system]

            # loop over target systems
            i_target_start = 1
            for i_target_system in eachindex(target_buffers)

                # unpack target buffer
                target_buffer = target_buffers[i_target_system]
                derivatives_switch = derivatives_switches[i_target_system]
                target_bodies_index = branch.bodies_index[i_target_system]

                # view of influence corresponding to these targets
                n_targets = length(target_bodies_index)
                this_influence = view(influence, i_target_start:i_target_start + n_targets - 1)

                # loop over bodies in this branch
                for (isb, i_source_body) in enumerate(source_bodies_index)

                    # reset targets
                    reset!(target_buffer, target_bodies_index)

                    # direct influence on these targets
                    direct!(target_buffer, target_bodies_index, derivatives_switch, source_system, source_buffer, i_source_body:i_source_body)

                    # compute influences
                    influence!(this_influence, view(target_buffer, :, target_bodies_index), source_system, view(source_buffer, :, source_bodies_index))

                    # update matrix
                    matrix[i_target_start:i_target_start + n_targets - 1, isb + i_source_start] .= this_influence

                end

                # update target starting index
                i_target_start += length(branch.bodies_index[i_target_system])
            end
            
            # update source starting index
            i_source_start += length(source_bodies_index)
        end
    end
    
    # restore old strengths
    restore_strengths!(source_buffers, source_systems, old_strengths)

    return matrices
end

function FastGaussSeidel(target_systems::Tuple, source_systems::Tuple; 
    expansion_order=4, multipole_threshold=0.5, leaf_size=30, lamb_helmholtz=true,
    interaction_list_method=Barba(), shrink_recenter=true,
    derivatives_switches=DerivativesSwitch(true, true, false, target_systems)
)

    #--- identical source and target trees ---#

    @assert target_systems === source_systems "different sources and targets are not yet supported for FastGaussSeidel"
    @assert interaction_list_method == Barba() "only the Barba() interaction_list_method is currently supported for FastGaussSeidel"

    #--- generate octree ---#

    # promote leaf_size to vector
    leaf_size = to_vector(leaf_size, length(source_systems))

    # create trees
    TF = promote_type(eltype.(target_systems)...)
    target_tree = Tree(target_systems, true; expansion_order, leaf_size, shrink_recenter)
    source_tree = Tree(source_systems, false; expansion_order, leaf_size, shrink_recenter)

    #--- build interaction lists ---#

    farfield, nearfield, self_induced = true, true, false # self-induced interactions accounted for in self-influence matrices
    m2l_list, direct_list = build_interaction_lists(target_tree.branches, source_tree.branches, leaf_size, multipole_threshold, farfield, nearfield, self_induced, interaction_list_method)

    #--- build non-self influence matrices ---#

    nonself_matrices, sorted_list = nonself_influence_matrices(target_tree.buffers, source_tree.buffers, source_systems, target_tree, source_tree, direct_list, derivatives_switches)
    
    #--- index by source ---#

    index_map = index_by_source(sorted_list, source_tree.leaf_index)

    #--- build self-influence matrices ---#

    self_matrices = self_influence_matrices(target_tree.buffers, source_tree.buffers, source_systems, target_tree, source_tree, derivatives_switches)

    #--- source strength vector ---#

    strengths = zeros(TF, get_n_bodies(source_systems))
    strengths_by_leaf = Vector{UnitRange{Int}}(undef, length(source_tree.leaf_index))
    i_start = 1
    for (i_leaf, i_branch) in enumerate(source_tree.leaf_index)

        # get bodies index
        bodies_index = source_tree.branches[i_branch].bodies_index
        
        # get number of bodies in this branch
        n_bodies = get_n_bodies(bodies_index)
        
        # store strength index
        strengths_by_leaf[i_leaf] = i_start:i_start + n_bodies - 1
        i_start += n_bodies

    end
    
    #--- right-hand side vector ---#
    
    right_hand_side = Vector{TF}(undef, get_n_bodies(target_systems))
    
    #--- external right-hand side vector ---#
    
    external_right_hand_side = Vector{TF}(undef, get_n_bodies(target_systems))
    
    #--- influences per system ---#
    
    influences_per_system = Vector{Vector{TF}}(undef, length(target_systems))
    for i_target_system in eachindex(target_systems)
        # get number of bodies in this target system
        n_bodies = get_n_bodies(target_systems[i_target_system])
        
        # create vector for influences
        influences_per_system[i_target_system] = zeros(TF, n_bodies)
    end

    return FastGaussSeidel{TF,length(source_systems)}(
        self_matrices,
        nonself_matrices,
        index_map,
        m2l_list,
        sorted_list,
        multipole_threshold,
        lamb_helmholtz,
        strengths,
        strengths_by_leaf,
        source_tree,
        target_tree,
        right_hand_side,
        external_right_hand_side,
        influences_per_system
    )
end

function reset!(v::Vector{TF}) where TF
    # reset vector to zero
    v .= zero(TF)
end

function reset!(influences_per_system::Vector{<:AbstractVector})
    for i in eachindex(influences_per_system)
        reset!(influences_per_system[i])
    end
end

function update_strengths_by_leaf!(strengths::Vector{TF}, strengths_by_leaf::Vector{UnitRange{Int}}, source_systems::Tuple, source_buffers::Tuple{<:Matrix}, source_tree::Tree) where TF
    # update strengths by leaf
    for (i_leaf, i_branch) in enumerate(source_tree.leaf_index)
        # get bodies index
        bodies_index = source_tree.branches[i_branch].bodies_index

        # get strengths for this leaf
        strength_index = strengths_by_leaf[i_leaf]
        strengths = view(strengths, strength_index)

        # loop over source systems
        i_strength = 1
        for i_source_system in eachindex(source_systems)
            # unpack source system and buffer
            source_system = source_systems[i_source_system]
            source_buffer = source_buffers[i_source_system]

            # update strengths
            for i_body in bodies_index[i_source_system]
                # set strength value
                strengths[i_strength] = strength_to_value(source_buffer, source_system, i_body)

                # increment index
                i_strength += 1
            end
        end
        @assert i_strength - 1 == get_n_bodies(bodies_index) "number of strengths does not match number of bodies in leaf $(i_leaf)"
    end
end

function solve!(target_systems::Tuple, source_systems::Tuple, solver::FastGaussSeidel{TF,N}; 
    derivatives_switches=DerivativesSwitch(true, true, false, target_systems),
    max_iterations=10, tolerance=1e-3,
) where {TF,N}

    #--- unpack containers ---#

    source_tree = solver.source_tree
    target_tree = solver.target_tree
    source_buffers = source_tree.buffers
    target_buffers = target_tree.buffers
    self_matrices = solver.self_matrices
    nonself_matrices = solver.nonself_matrices
    index_map = solver.index_map
    m2l_list = solver.m2l_list
    direct_list = solver.direct_list
    multipole_threshold = solver.multipole_threshold
    lamb_helmholtz = solver.lamb_helmholtz
    strengths = solver.strengths
    strengths_by_leaf = solver.strengths_by_leaf
    influences_per_system = solver.influences_per_system
    right_hand_side = self_matrices.rhs
    external_right_hand_side = solver.external_right_hand_side

    #--- external right-hand side based on current influence ---#

    # reset and update buffers
    target_influence_to_buffer!(target_buffers, target_systems, derivatives_switches, target_tree.sort_index_list)

    # run influence function on buffers
    reset!(external_right_hand_side)
    influence!(external_right_hand_side, influences_per_system, target_buffers, source_systems, source_buffers, source_tree)

    #--- right-hand side based on non-self influence (current guess) ---#
    
    for i_leaf in source_tree.leaf_index
        # unpack influence matrix and right-hand side
        mat, rhs = get_matrix_vector(nonself_matrices, i_leaf)

        if length(rhs) > 0
            # unpack strengths
            leaf_strengths = view(strengths, strengths_by_leaf[i_leaf])

            # compute the influence
            mul!(rhs, mat, leaf_strengths)

            # move to right-hand side

        end
    end

    #--- update strengths ---#

    update_strengths_by_leaf!(strengths, strengths_by_leaf, source_systems, source_buffers, source_tree)

    #--- fast gauss seidel iterations ---#

    # prepare inputs
    empty_direct_list = Vector{Tuple{Int,Int}}(undef, 0)

    # begin iterations
    for iteration in 1:max_iterations
        
        #--- external influence on the right-hand side ---#

        right_hand_side .= external_right_hand_side

        #--- farfield influence ---#
        
        # fmm call
        reset!(target_buffers)
        fmm!(target_systems, target_tree, source_systems, source_tree, leaf_size_source, m2l_list, empty_direct_list, derivatives_switches, interaction_list_method;
            expansion_order, ε_tol=nothing, lamb_helmholtz,
            upward_pass=true, horizontal_pass=true, downward_pass=true,
            # horizontal_pass_verbose::Bool=false,
            reset_target_tree=true, reset_source_tree=true,
            # nearfield_device::Bool=false,
            tune=false, update_target_systems=false, multipole_threshold,
            # t_source_tree=0.0, t_target_tree=0.0, t_lists=0.0,
            # silence_warnings=false,
        )

        # move to right-hand side
        influence!(right_hand_side, influences_per_system, target_buffers, source_systems, source_buffers, source_tree)

        #--- self influence (solve for strengths) ---#

        for (i_leaf, i_branch) in enumerate(source_tree.leaf_index)
            
            # unpack influence matrix and right-hand side
            mat, rhs = get_matrix_vector(self_matrices, i_leaf)

            # unpack strengths
            leaf_strengths = view(strengths, strengths_by_leaf[i_leaf])

            # solve for strengths
            leaf_strengths .= mat \ rhs

            #--- update non-self influence ---#

            # unpack non-self influence matrix and target influence
            mat, inf = get_matrix_vector(nonself_matrices, i_leaf)

            # remove old influence from right-hand side

            # compute the updated influence
            mul!(inf, mat, leaf_strengths)

            # move to right-hand side

        end

    end

end

"""
    influence!(sorted_influences, influences_per_system, target_buffers, source_systems, source_buffers, source_tree)

Evaluate the influence as pertains to the boundary element influence matrix and subtracts it from `sorted_influences` (which would act like the RHS of a linear system). Based on the current state of the `target_buffers` and `source_buffers`. Note that `source_systems` is provided solely for dispatch. Note also that `influences_per_system` is overwritten each time.

* `sorted_influences::Vector{Float64}`: single vector containing the influence for every body in the target buffers, sorted by source branch in the direct interaction list
* `influences_per_system::Vector{Vector{Float64}}`: vector of vectors containing the influence for each target system, sorted the same way as the buffers
* `target_buffers::NTuple{N,Matrix{Float64}}`: target buffers used to compute the influence
* `source_systems::NTuple{N,<:{UserDefinedSystem}}`: system objects used for dispatch
* `source_buffers::NTuple{N,Matrix{Float64}}`: source buffers used to compute the influence

"""
function influence!(sorted_influences::Vector{TF}, influences_per_system::Vector{Vector{TF}}, target_buffers::Tuple, source_systems::Tuple, source_buffers::Tuple, source_tree::Tree) where TF
    @assert length(target_buffers) == length(source_buffers) == length(source_systems)

    #--- evaluate influences ---#

    for i_system in eachindex(target_buffers)
        # unpack containers
        influence = influences_per_system[i_system]
        target_buffer = target_buffers[i_system]
        source_buffer = source_buffers[i_system]
        source_system = source_systems[i_system]

        # evaluate influence
        influence!(influence, target_buffer, source_system, source_buffer)
    end

    #--- sort by source ---#

    i_influence = 1
    for i_leaf in source_tree.leaf_index
        # unpack bodies index
        bodies_index = source_tree.branches[i_leaf].bodies_index

        for i_system in eachindex(source_systems)
            # unpack containers
            index = bodies_index[i_system]
            influences = influences_per_system[i_system]

            # update influences
            sorted_influences[i_influence:i_influence + length(index) - 1] .-= view(influences, index)
            i_influence += length(index)
        end
    end
end

"""
    influence!(influence, target_buffer, source_system, source_buffer)

Evaluate the influence as pertains to the boundary element influence matrix and overwrites it to `influence` (which would need to be subtracted for it to act like the RHS of a linear system). Based on the current state of the `target_buffer` and `source_buffer`. Should be overloaded for each system type that is used in the boundary element solver.

Note that `source_system` is provided solely for dispatch. It's member bodies will be out of order and should not be referenced.

**Arguments:**

* `influence::AbstractVector{TF}`: vector containing the influence for every body in the target buffer
* `target_buffer::Matrix{TF}`: target buffer used to compute the influence
* `source_system::{UserDefinedSystem}`: system object used solely for dispatch
* `source_buffer::Matrix{TF}`: source buffer used to compute the influence

"""
function influence!(influence, target_buffer, source_system, source_buffer)
    error("influence! not overloaded for systems of type $(typeof(source_system))")
end
