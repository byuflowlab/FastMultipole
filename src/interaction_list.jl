#--- create interaction lists ---#

function build_interaction_lists(target_branches, source_branches, source_leaf_size, multipole_threshold, farfield, nearfield, self_induced)

    # prepare containers
    m2l_list = Vector{SVector{2,Int32}}(undef,0)
    direct_list = Vector{SVector{2,Int32}}(undef,0)

    # populate lists
    build_interaction_lists!(m2l_list, direct_list, Int32(1), Int32(1), target_branches, source_branches, source_leaf_size, multipole_threshold, Val(farfield), Val(nearfield), Val(self_induced))

    return m2l_list, direct_list
end

mean(x) = sum(x) / length(x)

function build_interaction_lists!(m2l_list, direct_list, i_target, j_source, target_branches, source_branches, source_leaf_size, multipole_threshold, farfield::Val{ff}, nearfield::Val{nf}, self_induced::Val{si}) where {ff,nf,si}
    # unpack
    source_branch = source_branches[j_source]

    # if source_branch.source

        target_branch = target_branches[i_target]

        # if target_branch.target

            # branch center separation distance
            Δx, Δy, Δz = target_branch.target_center - source_branch.source_center
            separation_distance_squared = Δx*Δx + Δy*Δy + Δz*Δz

            # get r_min, r_max, ρ_min, ρ_max here, based on error method
            # r_min, r_max, ρ_min, ρ_max = get_r_ρ(target_branch, source_branch, separation_distance_squared)

            #    # Barba's multipole acceptance test- slightly different than mine
            #    summed_radii = target_branch.target_radius + source_branch.source_radius
            #    mac = separation_distance * multipole_threshold >= summed_radii # meet M2L criteria

            # decide whether or not to accept the multipole expansion
            summed_radii = source_branch.source_radius + target_branch.target_radius
            # summed_radii = sqrt(3) * mean(source_branch.source_box) + sqrt(3) * mean(target_branch.target_box)

            if separation_distance_squared * multipole_threshold * multipole_threshold > summed_radii * summed_radii
            #if ρ_max <= multipole_threshold * r_min && r_max <= multipole_threshold * ρ_min # exploring a new criterion
                if ff
                    push!(m2l_list, SVector{2}(i_target, j_source))
                end

                return nothing
            end

            fraction = 0.0
            for i_sys in 1:length(source_branch.bodies_index)
                fraction += length(source_branch.bodies_index[i_sys]) / (source_leaf_size[i_sys] * source_leaf_size[i_sys])
            end
            n_targets = 0
            for i_sys in eachindex(target_branch.bodies_index)
                n_targets += length(target_branch.bodies_index[i_sys])
            end
            fraction *= n_targets

            if fraction < 4.0 || source_branch.n_branches == target_branch.n_branches == 0
            # if source_branch.n_branches == target_branch.n_branches == 0 # both leaves
            # elseif source_branch.n_branches == target_branch.n_branches == 0 # both leaves
                nf && (i_target!=j_source || si) && push!(direct_list, SVector{2}(i_target, j_source))

            elseif source_branch.n_branches == 0 || (target_branch.target_radius >= source_branch.source_radius && target_branch.n_branches != 0) # source is a leaf OR target is not a leaf and is bigger or the same size

                for i_child in target_branch.branch_index
                    build_interaction_lists!(m2l_list, direct_list, i_child, j_source, target_branches, source_branches, source_leaf_size, multipole_threshold, farfield, nearfield, self_induced)
                end

            else # source is not a leaf AND target is a leaf or is smaller

                for j_child in source_branch.branch_index
                    build_interaction_lists!(m2l_list, direct_list, i_target, j_child, target_branches, source_branches, source_leaf_size, multipole_threshold, farfield, nearfield, self_induced)
                end

            end
        # end
    # end
end

@inline preallocate_bodies_index(T::Type{<:Branch{<:Any,NT}}, n) where NT = Tuple(Vector{UnitRange{Int64}}(undef, n) for _ in 1:NT)
# @inline preallocate_bodies_index(T::Type{<:SingleBranch}, n) = Vector{UnitRange{Int64}}(undef, n)

function sort_by_target(direct_list, target_branches::Vector{<:Branch})
    # count cardinality of each target leaf in direct_list
    target_counter = zeros(Int32, 2, length(target_branches))
    for (i,j) in direct_list
        target_counter[1,i] += 1
    end

    # cumsum cardinality to obtain an index map
    target_counter[2,1] = Int32(1)
    for i in 2:size(target_counter,2)
        target_counter[2,i] = target_counter[2,i-1] + target_counter[1,i-1]
    end

    # preallocate sorted direct_list
    sorted_direct_list = similar(direct_list)

    # sort direct_list by source
    for ij in direct_list
        # get source branch index
        i = ij[1]

        # get and update target destination index for this branch
        i_dest = target_counter[2,i]
        target_counter[2,i] += Int32(1)

        # place target-source pair in the sorted list
        sorted_direct_list[i_dest] = ij
    end

    return sorted_direct_list
end

function sort_by_source(direct_list, source_branches::Vector{<:Branch})
    # count cardinality of each source leaf in direct_list
    source_counter = zeros(Int32, 2, length(source_branches))
    for (i,j) in direct_list
        source_counter[1,j] += Int32(1)
    end

    # cumsum cardinality to obtain an index map
    source_counter[2,1] = Int32(1)
    for i in 2:size(source_counter,2)
        source_counter[2,i] = source_counter[2,i-1] + source_counter[1,i-1]
    end

    # preallocate sorted direct_list
    sorted_direct_list = similar(direct_list)

    # sort direct_list by source
    for ij in direct_list
        # get source branch index
        j = ij[2]

        # get and update target destination index for this branch
        i_dest = source_counter[2,j]
        source_counter[2,j] += Int32(1)

        # place target-source pair in the sorted list
        sorted_direct_list[i_dest] = ij
    end

    return sorted_direct_list
end

@inline function update_direct_bodies!(direct_bodies::Vector{<:UnitRange}, leaf_index, bodies_index::UnitRange)
    direct_bodies[leaf_index] = bodies_index
end

@inline function update_direct_bodies!(direct_bodies_list, leaf_index, bodies_indices::AbstractVector{<:UnitRange})
    for (direct_bodies, bodies_index) in zip(direct_bodies_list, bodies_indices)
        update_direct_bodies!(direct_bodies, leaf_index, bodies_index)
    end
end

function InteractionList(direct_list, target_systems, target_tree::Tree, source_systems, source_tree::Tree{TF}, derivatives_switches) where TF
    # unpack tree
    leaf_index = source_tree.leaf_index

    # preallocate containers
    influence_matrices = Vector{Matrix{TF}}(undef, length(leaf_index))

    # determine strength dimensions
    strength_dims = get_strength_dims(source_systems)

    # add influence matrices
    for (i_matrix,i_source_branch) in enumerate(leaf_index)
        add_influence_matrix!(influence_matrices, i_matrix, target_systems, target_tree.branches, source_systems, source_tree.branches, i_source_branch, strength_dims, direct_list, derivatives_switches)
    end

    # create largest needed storage strength and influence vectors
    n_cols_max = 0
    n_rows_max = 0
    for influence_matrix in influence_matrices
        n_rows, n_cols = size(influence_matrix)
        n_cols_max = max(n_cols_max, n_cols)
        n_rows_max = max(n_rows_max, n_rows)
    end
    strengths = zeros(TF,n_cols_max)
    influence = zeros(TF,n_rows_max)

    return InteractionList{TF}(influence_matrices, strengths, influence, direct_list)
end

