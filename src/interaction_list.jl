#--- create interaction lists ---#

function build_interaction_lists(target_branches, source_branches, source_leaf_index, multipole_threshold, farfield, nearfield, self_induced, error_method, expansion_order)

    # prepare containers
    m2l_list = Vector{Tuple{Int32,Int32,Int64}}(undef,0)
    direct_list = Vector{SVector{2,Int32}}(undef,0)

    # populate lists
    build_interaction_lists!(m2l_list, direct_list, Int32(1), Int32(1), target_branches, source_branches, source_leaf_index, multipole_threshold, Val(farfield), Val(nearfield), Val(self_induced), error_method, expansion_order)

    return m2l_list, direct_list
end

function get_P(r_min, r_max, ρ_min, ρ_max, expansion_order::Int)
    return expansion_order
end

function get_P(r_min, r_max, ρ_min, ρ_max, ::Dynamic{PMAX,RTOL}) where {PMAX,RTOL}
    return get_P(r_max / ρ_min, ρ_max / r_min, PMAX, RTOL)
end

"""
    get_P(r_max_over_ρ_min, ρ_max_over_r_min, Pmax, ε_rel)

Returns the smallest expansion order not greater than `Pmax` and satisfying the specified relative error tolerance.

# Inputs

* `r_max_over_ρ_min::Float64`: ratio `r_max / ρ_min`
* `ρ_max_over_r_min::Float64`: ratio `ρ_max / r_min`
* `Pmax::Int64`: maximum allowable expansion order
* `ε_rel::Float64`: relative error tolerance

# where

* `r_max::Float64`: distance from the local expansion to the farthest target
* `r_min::Flaot64`: distance from the multipole expansion to the closest target
* `ρ_max::Float64`: distance from the multipole expansion to the farthest source
* `ρ_min::Float64`: distance from the local expansion to the closest source

# Ouputs

* `P::Int`: the smallest expansion order to satisfy the error tolerance

"""
function get_P(r_max_over_ρ_min, ρ_max_over_r_min, Pmax, ε_rel)
    t1 = ρ_max_over_r_min
    t2 = r_max_over_ρ_min
    for P in 0:Pmax-1
        t1 + t2 < ε_rel && (return P)
        t1 *= ρ_max_over_r_min
        t2 *= r_max_over_ρ_min
    end

    @warn "error tolerance not met with dynamic expansion order; try increasing `Pmax`"

    return Pmax
end

function minimum_distance(dx, rx_minus, rx_plus)
    right = dx + rx_plus
    left = dx + rx_minus
    same_sign = right * left >= 0
    return same_sign * min(right, left)
end

function minimum_distance(center1, center2, box2::SVector{6,<:Any})
    x1, y1, z1 = center1
    x2, y2, z2 = center2
    bx_min, bx_max, by_min, by_max, bz_min, bz_max = box2
    Δx = minimum_distance(x2-x1, bx_min, bx_max)
    Δy = minimum_distance(y2-y1, by_min, by_max)
    Δz = minimum_distance(z2-z1, bz_min, bz_max)
    return sqrt(Δx*Δx + Δy*Δy + Δz*Δz)
end

function minimum_distance(center1, center2, box2::SVector{3,<:Any})
    x1, y1, z1 = center1
    x2, y2, z2 = center2
    bx, by, bz = box2
    Δx = minimum_distance(x2-x1, -bx, bx)
    Δy = minimum_distance(y2-y1, -by, by)
    Δz = minimum_distance(z2-z1, -bz, bz)
    return sqrt(Δx*Δx + Δy*Δy + Δz*Δz)
end

@inline function get_r_ρ(target_branch, source_branch, separation_distance, ::UnequalBoxes)
    r_max = target_branch.target_radius
    r_min = minimum_distance(source_branch.center, target_branch.center, target_branch.target_box)
    ρ_max = source_branch.source_radius
    ρ_min = minimum_distance(target_branch.center, source_branch.center, source_branch.source_box)

    return r_min, r_max, ρ_min, ρ_max
end

@inline function get_r_ρ(target_branch, source_branch, separation_distance, ::UnequalSpheres)
    r_max = target_branch.target_radius
    r_min = separation_distance - r_max
    ρ_max = source_branch.source_radius
    ρ_min = separation_distance - ρ_max

    return r_min, r_max, ρ_min, ρ_max
end

@inline function get_r_ρ(target_branch, source_branch, separation_distance, ::EqualSpheres)
    # see Section 3.5.2, Pringle, 1994
    # but normalized by a different mean potential
    # to remove the singularity
    r_max = source_branch.source_radius # should I use source or target radius?
    r_min = separation_distance - r_max
    ρ_max = r_max
    ρ_min = r_min

    return r_min, r_max, ρ_min, ρ_max
end

function build_interaction_lists!(m2l_list, direct_list, i_target, j_source, target_branches, source_branches, source_leaf_index, multipole_threshold, farfield::Val{ff}, nearfield::Val{nf}, self_induced::Val{si}, error_method::ErrorMethod, expansion_order) where {ff,nf,si}
    # unpack
    source_branch = source_branches[j_source]
    target_branch = target_branches[i_target]

    # branch center separation distance
    Δx, Δy, Δz = source_branch.center - target_branch.center
    separation_distance = sqrt(Δx*Δx + Δy*Δy + Δz*Δz)

    # get r_min, r_max, ρ_min, ρ_max here, based on error method
    r_min, r_max, ρ_min, ρ_max = get_r_ρ(target_branch, source_branch, separation_distance, error_method)

    #    # Barba's multipole acceptance test- slightly different than mine
    #    summed_radii = target_branch.target_radius + source_branch.source_radius
    #    mac = separation_distance * multipole_threshold >= summed_radii # meet M2L criteria

    # decide whether or not to accept the multipole expansion
    if ρ_max <= multipole_threshold * r_min && r_max <= multipole_threshold * ρ_min
        if ff
            P = get_P(r_min, r_max, ρ_min, ρ_max, expansion_order)
            push!(m2l_list, (i_target, j_source, P))
        end
    elseif source_branch.n_branches == target_branch.n_branches == 0 # both leaves
        nf && (i_target!=j_source || si) && push!(direct_list, SVector{2}(i_target, j_source))
    elseif source_branch.n_branches == 0 || (target_branch.target_radius >= source_branch.source_radius && target_branch.n_branches != 0) # source is a leaf OR target is not a leaf and is bigger or the same size
        for i_child in target_branch.branch_index
            build_interaction_lists!(m2l_list, direct_list, i_child, j_source, target_branches, source_branches, source_leaf_index, multipole_threshold, farfield, nearfield, self_induced, error_method, expansion_order)
        end
    else # source is not a leaf AND target is a leaf or is smaller
        for j_child in source_branch.branch_index
            build_interaction_lists!(m2l_list, direct_list, i_target, j_child, target_branches, source_branches, source_leaf_index, multipole_threshold, farfield, nearfield, self_induced, error_method, expansion_order)
        end
    end
end

@inline preallocate_bodies_index(T::Type{<:MultiBranch{<:Any,NT}}, n) where NT = Tuple(Vector{UnitRange{Int64}}(undef, n) for _ in 1:NT)
@inline preallocate_bodies_index(T::Type{<:SingleBranch}, n) = Vector{UnitRange{Int64}}(undef, n)

function sort_list_by_target(direct_list, target_branches::Vector{TT}, source_branches::Vector{TS}, n_leaves) where {TT,TS}
    target_counter = zeros(Int32, n_leaves)
    place_counter = zeros(Int32, n_leaves)
    direct_target_bodies = preallocate_bodies_index(TT, length(direct_list))
    direct_source_bodies = preallocate_bodies_index(TS, length(direct_list))

    # tally the contributions of each source
    for (i_target, j_source) in direct_list
        i_leaf = source_branches[i_target].i_leaf
        target_counter[i_leaf] += 1
    end

    # prepare place counter
    i_cum = 1
    for (i,n) in enumerate(target_counter)
        place_counter[i] = i_cum
        i_cum += n
    end

    # place interactions
    for (i,(i_target, j_source)) in enumerate(direct_list)
        i_leaf = target_branches[i_target].i_leaf
        update_direct_bodies!(direct_target_bodies, place_counter[i_leaf], target_branches[i_target].bodies_index)
        update_direct_bodies!(direct_source_bodies, place_counter[i_leaf], source_branches[j_source].bodies_index)
        place_counter[i_leaf] += 1
    end

    return direct_target_bodies, direct_source_bodies
end

function sort_list_by_source(direct_list, target_branches::Vector{TT}, source_branches::Vector{TS}, n_leaves) where {TT,TS}
    source_counter = zeros(Int32, n_leaves)
    place_counter = zeros(Int32, n_leaves)
    direct_target_bodies = preallocate_bodies_index(TT, length(direct_list))
    direct_source_bodies = preallocate_bodies_index(TS, length(direct_list))

    # tally the contributions of each source
    for (i_target, j_source) in direct_list
        j_leaf = source_branches[j_source].i_leaf
        source_counter[j_leaf] += 1
    end

    # prepare place counter
    i_cum = 1
    for (i,n) in enumerate(source_counter)
        place_counter[i] = i_cum
        i_cum += n
    end

    # place interactions
    for (i, (i_target, j_source)) in enumerate(direct_list)
        i_leaf = target_branches[i_target].i_leaf
        j_leaf = source_branches[j_source].i_leaf
        update_direct_bodies!(direct_target_bodies, place_counter[j_leaf], target_branches[i_target].bodies_index)
        update_direct_bodies!(direct_source_bodies, place_counter[j_leaf], source_branches[j_source].bodies_index)
        place_counter[j_leaf] += 1
    end

    return direct_target_bodies, direct_source_bodies
end


@inline function update_direct_bodies!(direct_bodies::Vector{<:UnitRange}, leaf_index, bodies_index::UnitRange)
    direct_bodies[leaf_index] = bodies_index
end

@inline function update_direct_bodies!(direct_bodies_list, leaf_index, bodies_indices::AbstractVector{<:UnitRange})
    for (direct_bodies, bodies_index) in zip(direct_bodies_list, bodies_indices)
        update_direct_bodies!(direct_bodies, leaf_index, bodies_index)
    end
end

function InteractionList(direct_list, target_systems, target_tree::Tree, source_systems, source_tree::Tree{TF,<:Any}, derivatives_switches) where TF
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

#--- dispatching dynamic expansion order ---#

@inline function get_Pmax(expansion_order::Int)
    return expansion_order
end

@inline function get_Pmax(expansion_order::Dynamic{PMAX,<:Any}) where PMAX
    return PMAX
end

