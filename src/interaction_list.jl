#--- create interaction lists ---#

function build_interaction_lists(target_branches, source_branches, source_leaf_index, multipole_threshold, farfield, nearfield, self_induced, error_method, expansion_order)

    # prepare containers
    m2l_list = Vector{Tuple{Int32,Int32,Int64}}(undef,0)
    direct_list = Vector{SVector{2,Int32}}(undef,0)

    # populate lists
    build_interaction_lists!(m2l_list, direct_list, Int32(1), Int32(1), target_branches, source_branches, source_leaf_index, multipole_threshold, Val(farfield), Val(nearfield), Val(self_induced), error_method, expansion_order)

    return m2l_list, direct_list
end

function get_P(r_min, r_max, ρ_min, ρ_max, ΔC2, expansion_order::Int, error_method::ErrorMethod)
    return expansion_order
end

function get_P(r_min, r_max, ρ_min, ρ_max, ΔC2, ::Dynamic{PMAX,RTOL}, error_method::ErrorMethod) where {PMAX,RTOL}
    return get_P(r_min, r_max, ρ_min, ρ_max, ΔC2, PMAX, RTOL, error_method)
end

@inline function warn_Pmax()
    if WARNING_FLAG_PMAX[]
        @warn "error tolerance not met with dynamic expansion order; try increasing `Pmax`"
        WARNING_FLAG_PMAX[] = false
    end
end


"""
    get_P(r_min, r_max, ρ_min, ρ_max, ΔC2, Pmax, ε_rel, error_method)

Returns the smallest expansion order not greater than `Pmax` and satisfying the specified relative error tolerance.

# Inputs

* `r_min::Flaot64`: distance from the multipole expansion to the closest target
* `r_max::Float64`: distance from the local expansion to the farthest target
* `ρ_min::Float64`: distance from the local expansion to the closest source
* `ρ_max::Float64`: distance from the multipole expansion to the farthest source
* `ΔC2::Float64`: distance squared between multipole and local centers
* `Pmax::Int64`: maximum allowable expansion order
* `ε_rel::Float64`: relative error tolerance
* `error_method::ErrorMethod`: type used to dispatch on the desired error method

# Ouputs

* `P::Int`: the smallest expansion order to satisfy the error tolerance

"""
function get_P(r_min, r_max, ρ_min, ρ_max, ΔC2, Pmax, ε_rel, ::Union{EqualSpheres, UnequalSpheres, UnequalBoxes})
    ρ_max_over_r_min = ρ_max / r_min
    r_max_over_ρ_min = r_max / ρ_min
    t1 = ρ_max_over_r_min
    t2 = r_max_over_ρ_min
    for P in 0:Pmax-1
        t1 + t2 < ε_rel && (return P)
        t1 *= ρ_max_over_r_min
        t2 *= r_max_over_ρ_min
    end

    warn_Pmax()

    return Pmax
end

function get_P(r_min, r_max, ρ_min, ρ_max, ΔC2, Pmax, ε_rel, ::Union{UniformUnequalSpheres, UniformUnequalBoxes})

    ρ_max_over_r_min = ρ_max / r_min
    r_max_over_ρ_min = r_max / ρ_min

    # multipole error
    ε_multipole = 1.5 * ρ_max / (r_min - ρ_max)

    # local error
    ρ_max2 = ρ_max * ρ_max
    Γ = 3 * r_max * r_min / (2 * ρ_max2 * ρ_max)

    # distance from multipole center to point closest to the local expansion
    ΔC = sqrt(ΔC2)
    ρ_max_line = ΔC - ρ_min

    # distance from local center to farthest multipole point
    ΔC_plus = ρ_min + ρ_max_line + ρ_max_line

    # distance from local center to nearest multipole point
    ΔC_minus = ρ_min

    t_plus = r_max / ΔC_plus
    t_minus = r_max / ΔC_minus
    L = log((1-t_plus)/(1-t_minus))
    ΔC2_inv = 1/(2*ΔC)
    ρ2_ΔC2_2ΔC = (ρ_max2 - ΔC2) * ΔC2_inv
    r2_ΔC2 = r_max * r_max * ΔC2_inv

    # recursively tracked quantities
    Lζ_sum_pm1 = L
    Lζ_sum_pm2 = L
    Lζ_sum_pm3 = L
    t_plus_n = t_plus
    t_minus_n = t_minus
    η_p = log(ΔC_plus / ΔC_minus)
    r_inv = 1/r_max
    η_pm1 = η_p + (ΔC_plus - ΔC_minus) * r_inv
    η_pm2 = η_pm1 + (ΔC_plus*ΔC_plus - ΔC_minus*ΔC_minus) * r_inv * r_inv * 0.5

    #--- test p=0 ---#

    # test error
    ε_local = Γ * (r_max * (η_pm1 + Lζ_sum_pm2) + ρ2_ΔC2_2ΔC * (η_p + Lζ_sum_pm1) - r2_ΔC2 * (η_pm2 + Lζ_sum_pm3))
    ε_multipole * 0.25 + ε_local < ε_rel && (return 0)

    # recurse multipole
    ε_multipole *= ρ_max_over_r_min

    # recurse local
    η_pm2 = η_pm1
    η_pm1 = η_p
    η_p = zero(r_min)

    #--- test p>0 ---#

    for P in 1:Pmax-1
        # get local error
        ε_local = Γ * (r_max * (η_pm1 + Lζ_sum_pm2) + ρ2_ΔC2_2ΔC * (η_p + Lζ_sum_pm1) - r2_ΔC2 * (η_pm2 + Lζ_sum_pm3))

        # check error
        ε_multipole / (P+4) + ε_local < ε_rel && (return P)

        # recurse multipole
        ε_multipole *= ρ_max_over_r_min

        # recurse local
        Lζ_sum_pm3 = Lζ_sum_pm2
        Lζ_sum_pm2 = Lζ_sum_pm1
        Lζ_sum_pm1 += (t_plus_n - t_minus_n) / P
        t_plus_n *= t_plus
        t_minus_n *= t_minus
        η_pm2 = η_pm1
        η_pm1 = η_p
        η_p = zero(r_min)
    end

    warn_Pmax()

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

@inline function get_r_ρ(target_branch, source_branch, separation_distance_squared, ::Union{UnequalBoxes, UniformUnequalBoxes})
    r_max = target_branch.target_radius
    r_min = minimum_distance(source_branch.center, target_branch.center, target_branch.target_box)
    ρ_max = source_branch.source_radius
    ρ_min = minimum_distance(target_branch.center, source_branch.center, source_branch.source_box)

    return r_min, r_max, ρ_min, ρ_max
end

@inline function get_r_ρ(target_branch, source_branch, separation_distance_squared, ::Union{UnequalSpheres, UniformUnequalSpheres})
    separation_distance = sqrt(separation_distance_squared)
    r_max = target_branch.target_radius
    r_min = separation_distance - r_max
    ρ_max = source_branch.source_radius
    ρ_min = separation_distance - ρ_max

    return r_min, r_max, ρ_min, ρ_max
end

@inline function get_r_ρ(target_branch, source_branch, separation_distance_squared, ::EqualSpheres)
    # see Section 3.5.2, Pringle, 1994
    # but normalized by a different mean potential
    # to remove the singularity
    separation_distance = sqrt(separation_distance_squared)
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
    separation_distance_squared = Δx*Δx + Δy*Δy + Δz*Δz

    # get r_min, r_max, ρ_min, ρ_max here, based on error method
    r_min, r_max, ρ_min, ρ_max = get_r_ρ(target_branch, source_branch, separation_distance_squared, error_method)

    #    # Barba's multipole acceptance test- slightly different than mine
    #    summed_radii = target_branch.target_radius + source_branch.source_radius
    #    mac = separation_distance * multipole_threshold >= summed_radii # meet M2L criteria

    # decide whether or not to accept the multipole expansion
    summed_radii = r_max + ρ_max
    if separation_distance_squared * multipole_threshold * multipole_threshold > summed_radii * summed_radii
    #if ρ_max <= multipole_threshold * r_min && r_max <= multipole_threshold * ρ_min # exploring a new criterion
        if ff
            P = get_P(r_min, r_max, ρ_min, ρ_max, separation_distance_squared, expansion_order, error_method)
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

function sort_list_by_target(direct_list, target_branches::Vector{<:Branch})
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

function sort_list_by_source(direct_list, source_branches::Vector{<:Branch})
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

