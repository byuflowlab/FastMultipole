#--- determine distances for error formulae ---#

@inline minabs(x, y) = x * (abs(x) <= abs(y)) + y * (abs(y) < abs(x))

function minimum_distance(dx, rx_minus, rx_plus)
    right = dx + rx_plus
    left = dx + rx_minus
    same_sign = right * left >= 0
    return same_sign * minabs(right, left)
end

function minimum_distance(center1, center2, box2::SVector{3,<:Any})
    x1, y1, z1 = center1
    x2, y2, z2 = center2
    bx, by, bz = box2
    Δx = minimum_distance(x2-x1, -bx, bx)
    Δy = minimum_distance(y2-y1, -by, by)
    Δz = minimum_distance(z2-z1, -bz, bz)
    return SVector{3}(Δx, Δy, Δz)
end

function closest_corner(center1, center2, box2::SVector{3,<:Any})
    Δxc, Δyc, Δzc = center1 - center2
    bx, by, bz = box2
    Δx = minabs(Δxc-bx,Δxc+bx)
    Δy = minabs(Δyc-by,Δyc+by)
    Δz = minabs(Δzc-bz,Δzc+bz)
    return SVector{3}(Δxc-Δx, Δyc-Δy, Δzc-Δz)
end

@inline function get_r_ρ(local_branch, multipole_branch, separation_distance_squared, ::UnequalBoxes)
    r_max = local_branch.target_radius
    r_min = norm(minimum_distance(multipole_branch.source_center, local_branch.target_center, local_branch.target_box))
    ρ_max = multipole_branch.source_radius
    ρ_min = norm(minimum_distance(local_branch.target_center, multipole_branch.source_center, multipole_branch.source_box))

    return r_min, r_max, ρ_min, ρ_max
end

@inline function get_r_ρ(local_branch, multipole_branch, separation_distance_squared, ::UniformUnequalBoxes)
    r_max = local_branch.target_radius
    r_min = norm(minimum_distance(multipole_branch.source_center, local_branch.target_center, local_branch.target_box))
    ρ_max = multipole_branch.source_radius
    ρ_min = sqrt(separation_distance_squared) - ρ_max

    return r_min, r_max, ρ_min, ρ_max
end

@inline function get_r_ρ(local_branch, multipole_branch, separation_distance_squared, ::Union{UnequalSpheres, UniformUnequalSpheres})
    separation_distance = sqrt(separation_distance_squared)
    r_max = local_branch.target_radius
    r_min = separation_distance - r_max
    ρ_max = multipole_branch.source_radius
    ρ_min = separation_distance - ρ_max

    return r_min, r_max, ρ_min, ρ_max
end

@inline function get_r_ρ(local_branch, multipole_branch, separation_distance_squared, ::EqualSpheres)
    # see Section 3.5.2, Pringle, 1994
    # but normalized by a different mean potential
    # to remove the singularity
    separation_distance = sqrt(separation_distance_squared)
    r_max = multipole_branch.source_radius # should I use source or target radius?
    r_min = separation_distance - r_max
    ρ_max = r_max
    ρ_min = r_min

    return r_min, r_max, ρ_min, ρ_max
end

function dipole_from_multipole(expansion)
    # dipole term
    qx = -2 * expansion[2,1,3]
    qy = -2 * expansion[1,1,3]
    qz = expansion[1,1,2]
    return qx, qy, qz
end

function get_A_local(multipole_branch, r_vec, expansion_order)
    # extract expansion
    expansion = multipole_branch.multipole_expansion

    # monopole term
    q_monopole = expansion[1,1,1]

    # dipole term
    q_dipole_x, q_dipole_y, q_dipole_z = dipole_from_multipole(expansion)

    # scaling factor
    rx, ry, rz = r_vec
    r2 = rx*rx + ry*ry + rz*rz
    Q = max(abs(q_monopole), (expansion_order+1) * abs(q_dipole_x*r_vec[1] + q_dipole_y*r_vec[2] + q_dipole_z*r_vec[3]) / r2)

    Q = q_monopole
    return Q
end

#--- actual error prediction functions ---#

function multipole_error(local_branch, multipole_branch, P, error_method::Union{EqualSpheres, UnequalSpheres, UnequalBoxes})
    # get distances
    Δx, Δy, Δz = local_branch.target_center - multipole_branch.source_center
    separation_distance_squared = Δx*Δx + Δy*Δy + Δz*Δz
    r_min, _, _, ρ_max = get_r_ρ(local_branch, multipole_branch, separation_distance_squared, error_method)

    # calculate error upper bound
    A = multipole_branch.multipole_expansion[1,1,1]
    ρ_max_over_r_min = ρ_max / r_min
    ε = ρ_max_over_r_min * A / (r_min - ρ_max)
    for p in 1:P
        ε *= ρ_max_over_r_min
    end

    return abs(ε) * ONE_OVER_4π
end

function multipole_error(local_branch, multipole_branch, P, error_method::Union{UniformUnequalSpheres, UniformUnequalBoxes})
    return 3 / (2*P+8) * multipole_error(local_branch, multipole_branch, P, UnequalSpheres())
    #return 3 / (2*P+8) * sqrt(2/(2*P+3)) * multipole_error(local_branch, multipole_branch, P, UnequalSpheres())
end

function multipole_error(local_branch, multipole_branch, P, error_method::UniformCubesVelocity)
    # extract fields
    source_box = multipole_branch.source_box
    target_box = local_branch.target_box
    source_center = multipole_branch.source_center
    target_center = local_branch.target_center

    # choose the closest point in the local branch
    r_vec = minimum_distance(source_center, target_center, target_box)

    # length scale
    iszero(source_box) && (return zero(r_vec[1]))
    s = (source_box[1] + source_box[2] + source_box[3]) * 0.6666666666666666

    # monopole term
    q_monopole = multipole_branch.multipole_expansion[1,1,1]

    # dipole term
    q_dipole_x, q_dipole_y, q_dipole_z = dipole_from_multipole(multipole_branch.multipole_expansion)
    q_dipole_2 = q_dipole_x * q_dipole_x + q_dipole_y * q_dipole_y + q_dipole_z * q_dipole_z

    # determine whether monopole/dipole dominates
    rx, ry, rz = r_vec
    r2 = rx*rx + ry*ry + rz*rz
    r_inv_2 = 1/r2
    r_inv = sqrt(r_inv_2)
    q_dot_r̂ = (q_dipole_x*r_vec[1] + q_dipole_y*r_vec[2] + q_dipole_z*r_vec[3]) * r_inv
    Q = max(abs(q_monopole), r_inv * sqrt(q_dot_r̂ * q_dot_r̂ * (P+3) * (P+3) + r_inv_2 * (q_dipole_2 - q_dot_r̂ * q_dot_r̂)))

    # get polar/azimuth angles for accessing error database
    _, θr, ϕr = cartesian_to_spherical(r_vec)
    iθ = get_iθ(θr)
    iϕ = get_iϕ(ϕr)

    # get s^(P+4)/r^(P+1)
    s_rinv = s * r_inv
    s2 = s*s
    scalar = s_rinv / s
    for n in 1:P
        scalar *= s_rinv
    end

    # n_max = P+1

    # compute error
    # ε = 0.0
    # for n in P+1:n_max
        # ε += FastMultipole.MULTIPOLE_INTEGRALS[n,iθ,iϕ] * scalar
        # scalar *= s_rinv
    # end
    ε = FastMultipole.MULTIPOLE_INTEGRALS[P+1,iθ,iϕ] * scalar * Q * (P+2) * r_inv

    return abs(ε) * ONE_OVER_4π
end

function multipole_error(local_branch, multipole_branch, P, error_method::UniformCubes)
    # extract fields
    source_box = multipole_branch.source_box
    target_box = local_branch.target_box
    source_center = multipole_branch.source_center
    target_center = local_branch.target_center

    # choose the closest point in the local branch
    r_vec = minimum_distance(source_center, target_center, target_box)

    # length scale
    iszero(source_box) && (return zero(r_vec[1]))
    s = (source_box[1] + source_box[2] + source_box[3]) * 0.6666666666666666

    # monopole term
    q_monopole = multipole_branch.multipole_expansion[1,1,1]

    # dipole term
    q_dipole_x, q_dipole_y, q_dipole_z = dipole_from_multipole(multipole_branch.multipole_expansion)

    # determine whether monopole/dipole dominates
    rx, ry, rz = r_vec
    r2 = rx*rx + ry*ry + rz*rz
    r_inv = 1/sqrt(r2)
    Q = max(abs(q_monopole), (P+2) * abs(q_dipole_x*r_vec[1] + q_dipole_y*r_vec[2] + q_dipole_z*r_vec[3]) * r_inv * r_inv)
    # if Q == abs(q_monopole)
    #     println("monopole")
    # end

    # get polar/azimuth angles for accessing error database
    _, θr, ϕr = cartesian_to_spherical(r_vec)
    iθ = get_iθ(θr)
    iϕ = get_iϕ(ϕr)

    # get s^(P+4)/r^(P+1)
    s_rinv = s * r_inv
    s2 = s*s
    scalar = s_rinv / s
    for n in 1:P
        scalar *= s_rinv
    end

    # n_max = P+1

    # compute error
    # ε = 0.0
    # for n in P+1:n_max
        # ε += FastMultipole.MULTIPOLE_INTEGRALS[n,iθ,iϕ] * scalar
        # scalar *= s_rinv
    # end
    ε = FastMultipole.MULTIPOLE_INTEGRALS[P+1,iθ,iϕ] * scalar * Q

    return abs(ε) * ONE_OVER_4π
end

function local_error(local_branch, multipole_branch, P, error_method::Union{EqualSpheres, UnequalSpheres, UnequalBoxes})
    # get distances
    Δx, Δy, Δz = local_branch.target_center - multipole_branch.source_center
    separation_distance_squared = Δx*Δx + Δy*Δy + Δz*Δz
    r_min, r_max, ρ_min, ρ_max = get_r_ρ(local_branch, multipole_branch, separation_distance_squared, error_method)

    # calculate error upper bound
    r_max_over_ρ_min = r_max / ρ_min
    ε = r_max_over_ρ_min

    A = multipole_branch.multipole_expansion[1,1,1]
    ε *= A / (ρ_min - r_max)

    for p in 1:P
        ε *= r_max_over_ρ_min
    end

    return abs(ε) * ONE_OVER_4π
end

function local_error(local_branch, multipole_branch, P, error_method::Union{UniformUnequalSpheres, UniformUnequalBoxes})
    # get distances
    Δx, Δy, Δz = local_branch.target_center - multipole_branch.source_center
    separation_distance_squared = Δx*Δx + Δy*Δy + Δz*Δz
    r_min, r_max, ρ_min, ρ_max = get_r_ρ(local_branch, multipole_branch, separation_distance_squared, error_method)

    # local error
    ρ_max2 = ρ_max * ρ_max
    r_bar = r_min - ρ_max
    A = multipole_branch.multipole_expansion[1,1,1]
    Γ = 3 * r_max * A / (2 * ρ_max2 * ρ_max)

    # distance from multipole center to point closest to the local expansion
    ΔC = sqrt(separation_distance_squared)
    ρ_max_line = ΔC - ρ_min

    # distance from local center to farthest multipole point
    ΔC_plus = ρ_min + ρ_max_line + ρ_max_line

    # distance from local center to nearest multipole point
    ΔC_minus = ρ_min

    t_plus = r_max / ΔC_plus
    one_over_ΔC_minus = 1/ΔC_minus
    t_minus = r_max * one_over_ΔC_minus
    L = log((1-t_plus)/(1-t_minus))
    ΔC2_inv = 1/(2*ΔC)
    ρ2_ΔC2_2ΔC = (ρ_max2 - separation_distance_squared) * ΔC2_inv
    r2_ΔC2 = r_max * r_max * ΔC2_inv

    # recursively tracked quantities
    Lζ_sum_pm1 = L
    Lζ_sum_pm2 = L
    Lζ_sum_pm3 = L
    t_plus_n = t_plus
    t_minus_n = t_minus
    η_p = log(ΔC_plus * one_over_ΔC_minus)
    r_inv = 1/r_max
    η_pm1 = η_p + (ΔC_plus - ΔC_minus) * r_inv
    η_pm2 = η_pm1 + (ΔC_plus*ΔC_plus - ΔC_minus*ΔC_minus) * r_inv * r_inv * 0.5

    #--- test p=0 ---#

    # test error
    ε_local = Γ * (r_max * (η_pm1 + Lζ_sum_pm2) + ρ2_ΔC2_2ΔC * (η_p + Lζ_sum_pm1) - r2_ΔC2 * (η_pm2 + Lζ_sum_pm3))

    if P==0
        return ε_local
    end

    # recurse local
    η_pm2 = η_pm1
    η_pm1 = η_p
    η_p = zero(r_min)

    #--- test p>0 ---#
    for p in 1:P
        # get local error
        ε_local = Γ * (r_max * (η_pm1 + Lζ_sum_pm2) + ρ2_ΔC2_2ΔC * (η_p + Lζ_sum_pm1) - r2_ΔC2 * (η_pm2 + Lζ_sum_pm3))

        # recurse local
        Lζ_sum_pm3 = Lζ_sum_pm2
        Lζ_sum_pm2 = Lζ_sum_pm1
        Lζ_sum_pm1 += (t_plus_n - t_minus_n) / p
        t_plus_n *= t_plus
        t_minus_n *= t_minus
        η_pm2 = η_pm1
        η_pm1 = η_p
        η_p = zero(r_min)
    end

    ε_local = Γ * (r_max * (η_pm1 + Lζ_sum_pm2) + ρ2_ΔC2_2ΔC * (η_p + Lζ_sum_pm1) - r2_ΔC2 * (η_pm2 + Lζ_sum_pm3))

    return abs(ε_local) * ONE_OVER_4π
end

function local_error(local_branch, multipole_branch, expansion_order, error_method::UniformCubes)
    # extract fields
    source_box = multipole_branch.source_box
    target_box = local_branch.target_box
    source_center = multipole_branch.source_center
    target_center = local_branch.target_center
    P = expansion_order

    # choose the closest point in the local branch
    r_vec = closest_corner(source_center, target_center, target_box)
    rx, ry, rz = r_vec
    r = sqrt(rx*rx + ry*ry + rz*rz)
    r_inv = 1/r

    # monopole term
    q_monopole = multipole_branch.multipole_expansion[1,1,1]

    # dipole term
    q_dipole_x, q_dipole_y, q_dipole_z = dipole_from_multipole(multipole_branch.multipole_expansion)

    # determine whether monopole/dipole dominates
    q_dot_r̂ = (q_dipole_x * rx + q_dipole_y * ry + q_dipole_z * rz) * r_inv
    Q = max(abs(q_monopole), (P+1) * q_dot_r̂ * r_inv)

    # get polar/azimuth angles for accessing error database
    s = sum(target_box) * 0.666666666666666
    Rx, Ry, Rz = source_center - target_center
    R = sqrt(Rx*Rx + Ry*Ry + Rz*Rz)
    Rinv = 1/R
    ω = asin(s * 0.5 * Rinv)
    iω = get_iω(ω)
    cθ = (Rx*rx + Ry*ry + Rz*rz) * Rinv / r # cos(θ)
    γ = acos(cθ)
    iγ = get_iγ(γ)

    # prepare recursive factors
    r_over_R = r * Rinv
    scalar = Rinv * r_over_R
    for n in 1:P
        scalar *= r_over_R
    end

    # calculate error
    ε = 0.0
    for n in P+1:P+2
        ε += LOCAL_INTEGRALS[n,iω,iγ] * scalar
        scalar *= r_over_R
    end

    # scale by charge and return
    return abs(ε) * Q * ONE_OVER_4π
end

function local_error(local_branch, multipole_branch, expansion_order, error_method::UniformCubesVelocity)
    # extract fields
    source_box = multipole_branch.source_box
    target_box = local_branch.target_box
    source_center = multipole_branch.source_center
    target_center = local_branch.target_center
    P = expansion_order

    # choose the closest point in the local branch
    r_vec = closest_corner(source_center, target_center, target_box)
    rx, ry, rz = r_vec
    r = sqrt(rx*rx + ry*ry + rz*rz)
    r_inv = 1/r

    # monopole term
    q_monopole = multipole_branch.multipole_expansion[1,1,1]

    # dipole term
    q_dipole_x, q_dipole_y, q_dipole_z = dipole_from_multipole(multipole_branch.multipole_expansion)
    q_dipole_2 = q_dipole_x * q_dipole_x + q_dipole_y * q_dipole_y + q_dipole_z * q_dipole_z

    # determine whether monopole/dipole dominates
    q_dot_r̂ = (q_dipole_x * rx + q_dipole_y * ry + q_dipole_z * rz) * r_inv
    Q = max(abs(q_monopole), r_inv * sqrt((P*P-1) * q_dot_r̂ * q_dot_r̂ + q_dipole_2))

    # get polar/azimuth angles for accessing error database
    s = sum(target_box) * 0.666666666666666
    Rx, Ry, Rz = source_center - target_center
    R = sqrt(Rx*Rx + Ry*Ry + Rz*Rz)
    Rinv = 1/R
    ω = asin(s * 0.5 * Rinv)
    iω = get_iω(ω)
    cθ = (Rx*rx + Ry*ry + Rz*rz) * Rinv / r # cos(θ)
    γ = acos(cθ)
    iγ = get_iγ(γ)

    # prepare recursive factors
    r_over_R = r * Rinv
    scalar = Rinv * r_over_R
    for n in 1:P
        scalar *= r_over_R
    end

    # calculate error
    ε = 0.0
    for n in P+1:P+3
        ε += LOCAL_INTEGRALS[n,iω,iγ] * scalar * n
        scalar *= r_over_R
    end

    # scale by charge and return
    return abs(ε) * Q * ONE_OVER_4π * r_inv # * r # * r to account for lamb-helmholtz decomp
end

function error(local_branch, multipole_branch, expansion_order, error_method)
    ε_multipole = multipole_error(local_branch, multipole_branch, expansion_order, error_method)
    ε_local = local_error(local_branch, multipole_branch, expansion_order, error_method)
    return max(ε_multipole, ε_local)
end

