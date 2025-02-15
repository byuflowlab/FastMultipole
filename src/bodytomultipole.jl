#------- multipole generation function -------#

#=

Convenience functions are accessed by overloading as follows:

body_to_multipole!(system::<mysystemtype>, args...) = body_to_multipole!(<geometry>{<kernel>}, system, args...)

For example, if my system type were called VortexSheetSurface and consisted of vortex panels, I would overload

body_to_multipole!(system::VortexSheetSurface, args...) = body_to_multipole!(Panel{Vortex}, system, args...)

=#

# @inline function body_to_multipole!(branch::SingleBranch, system, harmonics, expansion_order)
#     body_to_multipole!(system, branch, branch.bodies_index, harmonics, expansion_order)
# end

function body_to_multipole!(branches::Vector{<:Branch}, leaf_index, system, i_system, expansion_order)
    for i_branch in leaf_index
        branch = branches[i_branch]
        body_to_multipole!(system, branch, branch.bodies_index[i_system], branch.harmonics, expansion_order)
    end
end

function body_to_multipole!(branch::Branch, systems::Tuple, harmonics, expansion_order)
    # iterate over systems
    for (system, bodies_index) in zip(systems, branch.bodies_index)
        body_to_multipole!(system, branch, bodies_index, harmonics, expansion_order)
    end
end

#------- multipole generation convenience functions -------#

#--- helper functions for integrating over simplices as in Gumerov 2023 ---#

function xyz_to_ξηz(x)
    ξ_real = x[1] * 0.5
    ξ_imag = x[2] * 0.5
    η_real = ξ_real
    η_imag = -ξ_imag
    z = x[3]
    return ξ_real, ξ_imag, η_real, η_imag, z
end

function get_n(harmonics, index, n, m, i, _1_m)
    if m > 0
        x_n_mm1_real = harmonics[1,index,i-1]
        x_n_mm1_imag = harmonics[2,index,i-1]
    elseif m == 0
        if n == 0
            x_n_mm1_real = zero(eltype(harmonics))
            x_n_mm1_imag = zero(eltype(harmonics))
        else
            x_n_mm1_real = -_1_m * harmonics[1,index,i+1]
            x_n_mm1_imag = _1_m * harmonics[2,index,i+1]
        end
    end

    x_n_m_real = harmonics[1,index,i]
    x_n_m_imag = harmonics[2,index,i]

    if m < n
        x_n_mp1_real = harmonics[1,index,i+1]
        x_n_mp1_imag = harmonics[2,index,i+1]
    else
        x_n_mp1_real = zero(eltype(harmonics))
        x_n_mp1_imag = zero(eltype(harmonics))
    end

    return x_n_mm1_real, x_n_mm1_imag, x_n_m_real, x_n_m_imag, x_n_mp1_real, x_n_mp1_imag
end


function get_nm1(harmonics, index, n, m, i_nm1_m, _1_m)
    if m > 0
        @assert 0 <= m-1 <= n-1
        x_nm1_mm1_real = harmonics[1,index,i_nm1_m-1]
        x_nm1_mm1_imag = harmonics[2,index,i_nm1_m-1]
    elseif n > 1
        @assert m-1==-1
        x_nm1_mm1_real = -_1_m * harmonics[1,index,i_nm1_m+1]
        x_nm1_mm1_imag = _1_m * harmonics[2,index,i_nm1_m+1]
    else
        x_nm1_mm1_real = zero(eltype(harmonics))
        x_nm1_mm1_imag = zero(eltype(harmonics))
    end

    if m < n
        @assert 0 <= m <= n-1
        x_nm1_m_real = harmonics[1,index,i_nm1_m]
        x_nm1_m_imag = harmonics[2,index,i_nm1_m]
    else
        x_nm1_m_real = zero(eltype(harmonics))
        x_nm1_m_imag = zero(eltype(harmonics))
    end

    if m+1 < n
        @assert 0 <= m+1 <= n-1
        x_nm1_mp1_real = harmonics[1,index,i_nm1_m+1]
        x_nm1_mp1_imag = harmonics[2,index,i_nm1_m+1]
    else
        x_nm1_mp1_real = zero(eltype(harmonics))
        x_nm1_mp1_imag = zero(eltype(harmonics))
    end

    return x_nm1_mm1_real, x_nm1_mm1_imag, x_nm1_m_real, x_nm1_m_imag, x_nm1_mp1_real, x_nm1_mp1_imag
end

function calculate_q!(harmonics, ξ_real, ξ_imag, η_real, η_imag, z, expansion_order)

    #--- n = 0, m = 0 ---#

    i = 1
    harmonics[1,1,i] = one(eltype(harmonics))
    harmonics[2,1,i] = zero(eltype(harmonics))
    i += 1

    #--- n > 0 ---#
    for n in 1:expansion_order
        i_nm1_m = i - n
        _1_m = 1.0
        for m in 0:n
            # get n-1 values
            q_nm1_mm1_real, q_nm1_mm1_imag, q_nm1_m_real, q_nm1_m_imag, q_nm1_mp1_real, q_nm1_mp1_imag = get_nm1(harmonics, 1, n, m, i_nm1_m, _1_m)

            # set coefficient
            harmonics[1,1,i] = (-(ξ_real * q_nm1_mm1_imag + ξ_imag * q_nm1_mm1_real) - (η_real * q_nm1_mp1_imag + η_imag * q_nm1_mp1_real) - z * q_nm1_m_real) / n
            harmonics[2,1,i] = ((ξ_real * q_nm1_mm1_real - ξ_imag * q_nm1_mm1_imag) + (η_real * q_nm1_mp1_real - η_imag * q_nm1_mp1_imag) - z * q_nm1_m_imag) / n

            # recurse
            i += 1
            i_nm1_m += 1
            _1_m = -_1_m
        end
    end
end

"""
assumes q has already been calculated

If X0_real is replaced with X0_real + Xv_real, etc., the result becomes jnm instead of pnm, as used for panels.
"""
function calculate_pj!(harmonics, ξ0_real, ξ0_imag, η0_real, η0_imag, z0, expansion_order)
    #--- n = 0, m = 0 ---#

    i = 1
    harmonics[1,2,1] = one(eltype(harmonics))
    harmonics[2,2,1] = zero(eltype(harmonics))
    i += 1

    #--- n > 0 ---#
    for n in 1:expansion_order
        i_nm1_m = i - n
        _1_m = 1
        for m in 0:n
            # get n-1 values
            p_nm1_mm1_real, p_nm1_mm1_imag, p_nm1_m_real, p_nm1_m_imag, p_nm1_mp1_real, p_nm1_mp1_imag = get_nm1(harmonics, 2, n, m, i_nm1_m, _1_m)

            # get qnm
            qnm_real = harmonics[1,1,i]
            qnm_imag = harmonics[2,1,i]

            # set coefficient
            harmonics[1,2,i] = (-(ξ0_real * p_nm1_mm1_imag + ξ0_imag * p_nm1_mm1_real) - (η0_real * p_nm1_mp1_imag + η0_imag * p_nm1_mp1_real) - z0 * p_nm1_m_real + qnm_real) / (n+1)
            harmonics[2,2,i] = ((ξ0_real * p_nm1_mm1_real - ξ0_imag * p_nm1_mm1_imag) + (η0_real * p_nm1_mp1_real - η0_imag * p_nm1_mp1_imag) - z0 * p_nm1_m_imag + qnm_imag) / (n+1)

            # recurse
            i += 1
            i_nm1_m += 1
            _1_m = -_1_m
        end
    end
end

"""
assumes j has already been calculated

Note: if X0_real is replaced with X0_real + Xw_real, etc., this becomes bnm (as used for volumes) instead of inm (as used for panels)
"""
function calculate_ib!(harmonics, ξ0_real, ξ0_imag, η0_real, η0_imag, z0, expansion_order)
    #--- n = 0, m = 0 ---#

    i = 1
    harmonics[1,1,1] = eltype(harmonics)(0.5)
    harmonics[2,1,1] = zero(eltype(harmonics))
    i += 1

    #--- n > 0 ---#
    for n in 1:expansion_order
        i_nm1_m = i - n
        _1_m = 1
        for m in 0:n
            # get n-1 values
            j_nm1_mm1_real, j_nm1_mm1_imag, j_nm1_m_real, j_nm1_m_imag, j_nm1_mp1_real, j_nm1_mp1_imag = get_nm1(harmonics, 1, n, m, i_nm1_m, _1_m)

            # get jnm
            jnm_real = harmonics[1,2,i]
            jnm_imag = harmonics[2,2,i]

            # set coefficient
            harmonics[1,1,i] = (-(ξ0_real * j_nm1_mm1_imag + ξ0_imag * j_nm1_mm1_real) - (η0_real * j_nm1_mp1_imag + η0_imag * j_nm1_mp1_real) - z0 * j_nm1_m_real + jnm_real) / (n+2)
            harmonics[2,1,i] = ((ξ0_real * j_nm1_mm1_real - ξ0_imag * j_nm1_mm1_imag) + (η0_real * j_nm1_mp1_real - η0_imag * j_nm1_mp1_imag) - z0 * j_nm1_m_imag + jnm_imag) / (n+2)

            # recurse
            i += 1
            i_nm1_m += 1
            _1_m = -_1_m
        end
    end
end

function calculate_a!(harmonics, ξ0_real, ξ0_imag, η0_real, η0_imag, z0, expansion_order)
    #--- n = 0, m = 0 ---#

    i = 1
    harmonics[1,2,1] = eltype(harmonics)(1/6)
    harmonics[2,2,1] = zero(eltype(harmonics))
    i += 1

    #--- n > 0 ---#
    for n in 1:expansion_order
        i_nm1_m = i - n
        _1_m = 1
        for m in 0:n
            # get n-1 values
            a_nm1_mm1_real, a_nm1_mm1_imag, a_nm1_m_real, a_nm1_m_imag, a_nm1_mp1_real, a_nm1_mp1_imag = get_nm1(harmonics, 2, n, m, i_nm1_m, _1_m)

            # get bnm
            bnm_real = harmonics[1,1,i]
            bnm_imag = harmonics[2,1,i]

            # set coefficient
            harmonics[1,2,i] = (-(ξ0_real * a_nm1_mm1_imag + ξ0_imag * a_nm1_mm1_real) - (η0_real * a_nm1_mp1_imag + η0_imag * a_nm1_mp1_real) - z0 * a_nm1_m_real + bnm_real) / (n+1)
            harmonics[2,2,i] = ((ξ0_real * a_nm1_mm1_real - ξ0_imag * a_nm1_mm1_imag) + (η0_real * a_nm1_mp1_real - η0_imag * a_nm1_mp1_imag) - z0 * a_nm1_m_imag + bnm_imag) / (n+1)

            # recurse
            i += 1
            i_nm1_m += 1
            _1_m = -_1_m
        end
    end
end


"""
assumes source expansion coefficients have already been calculated
"""
function source_to_dipole!(harmonics, i_dipole, i_source, qx, qy, qz, expansion_order)

    #--- n = 0 ---#

    harmonics[1,i_dipole,1] = 0.0
    harmonics[2,i_dipole,1] = 0.0

    #--- n > 0 ---#

    i = 2
    for n in 1:expansion_order
        i_nm1_m = i - n
        _1_m = 1
        for m in 0:n
            # get n-1 values
            p_nm1_mm1_real, p_nm1_mm1_imag, p_nm1_m_real, p_nm1_m_imag, p_nm1_mp1_real, p_nm1_mp1_imag = get_nm1(harmonics, i_source, n, m, i_nm1_m, _1_m)

            # set coefficient
            harmonics[1,i_dipole,i] = -qx * 0.5 * (p_nm1_mp1_imag + p_nm1_mm1_imag) + qy * 0.5 * (p_nm1_mp1_real - p_nm1_mm1_real) - qz * p_nm1_m_real
            harmonics[2,i_dipole,i] = qx * 0.5 * (p_nm1_mp1_real + p_nm1_mm1_real) + qy * 0.5 * (p_nm1_mp1_imag - p_nm1_mm1_imag) - qz * p_nm1_m_imag

            # recurse
            i += 1
            i_nm1_m += 1
            _1_m = -_1_m
        end
    end
end

"""
my novel derivation
"""
function mirrored_source_to_vortex!(multipole_coefficients, harmonics, strength, i_source, multiplier, expansion_order)

    #--- first potential ϕ ---#

    vx, vy, vz = strength

    i = 1
    for n in 0:expansion_order
        _1_m = 1.0
        _1_np1 = 1.0/(n+1)
        for m in 0:n
            # get source coefficient values
            q_n_mm1_real, q_n_mm1_imag, q_n_m_real, q_n_m_imag, q_n_mp1_real, q_n_mp1_imag = get_n(harmonics, i_source, n, m, i, _1_m)

            # set coefficient
            nmmp1_2 = Float64(n-m+1)*0.5
            npmp1_2 = Float64(n+m+1)*0.5
            multipole_coefficients[1,1,i] += multiplier * _1_m * ((-vx * q_n_mm1_real + vy * q_n_mm1_imag) * nmmp1_2 + (vx * q_n_mp1_real + vy * q_n_mp1_imag) * npmp1_2 - vz * m * q_n_m_imag) * _1_np1
            multipole_coefficients[2,1,i] += multiplier * _1_m * ((vx * q_n_mm1_imag + vy * q_n_mm1_real) * nmmp1_2 + (-vx * q_n_mp1_imag + vy * q_n_mp1_real) * npmp1_2 - vz * m * q_n_m_real) * _1_np1

            if n == 0 && m == 0 && (abs(multiplier * _1_m * ((-vx * q_n_mm1_real + vy * q_n_mm1_imag) * nmmp1_2 + (vx * q_n_mp1_real + vy * q_n_mp1_imag) * npmp1_2 - vz * m * q_n_m_imag) * _1_np1) > 0.0 || (multiplier * _1_m * ((vx * q_n_mm1_imag + vy * q_n_mm1_real) * nmmp1_2 + (-vx * q_n_mp1_imag + vy * q_n_mp1_real) * npmp1_2 - vz * m * q_n_m_real) * _1_np1) > 0.0)
                throw("vortex element not zero for n=m=0")
            end

            # recurse
            i += 1
            _1_m = -_1_m
        end
    end

    #--- second potential χ ---#

    i = 2
    for n in 1:expansion_order
        i_nm1_m = i - n
        _1_m = 1.0
        _1_over_n = 1.0/n
        for m in 0:n
            # get n-1 values
            q_nm1_mm1_real, q_nm1_mm1_imag, q_nm1_m_real, q_nm1_m_imag, q_nm1_mp1_real, q_nm1_mp1_imag = get_nm1(harmonics, i_source, n, m, i_nm1_m, _1_m)

            # set coefficient
            multipole_coefficients[1,2,i] -= multiplier * _1_m * _1_over_n * (0.5 * (-vy * q_nm1_mm1_real - vx * q_nm1_mm1_imag + vy * q_nm1_mp1_real - vx * q_nm1_mp1_imag) - vz * q_nm1_m_real)
            multipole_coefficients[2,2,i] -= multiplier * _1_m * _1_over_n * (0.5 * (vy * q_nm1_mm1_imag - vx * q_nm1_mm1_real - vy * q_nm1_mp1_imag - vx * q_nm1_mp1_real) + vz * q_nm1_m_imag)

            # recurse
            i += 1
            i_nm1_m += 1
            _1_m = -_1_m
        end
    end
end

"""
might be slightly faster than mirrored_source_to_vortex
"""
function source_to_vortex_point!(multipole_coefficients, harmonics, strength, Δx, i_source, i_dipole, expansion_order)

    #--- first scalar potential ---#

    # analogous dipole
    qx, qy, qz = cross(Δx, strength)
    source_to_dipole!(harmonics, i_dipole, i_source, qx, qy, qz, expansion_order)

    # update coefficients
    _1_n = -1.0
    i = 2
    for n in 1:expansion_order
        _1_n_m = _1_n
        _1_np1 = -1.0/(n+1)
        for m in 0:n
            multipole_coefficients[1,1,i] += harmonics[1,i_dipole,i] * _1_n_m * _1_np1
            multipole_coefficients[2,1,i] -= harmonics[2,i_dipole,i] * _1_n_m * _1_np1 # Rnm*
            _1_n_m = -_1_n_m
            i += 1
        end
        _1_n = -_1_n
    end

    #--- second scalar potential ---#

    # analogous dipole
    qx, qy, qz = -strength
    source_to_dipole!(harmonics, i_dipole, i_source, qx, qy, qz, expansion_order)

    # update coefficients
    _1_n = -1.0
    i = 2
    for n in 1:expansion_order
        _1_n_m = _1_n
        _1_over_n = -1.0/n
        for m in 0:n
            multipole_coefficients[1,2,i] += harmonics[1,i_dipole,i] * _1_n_m * _1_over_n
            multipole_coefficients[2,2,i] -= harmonics[2,i_dipole,i] * _1_n_m * _1_over_n # Rnm*
            _1_n_m = -_1_n_m
            i += 1
        end
        _1_n = -_1_n
    end

end

#--- point elements ---#

@inline function body_to_multipole!(element::Type{<:Point}, system, branch, bodies_index, harmonics, expansion_order)

    # extract containers
    center = branch.source_center
    multipole_coefficients = branch.multipole_expansion

    # loop over bodies
    for i_body in bodies_index
        # relative body position
        x = system[i_body,POSITION]
        Δx = x - center

        # update values
        body_to_multipole_point!(element, multipole_coefficients, harmonics, Δx, system[i_body, Strength()], expansion_order)
    end

end

@inline function body_to_multipole_point!(::Type{Point{Source}}, multipole_coefficients, harmonics, Δx, strength, expansion_order)

    # transform to ξ, η, z
    ξ_real, ξ_imag, η_real, η_imag, z = xyz_to_ξηz(Δx)

    # invert the sign of the strength so v=∇ϕ instead of v=-∇ϕ
    # this will allow induced velocities to be added to those induced by vortex elements
    # we'll invert the potential again at the end
    strength = -strength

    # calculate q (regular harmonics)
    calculate_q!(harmonics, ξ_real, ξ_imag, η_real, η_imag, z, expansion_order)

    # update coefficients
    _1_n = 1.0
    i = 1
    for n in 0:expansion_order
        _1_n_m = _1_n
        for m in 0:n
            multipole_coefficients[1,1,i] += harmonics[1,1,i] * _1_n_m * strength
            multipole_coefficients[2,1,i] -= harmonics[2,1,i] * _1_n_m * strength # Rnm*
            _1_n_m = -_1_n_m
            i += 1
        end
        _1_n = -_1_n
    end
end

@inline function body_to_multipole_point!(::Type{Point{Dipole}}, multipole_coefficients, harmonics, Δx, strength, expansion_order)

    # transform to ξ, η, z
    ξ_real, ξ_imag, η_real, η_imag, z = xyz_to_ξηz(Δx)

    # calculate q
    calculate_q!(harmonics, ξ_real, ξ_imag, η_real, η_imag, z, expansion_order)

    # invert the sign of the strength so v=∇ϕ instead of v=-∇ϕ
    # this will allow induced velocities to be added to those induced by vortex elements
    # we'll invert the potential again at the end
    strength = -strength

    # convert to dipole
    i_source, i_dipole = 1, 2
    qx, qy, qz = strength
    source_to_dipole!(harmonics, i_dipole, i_source, qx, qy, qz, expansion_order)

    # update coefficients
    _1_n = 1.0
    i = 1
    for n in 0:expansion_order
        _1_n_m = _1_n
        for m in 0:n
            multipole_coefficients[1,1,i] += harmonics[1,i_dipole,i] * _1_n_m
            multipole_coefficients[2,1,i] -= harmonics[2,i_dipole,i] * _1_n_m # Rnm*
            _1_n_m = -_1_n_m
            i += 1
        end
        _1_n = -_1_n
    end

end

@inline function body_to_multipole_point!(::Type{Point{Vortex}}, multipole_coefficients, harmonics, Δx, strength, expansion_order)

    ## transform to ξ, η, z
    #ξ_real, ξ_imag, η_real, η_imag, z = xyz_to_ξηz(Δx)

    ## calculate q (regular harmonics for source point)
    #calculate_q!(harmonics, ξ_real, ξ_imag, η_real, η_imag, z, expansion_order)

    ## convert to vortex
    #i_source, i_dipole = 1, 2
    #source_to_vortex_point!(multipole_coefficients, harmonics, strength, Δx, i_source, i_dipole, expansion_order)

    # not sure if above or below is more efficient, but below generalizes to filaments/panels/volumes

    # transform mirrored geometry to ξ, η, z
    ξ_real, ξ_imag, η_real, η_imag, z = xyz_to_ξηz(-Δx)

    # calculate q (regular harmonics for source point)
    calculate_q!(harmonics, ξ_real, ξ_imag, η_real, η_imag, z, expansion_order)

    # convert to vortex
    i_source = 1
    multiplier = 1.0
    mirrored_source_to_vortex!(multipole_coefficients, harmonics, strength, i_source, multiplier, expansion_order)
end

@inline function body_to_multipole_point!(::Type{Point{SourceVortex}}, multipole_coefficients, harmonics, Δx, strength, expansion_order)
    scalar_strength = strength[1]
    vector_strength = SVector{3}(strength[2], strength[3], strength[4])
    body_to_multipole_point!(Point{Source}, multipole_coefficients, harmonics, Δx, scalar_strength, expansion_order)
    body_to_multipole_point!(Point{Vortex}, multipole_coefficients, harmonics, Δx, vector_strength, expansion_order)
end

#--- filament elements ---#

@inline function body_to_multipole!(element::Type{<:Filament}, system, branch, bodies_index, harmonics, expansion_order)
    # extract containers
    center = branch.source_center
    multipole_coefficients = branch.multipole_expansion

    # loop over bodies
    for i_body in bodies_index
        # relative body position
        x1 = system[i_body,VERTEX,1]
        x2 = system[i_body,VERTEX,2]
        x0 = x1 - center
        xu = x2 - x1

        # update values
        body_to_multipole_filament!(element, multipole_coefficients, harmonics, x0, xu, system[i_body, Strength()], expansion_order)
    end
end

function body_to_multipole_filament!(::Type{Filament{Source}}, multipole_coefficients, harmonics, x0, xu, strength, expansion_order)
    # transform to ξ, η, z
    ξ0_real, ξ0_imag, η0_real, η0_imag, z0 = xyz_to_ξηz(x0)
    ξu_real, ξu_imag, ηu_real, ηu_imag, zu = xyz_to_ξηz(xu)

    ξ0ξu_real = ξ0_real + ξu_real
    ξ0ξu_imag = ξ0_imag + ξu_imag
    η0ηu_real = η0_real + ηu_real
    η0ηu_imag = η0_imag + ηu_imag
    z0zu = z0 + zu

    # calculate q
    calculate_q!(harmonics, ξ0ξu_real, ξ0ξu_imag, η0ηu_real, η0ηu_imag, z0zu, expansion_order)

    # calculate p
    calculate_pj!(harmonics, ξ0_real, ξ0_imag, η0_real, η0_imag, z0, expansion_order)

    # invert the sign of the strength so v=∇ϕ instead of v=-∇ϕ
    # this will allow induced velocities to be added to those induced by vortex elements
    # we'll invert the potential again at the end
    strength = -strength

    # jacobian of the transformation to simplex
    Δx = norm(xu)
    Jq = Δx * strength

    # update expansion coefficients
    i = 1
    _1_n = 1.0
    for n in 0:expansion_order
        _1_n_m = _1_n
        for m in 0:n
            # set coefficients
            multipole_coefficients[1,1,i] += Jq * _1_n_m * harmonics[1,2,i]
            multipole_coefficients[2,1,i] -= Jq * _1_n_m * harmonics[2,2,i]

            # recurse
            i += 1
            _1_n_m = -_1_n_m
        end
        _1_n = -_1_n
    end
end

function body_to_multipole_filament!(::Type{Filament{Dipole}}, multipole_coefficients, harmonics, x0, xu, strength, expansion_order)
    # transform to ξ, η, z
    ξ0_real, ξ0_imag, η0_real, η0_imag, z0 = xyz_to_ξηz(x0)
    ξu_real, ξu_imag, ηu_real, ηu_imag, zu = xyz_to_ξηz(xu)

    ξ0ξu_real = ξ0_real + ξu_real
    ξ0ξu_imag = ξ0_imag + ξu_imag
    η0ηu_real = η0_real + ηu_real
    η0ηu_imag = η0_imag + ηu_imag
    z0zu = z0 + zu

    # calculate q
    calculate_q!(harmonics, ξ0ξu_real, ξ0ξu_imag, η0ηu_real, η0ηu_imag, z0zu, expansion_order)

    # calculate p
    calculate_pj!(harmonics, ξ0_real, ξ0_imag, η0_real, η0_imag, z0, expansion_order)

    # invert the sign of the strength so v=∇ϕ instead of v=-∇ϕ
    # this will allow induced velocities to be added to those induced by vortex elements
    # we'll invert the potential again at the end
    strength = -strength

    # calculate s
    qx, qy, qz = strength
    i_source, i_dipole = 2, 1
    source_to_dipole!(harmonics, i_dipole, i_source, qx, qy, qz, expansion_order)

    # jacobian of the transformation to simplex
    J = norm(xu)

    # update expansion coefficients
    i = 1
    _1_n = 1.0
    for n in 0:expansion_order
        _1_n_m = _1_n
        for m in 0:n
            # set coefficients
            multipole_coefficients[1,1,i] += J * _1_n_m * harmonics[1,i_dipole,i]
            multipole_coefficients[2,1,i] -= J * _1_n_m * harmonics[2,i_dipole,i]

            # recurse
            i += 1
            _1_n_m = -_1_n_m
        end
        _1_n = -_1_n
    end
end

function body_to_multipole_filament!(::Type{Filament{Vortex}}, multipole_coefficients, harmonics, x0, xu, strength, expansion_order)

    # transform to ξ, η, z
    ξ0_real, ξ0_imag, η0_real, η0_imag, z0 = xyz_to_ξηz(-x0)
    ξu_real, ξu_imag, ηu_real, ηu_imag, zu = xyz_to_ξηz(-xu)

    ξ0ξu_real = ξ0_real + ξu_real
    ξ0ξu_imag = ξ0_imag + ξu_imag
    η0ηu_real = η0_real + ηu_real
    η0ηu_imag = η0_imag + ηu_imag
    z0zu = z0 + zu

    # calculate q
    calculate_q!(harmonics, ξ0ξu_real, ξ0ξu_imag, η0ηu_real, η0ηu_imag, z0zu, expansion_order)

    # calculate p
    calculate_pj!(harmonics, ξ0_real, ξ0_imag, η0_real, η0_imag, z0, expansion_order)

    # jacobian of the transformation to simplex
    J = norm(xu)

    # convert to vortex
    i_source = 2
    mirrored_source_to_vortex!(multipole_coefficients, harmonics, strength, i_source, J, expansion_order)

end

#--- panel elements ---#

@inline function body_to_multipole!(element::Type{<:Panel}, system, branch, bodies_index, harmonics, expansion_order)
    # extract containers
    center = branch.source_center
    multipole_coefficients = branch.multipole_expansion

    # loop over bodies
    for i_body in bodies_index
        # relative body position
        x1 = system[i_body,VERTEX,1]
        x2 = system[i_body,VERTEX,2]
        x3 = system[i_body,VERTEX,3]
        x0 = x1 - center
        xu = x2 - x1
        xv = x3 - x1

        # get normal
        normal = system[i_body,Normal()]

        # update values
        body_to_multipole_panel!(element, multipole_coefficients, harmonics, x0, xu, xv, normal, system[i_body, STRENGTH], expansion_order)
    end
end

@inline function body_to_multipole!(element::Type{Panel{SourceDipole}}, multipole_coefficients, harmonics, x0, xu, xv, normal, strength, expansion_order)
    body_to_multipole_panel!(Panel{Source}, multipole_coefficients, harmonics, x0, xu, xv, normal, strength[1], expansion_order)
    body_to_multipole_panel!(Panel{Dipole}, multipole_coefficients, harmonics, x0, xu, xv, normal, strength[2] * normal, expansion_order)
end

function body_to_multipole_panel!(::Type{Panel{Source}}, multipole_coefficients, harmonics, x0, xu, xv, normal, strength, expansion_order)
    # transform to ξ, η, z
    ξ0_real, ξ0_imag, η0_real, η0_imag, z0 = xyz_to_ξηz(x0)
    ξu_real, ξu_imag, ηu_real, ηu_imag, zu = xyz_to_ξηz(xu)
    ξv_real, ξv_imag, ηv_real, ηv_imag, zv = xyz_to_ξηz(xv)

    # calculate q
    ξ0ξu_real = ξ0_real + ξu_real
    ξ0ξu_imag = ξ0_imag + ξu_imag
    η0ηu_real = η0_real + ηu_real
    η0ηu_imag = η0_imag + ηu_imag
    z0zu = z0 + zu

    calculate_q!(harmonics, ξ0ξu_real, ξ0ξu_imag, η0ηu_real, η0ηu_imag, z0zu, expansion_order)

    # calculate j
    ξ0ξv_real = ξ0_real + ξv_real
    ξ0ξv_imag = ξ0_imag + ξv_imag
    η0ηv_real = η0_real + ηv_real
    η0ηv_imag = η0_imag + ηv_imag
    z0zv = z0 + zv

    calculate_pj!(harmonics, ξ0ξv_real, ξ0ξv_imag, η0ηv_real, η0ηv_imag, z0zv, expansion_order)

    # calculate i
    calculate_ib!(harmonics, ξ0_real, ξ0_imag, η0_real, η0_imag, z0, expansion_order)

    # invert the sign of the strength so v=∇ϕ instead of v=-∇ϕ
    # this will allow induced velocities to be added to those induced by vortex elements
    # we'll invert the potential again at the end
    strength = -strength

    # jacobian of the transformation to simplex
    Δx = norm(cross(xu,xv))
    Jq = Δx * strength

    # update expansion coefficients
    i = 1
    _1_n = 1.0
    for n in 0:expansion_order
        _1_n_m = _1_n
        for m in 0:n
            # set coefficients
            multipole_coefficients[1,1,i] += Jq * _1_n_m * harmonics[1,1,i]
            multipole_coefficients[2,1,i] -= Jq * _1_n_m * harmonics[2,1,i]

            # recurse
            i += 1
            _1_n_m = -_1_n_m
        end
        _1_n = -_1_n
    end
end

function body_to_multipole_panel!(::Type{Panel{Dipole}}, multipole_coefficients, harmonics, x0, xu, xv, normal, strength, expansion_order)
    # transform to ξ, η, z
    ξ0_real, ξ0_imag, η0_real, η0_imag, z0 = xyz_to_ξηz(x0)
    ξu_real, ξu_imag, ηu_real, ηu_imag, zu = xyz_to_ξηz(xu)
    ξv_real, ξv_imag, ηv_real, ηv_imag, zv = xyz_to_ξηz(xv)

    # calculate q
    ξ0ξu_real = ξ0_real + ξu_real
    ξ0ξu_imag = ξ0_imag + ξu_imag
    η0ηu_real = η0_real + ηu_real
    η0ηu_imag = η0_imag + ηu_imag
    z0zu = z0 + zu

    calculate_q!(harmonics, ξ0ξu_real, ξ0ξu_imag, η0ηu_real, η0ηu_imag, z0zu, expansion_order)

    # calculate j
    ξ0ξv_real = ξ0_real + ξv_real
    ξ0ξv_imag = ξ0_imag + ξv_imag
    η0ηv_real = η0_real + ηv_real
    η0ηv_imag = η0_imag + ηv_imag
    z0zv = z0 + zv

    calculate_pj!(harmonics, ξ0ξv_real, ξ0ξv_imag, η0ηv_real, η0ηv_imag, z0zv, expansion_order)

    # calculate i
    calculate_ib!(harmonics, ξ0_real, ξ0_imag, η0_real, η0_imag, z0, expansion_order)

    # invert the sign of the strength so v=∇ϕ instead of v=-∇ϕ
    # this will allow induced velocities to be added to those induced by vortex elements
    # we'll invert the potential again at the end
    qx, qy, qz = -strength

    # convert to dipole coefficients
    i_source, i_dipole = 1, 2
    source_to_dipole!(harmonics, i_dipole, i_source, qx, qy, qz, expansion_order)

    # jacobian of the transformation to simplex
    J = norm(cross(xu,xv))

    # update expansion coefficients
    i = 1
    _1_n = 1.0
    for n in 0:expansion_order
        _1_n_m = _1_n
        for m in 0:n
            # set coefficients
            multipole_coefficients[1,1,i] += J * _1_n_m * harmonics[1,i_dipole,i]
            multipole_coefficients[2,1,i] -= J * _1_n_m * harmonics[2,i_dipole,i]

            # recurse
            i += 1
            _1_n_m = -_1_n_m
        end
        _1_n = -_1_n
    end
end

function body_to_multipole_panel!(::Type{Panel{Vortex}}, multipole_coefficients, harmonics, x0, xu, xv, normal, strength, expansion_order)
    # transform to ξ, η, z
    ξ0_real, ξ0_imag, η0_real, η0_imag, z0 = xyz_to_ξηz(-x0)
    ξu_real, ξu_imag, ηu_real, ηu_imag, zu = xyz_to_ξηz(-xu)

    # calculate q
    ξ0ξu_real = ξ0_real + ξu_real
    ξ0ξu_imag = ξ0_imag + ξu_imag
    η0ηu_real = η0_real + ηu_real
    η0ηu_imag = η0_imag + ηu_imag
    z0zu = z0 + zu

    calculate_q!(harmonics, ξ0ξu_real, ξ0ξu_imag, η0ηu_real, η0ηu_imag, z0zu, expansion_order)

    # calculate j
    ξv_real, ξv_imag, ηv_real, ηv_imag, zv = xyz_to_ξηz(-xv)
    ξ0ξv_real = ξ0_real + ξv_real
    ξ0ξv_imag = ξ0_imag + ξv_imag
    η0ηv_real = η0_real + ηv_real
    η0ηv_imag = η0_imag + ηv_imag
    z0zv = z0 + zv

    calculate_pj!(harmonics, ξ0ξv_real, ξ0ξv_imag, η0ηv_real, η0ηv_imag, z0zv, expansion_order)

    # calculate i
    calculate_ib!(harmonics, ξ0_real, ξ0_imag, η0_real, η0_imag, z0, expansion_order)

    # TODO
    # jacobian of the transformation to simplex
    J = norm(cross(xu,xv))

    # convert to vortex
    i_source = 1
    mirrored_source_to_vortex!(multipole_coefficients, harmonics, strength, i_source, J, expansion_order)

end

#=
#--- panel elements ---#

function B2M!(element::Type{Panel{NS,<:Any}}, system, branch, bodies_index, harmonics, expansion_order::Val{P}) where {NS,P}

    # identify containers
    if P > 2
        harmonics .= zero(TF)
        qnm_prev = view(harmonics,1:2,1:P+2)
        jnm_prev = view(harmonics,1:2,P+3:P<<1+4)
        inm_prev = view(harmonics,1:2,P<<1+5:3*P+6)
    else
        qnm_prev = zeros(TF, 2, P+2)
        jnm_prev = zeros(TF, 2, P+2)
        inm_prev = zeros(TF, 2, P+2)
    end

    # loop over panels
    for i_body in bodies_index
        strength = system[i_body, SCALAR_STRENGTH]
        R0 = system[i_body, VERTEX, 1] - branch.center
        normal = system[i_body,NORMAL]
        for i_side in 2:NS-1
            Ru = system[i_body,VERTEX,i_side] - R0_global
            Rv = system[i_body,VERTEX,i_side + 1] - R0_global
            if compute_normal
                normal = cross(Rv, Ru)
                normal /= norm(normal) * (-1)^invert_normal
            end
            B2M_panel!(element, branch.multipole_expansion, qnm_prev, jnm_prev, inm_prev, R0, Ru, Rv, normal, system, i_body, expansion_order)
        end
    end

end

@inline function B2M_panel!(element::Type{Panel{<:Any,<:Union{ConstantSource{scale}, ConstantNormalDipole{scale}, ConstantSourceNormalDipole{scale}}}}, multipole_expansion, qnm_prev, jnm_prev, inm_prev, R0, Ru, Rv, normal, system, i_body, expansion_order) where scale
    strength = system[i_body, SCALAR_STRENGTH] * scale
    coefficients = view(multipole_expansion, :, 1, :)
    _B2M!_panel(coefficients, qnm_prev, jnm_prev, inm_prev, R0, Ru, Rv, strength, normal, expansion_order, element)
end

@inline function B2M_panel!(element::Type{Panel{<:Any,Vortex{scale}}}, multipole_expansion, qnm_prev, jnm_prev, inm_prev, R0, Ru, Rv, normal, system, i_body, expansion_order) where scale
    qx, qy, qz = system[i_body, VECTOR_STRENGTH]

    # x component
    coefficients = view(multipole_expansion, :, 2, :)
    _B2M!_panel(coefficients, qnm_prev, jnm_prev, inm_prev, R0, Ru, Rv, qx * scale, normal, expansion_order, Panel{0,Source{scale}})

    # y component
    coefficients = view(multipole_expansion, :, 3, :)
    _B2M!_panel(coefficients, qnm_prev, jnm_prev, inm_prev, R0, Ru, Rv, qy * scale, normal, expansion_order, Panel{0,Source{scale}})

    # z component
    coefficients = view(multipole_expansion, :, 4, :)
    _B2M!_panel(coefficients, qnm_prev, jnm_prev, inm_prev, R0, Ru, Rv, qz * scale, normal, expansion_order, Panel{0,Source{scale}})
end

@inline function update_multipole_expansion_panel!(coefficients, i_compressed, strength_J, inm_real, inm_imag, iam_real, iam_imag, inm_prev_mm1_real, inm_prev_mm1_imag, inm_prev_m_real, inm_prev_m_imag, inm_prev_mp1_real, inm_prev_mp1_imag, n, panel_type::Type{<:Panel{<:Any,<:ConstantSource}}, negative_1_n)
    coefficients[1,i_compressed] += strength_J * complex_multiply_real(inm_real, inm_imag, iam_real, iam_imag)
    coefficients[2,i_compressed] -= strength_J * complex_multiply_imag(inm_real, inm_imag, iam_real, iam_imag)
end

@inline function update_multipole_expansion_panel!(coefficients, i_compressed, strength_J, inm_real, inm_imag, iam_real, iam_imag, inm_prev_mm1_real, inm_prev_mm1_imag, inm_prev_m_real, inm_prev_m_imag, inm_prev_mp1_real, inm_prev_mp1_imag, n, panel_type::Type{<:Panel{<:Any,<:ConstantNormalDipole}}, negative_1_n)
    sum_iterm_real = (inm_prev_mp1_real + inm_prev_mm1_real)/2
    sum_iterm_imag = (inm_prev_mp1_imag + inm_prev_mm1_imag)/2
    x1_real, x1_imag = complex_multiply(inm_prev_mm1_real + inm_prev_mp1_real, inm_prev_mm1_imag + inm_prev_mp1_imag, 0.0, n[1]/2)
    x2_real, x2_imag = n[2] * (inm_prev_mp1_real - inm_prev_mm1_real)/2, n[2] * (inm_prev_mp1_imag - inm_prev_mm1_imag)/2
    x3_real, x3_imag = n[3] * inm_prev_m_real, n[3] * inm_prev_m_imag
    lnm_real, lnm_imag = x1_real + x2_real - x3_real, x1_imag + x2_imag - x3_imag
    # coeff_gumerov = strength_J * (lnm_real - im*lnm_imag) * (-1)^m / (1.0im)^m
    Mnm_real, Mnm_imag = complex_multiply(lnm_real, -lnm_imag, iam_real, iam_imag, negative_1_n*strength_J, 0.0)

    coefficients[1,i_compressed] += Mnm_real
    coefficients[2,i_compressed] += Mnm_imag
end

@inline function update_multipole_expansion_panel!(coefficients, i_compressed, strength_J, inm_real, inm_imag, iam_real, iam_imag, inm_prev_mm1_real, inm_prev_mm1_imag, inm_prev_m_real, inm_prev_m_imag, inm_prev_mp1_real, inm_prev_mp1_imag, n, panel_type::Type{<:Panel{<:Any,<:ConstantSourceNormalDipole}}, negative_1_n)
    update_multipole_expansion_panel!(coefficients, i_compressed, strength_J[1], inm_real, inm_imag, iam_real, iam_imag, inm_prev_mm1_real, inm_prev_mm1_imag, inm_prev_m_real, inm_prev_m_imag, inm_prev_mp1_real, inm_prev_mp1_imag, n, Panel{<:Any,ConstantSource{scale}}, negative_1_n)
    update_multipole_expansion_panel!(coefficients, i_compressed, strength_J[2], inm_real, inm_imag, iam_real, iam_imag, inm_prev_mm1_real, inm_prev_mm1_imag, inm_prev_m_real, inm_prev_m_imag, inm_prev_mp1_real, inm_prev_mp1_imag, n, Panel{<:Any,ConstantNormalDipole{scale}}, negative_1_n)
end

@inline function update_multipole_expansion_panel!(coefficients, i_compressed, strength_J, inm_real, inm_imag, panel_type::Type{<:Panel{<:Any,<:ConstantSource}})
    coefficients[1,i_compressed] += strength_J * inm_real
    coefficients[2,i_compressed] -= strength_J * inm_imag
end

@inline function update_multipole_expansion_panel!(coefficients, i_compressed, strength_J, inm_real, inm_imag, panel_type::Type{<:Panel{<:Any,<:ConstantNormalDipole}})
    return nothing
end

@inline function update_multipole_expansion_panel!(coefficients, i_compressed, strength_J, inm_real, inm_imag, panel_type::Type{Panel{<:Any,ConstantSourceNormalDipole{scale}}}) where scale
    update_multipole_expansion_panel!(coefficients, i_compressed, strength_J[1], inm_real, inm_imag, Panel{0,ConstantSource{scale}})
    update_multipole_expansion_panel!(coefficients, i_compressed, strength_J[2], inm_real, inm_imag, Panel{0,ConstantNormalDipole{scale}})
end

@inline function _B2M!_panel(multipole_expansion, qnm_prev, jnm_prev, inm_prev, R0, Ru, Rv, strength, normal, expansion_order::Val{P}, panel_type::Type{<:Panel}) where P
    J = norm(cross(Ru,Rv))

    # transform into xi, eta, z
    xi0_real = R0[1]/2
    xi0_imag = R0[2]/2
    xiu_real = Ru[1]/2
    xiu_imag = Ru[2]/2
    xiv_real = Rv[1]/2
    xiv_imag = Rv[2]/2
    eta0_real = R0[1]/2
    eta0_imag = -R0[2]/2
    etau_real = Ru[1]/2
    etau_imag = -Ru[2]/2
    etav_real = Rv[1]/2
    etav_imag = -Rv[2]/2
    z0 = R0[3]
    zu = Ru[3]
    zv = Rv[3]

    # initialize recurrent values
    qnm_real = 1.0
    qnm_imag = 0.0
    qnm_m1_real = 0.0
    qnm_m1_imag = 0.0
    jnm_real = 1.0
    jnm_imag = 0.0
    jnm_m1_real = 0.0
    jnm_m1_imag = 0.0
    inm_real = 0.5
    inm_imag = 0.0
    inm_m1_real = 0.0
    inm_m1_imag = 0.0
    strength_J = J * strength

    # expansion coefficient for l=0, m=0
    # coefficients[1,1] += real(inm * strength_J)
    # coefficients[2,1] += imag(inm * strength_J) # always zero
    # coefficients[1,1] += inm_real * strength_J
    # coefficients[2,1] += inm_imag * strength_J # always zero

    update_multipole_expansion_panel!(multipole_expansion, 1, strength_J, inm_real, inm_imag, panel_type)

    qnm_prev[1,1] = qnm_real
    qnm_prev[2,1] = qnm_imag
    jnm_prev[1,1] = jnm_real
    jnm_prev[2,1] = jnm_imag
    inm_prev[1,1] = inm_real
    inm_prev[2,1] = inm_imag

    # i^(-1)
    one_over_i_real = 0
    one_over_i_imag = -1

    # recurse
    for n in 1:P
        # m=0
        iam_real = 1
        iam_imag = 0

        # qnm = (im*(xi0+xiu)*-conj(qnm_prev[2]) + im*(eta0+etau)*qnm_prev[2] - (z0+zu)*qnm_prev[1])/n
        x1_real, x1_imag = complex_multiply(xi0_real+xiu_real, xi0_imag+xiu_imag, -qnm_prev[1,2], qnm_prev[2,2], 0, 1)
        x2_real, x2_imag = complex_multiply(eta0_real+etau_real, eta0_imag+etau_imag, qnm_prev[1,2], qnm_prev[2,2], 0, 1)
        x3_real = (z0+zu)*qnm_prev[1,1]
        x3_imag = (z0+zu)*qnm_prev[2,1]
        qnm_real = (x1_real + x2_real - x3_real) / n
        qnm_imag = (x1_imag + x2_imag - x3_imag) / n

        # jnm = (im*(xi0+xiv)*-conj(jnm_prev[2]) + im*(eta0+etav)*jnm_prev[2] - (z0+zv)*jnm_prev[1] + qnm)/(n+1)
        x1_real, x1_imag = complex_multiply(xi0_real+xiv_real, xi0_imag+xiv_imag, -jnm_prev[1,2], jnm_prev[2,2], 0, 1)
        x2_real, x2_imag = complex_multiply(eta0_real+etav_real, eta0_imag+etav_imag, jnm_prev[1,2], jnm_prev[2,2], 0, 1)
        x3_real = (z0+zv) * jnm_prev[1,1]
        x3_imag = (z0+zv) * jnm_prev[2,1]
        jnm_real = (x1_real + x2_real - x3_real + qnm_real) / (n+1)
        jnm_imag = (x1_imag + x2_imag - x3_imag + qnm_imag) / (n+1)

        # inm = (im*xi0*-conj(inm_prev[2]) + im*eta0*inm_prev[2] - z0*inm_prev[1] + jnm)/(n+2)
        x1_real, x1_imag = complex_multiply(xi0_real, xi0_imag, -inm_prev[1,2], inm_prev[2,2], 0, 1)
        x2_real, x2_imag = complex_multiply(eta0_real, eta0_imag, inm_prev[1,2], inm_prev[2,2], 0, 1)
        x3_real, x3_imag = z0 * inm_prev[1,1], z0 * inm_prev[2,1]
        inm_real = (x1_real + x2_real - x3_real + jnm_real) / (n+2)
        inm_imag = (x1_imag + x2_imag - x3_imag + jnm_imag) / (n+2)

        i_compressed = 1 + (n * (n + 1)) >> 1
        # coefficients[1,1,i_compressed] += real(conj(inm * strength_J))
        # coefficients[2,1,i_compressed] += imag(conj(inm * strength_J))
        # coefficients[1,1,i_compressed] += inm_real * strength_J
        # coefficients[2,1,i_compressed] += -inm_imag * strength_J

        # update_multipole_expansion_panel!(multipole_expansion, i_compressed, strength_J, inm_real, inm_imag, panel_type)
        update_multipole_expansion_panel!(multipole_expansion, i_compressed, strength_J, inm_real, inm_imag, iam_real, iam_imag, -inm_prev[1,2], inm_prev[2,2], inm_prev[1,1], inm_prev[2,1], inm_prev[1,2], inm_prev[2,2], normal, panel_type, 1)

        qnm_m1_real = qnm_real
        qnm_m1_imag = qnm_imag
        jnm_m1_real = jnm_real
        jnm_m1_imag = jnm_imag
        inm_m1_real = inm_real
        inm_m1_imag = inm_imag

        # (-1)^m
        negative_1_m = -1

        for m in 1:n
            iam_real, iam_imag = complex_multiply(iam_real,iam_imag,one_over_i_real,one_over_i_imag)

            # qnm = (im*(xi0+xiu)*qnm_prev[m] + im*(eta0+etau)*qnm_prev[m+2] - (z0+zu)*qnm_prev[m+1])/n
            x1_real, x1_imag = complex_multiply(xi0_real+xiu_real, xi0_imag+xiu_imag, qnm_prev[1,m], qnm_prev[2,m], 0, 1)
            x2_real, x2_imag = complex_multiply(eta0_real+etau_real, eta0_imag+etau_imag, qnm_prev[1,m+2], qnm_prev[2,m+2], 0, 1)
            x3_real, x3_imag = (z0+zu)*qnm_prev[1,m+1], (z0+zu)*qnm_prev[2,m+1]
            qnm_real = (x1_real + x2_real - x3_real) / n
            qnm_imag = (x1_imag + x2_imag - x3_imag) / n

            # jnm = (im*(xi0+xiv)*jnm_prev[m] + im*(eta0+etav)*jnm_prev[m+2] - (z0+zv)*jnm_prev[m+1] + qnm)/(n+1)
            x1_real, x1_imag = complex_multiply(xi0_real+xiv_real, xi0_imag+xiv_imag, jnm_prev[1,m], jnm_prev[2,m], 0, 1)
            x2_real, x2_imag = complex_multiply(eta0_real+etav_real, eta0_imag+etav_imag, jnm_prev[1,m+2], jnm_prev[2,m+2], 0, 1)
            x3_real, x3_imag = (z0+zv)*jnm_prev[1,m+1], (z0+zv)*jnm_prev[2,m+1]
            jnm_real = (x1_real + x2_real - x3_real + qnm_real) / (n+1)
            jnm_imag = (x1_imag + x2_imag - x3_imag + qnm_imag) / (n+1)

            # inm = (im*xi0*inm_prev[m] + im*eta0*inm_prev[m+2] - z0*inm_prev[m+1] + jnm)/(n+2)
            x1_real, x1_imag = complex_multiply(xi0_real, xi0_imag, inm_prev[1,m], inm_prev[2,m], 0, 1)
            x2_real, x2_imag = complex_multiply(eta0_real, eta0_imag, inm_prev[1,m+2], inm_prev[2,m+2], 0, 1)
            x3_real, x3_imag = z0*inm_prev[1,m+1], z0*inm_prev[2,m+1]
            inm_real = (x1_real + x2_real - x3_real + jnm_real) / (n+2)
            inm_imag = (x1_imag + x2_imag - x3_imag + jnm_imag) / (n+2)

            i_compressed += 1
            # coefficients[1,i_compressed] += real(conj(inm * strength_J * iam))
            # coefficients[2,i_compressed] += imag(conj(inm * strength_J * iam))

            # coefficients[1,i_compressed] += strength_J * complex_multiply_real(inm_real, inm_imag, iam_real, iam_imag)
            # coefficients[2,i_compressed] -= strength_J * complex_multiply_imag(inm_real, inm_imag, iam_real, iam_imag)
            update_multipole_expansion_panel!(multipole_expansion, i_compressed, strength_J, inm_real, inm_imag, iam_real, iam_imag, inm_prev[1,m], inm_prev[2,m], inm_prev[1,m+1], inm_prev[2,m+1], inm_prev[1,m+2], inm_prev[2,m+2], normal, panel_type, negative_1_m)

            qnm_prev[1,m] = qnm_m1_real
            qnm_prev[2,m] = qnm_m1_imag
            jnm_prev[1,m] = jnm_m1_real
            jnm_prev[2,m] = jnm_m1_imag
            inm_prev[1,m] = inm_m1_real
            inm_prev[2,m] = inm_m1_imag
            qnm_m1_real = qnm_real
            qnm_m1_imag = qnm_imag
            jnm_m1_real = jnm_real
            jnm_m1_imag = jnm_imag
            inm_m1_real = inm_real
            inm_m1_imag = inm_imag

            # update (-1)^m
            negative_1_m *= -1
        end
        qnm_prev[1,n+1] = qnm_real
        qnm_prev[2,n+1] = qnm_imag
        jnm_prev[1,n+1] = jnm_real
        jnm_prev[2,n+1] = jnm_imag
        inm_prev[1,n+1] = inm_real
        inm_prev[2,n+1] = inm_imag
    end
end

# function B2M!_quadpanel(system, branch, bodies_index, harmonics::AbstractArray{TF}, expansion_order::Val{P}, panel_type::AbstractKernel{scale}; compute_normal=false, invert_normal=false) where {TF,P,scale}
#     if P > 2
#         harmonics .= zero(TF)
#         qnm_prev = view(harmonics,1:2,1:P+2)
#         jnm_prev = view(harmonics,1:2,P+3:P<<1+4)
#         inm_prev = view(harmonics,1:2,P<<1+5:3*P+6)
#     else
#         qnm_prev = zeros(TF, 2, P+2)
#         jnm_prev = zeros(TF, 2, P+2)
#         inm_prev = zeros(TF, 2, P+2)
#     end
#     for i_body in bodies_index
#         strength = system[i_body, SCALAR_STRENGTH]
#         R0_global = system[i_body, VERTEX, 1]
#         R0 = R0_global - branch.center
#         for i_side in 2:3
#             Ru = system[i_body,VERTEX,i_side] - R0_global
#             Rv = system[i_body,VERTEX,i_side + 1] - R0_global
#             if compute_normal
#                 normal = cross(Rv, Ru)
#                 normal /= norm(normal) * (-1)^invert_normal
#             else
#                 normal = system[i_body,NORMAL]
#             end
#             _B2M!_panel(branch.multipole_expansion, qnm_prev, jnm_prev, inm_prev, R0, Ru, Rv, strength, normal, expansion_order, panel_type)
#         end
#     end
# end
#
# function B2M!_tripanel(system, branch, bodies_index, harmonics::AbstractArray{TF}, expansion_order::Val{P}, panel_type::AbstractKernel) where {TF,P}
#     if P > 2
#         harmonics .= zero(TF)
#         qnm_prev = view(harmonics,1:2,1:P+2)
#         jnm_prev = view(harmonics,1:2,P+3:P<<1+4)
#         inm_prev = view(harmonics,1:2,P<<1+5:3*P+6)
#     else
#         qnm_prev = zeros(TF, 2, P+2)
#         jnm_prev = zeros(TF, 2, P+2)
#         inm_prev = zeros(TF, 2, P+2)
#     end
#     for i_body in bodies_index
#         strength = system[i_body, SCALAR_STRENGTH]
#         R0_global = system[i_body, VERTEX, 1]
#         R0 = R0_global - branch.center
#         Ru = system[i_body,VERTEX,2] - R0_global
#         Rv = system[i_body,VERTEX,3] - R0_global
#         normal = system[i_body,NORMAL]
#         _B2M!_panel(branch.multipole_expansion, qnm_prev, jnm_prev, inm_prev, R0, Ru, Rv, strength, normal, expansion_order, panel_type)
#     end
# end
=#
