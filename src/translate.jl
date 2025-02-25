#------- MULTIPOLE TO MULTIPOLE -------#

"""
Overwrites translated_weights
"""
function translate_multipole_z!(translated_weights, source_weights, t, expansion_order, ::Val{LH}) where LH
#function translate_multipole_z!(translated_weights, source_weights, t, expansion_order::Val{P}, ::Val{LH}) where {P,LH}
    _t = -t
    i = 1
    for n in 0:expansion_order
    	for m in 0:n
    		# preallocate recursive variables
            val1_real = zero(eltype(translated_weights))
            val1_imag = zero(eltype(translated_weights))
            if LH
                val2_real = zero(eltype(translated_weights))
                val2_imag = zero(eltype(translated_weights))
            end
    		n_np! = 1.0
    		n_np = 1
    		_t_n_np = one(t)

            for np in n:-1:m
                tmp = _t_n_np / n_np!
                i_weight = harmonic_index(np,m)
                val1_real += tmp * source_weights[1,1,i_weight]
                val1_imag += tmp * source_weights[2,1,i_weight]
                if LH
                    val2_real += tmp * source_weights[1,2,i_weight]
                    val2_imag += tmp * source_weights[2,2,i_weight]
                end
    			_t_n_np *= _t
    			n_np! *= n_np
    			n_np += 1
    		end

            # set translated coefficient
            translated_weights[1,1,i] = val1_real
            translated_weights[2,1,i] = val1_imag

            if LH
                translated_weights[1,2,i] = val2_real
                translated_weights[2,2,i] = val2_imag
            end

            # increment index
            i += 1
    	end
    end
end

function transform_lamb_helmholtz_multipole!(multipole_expansion, r, expansion_order)
    i_P_m = harmonic_index(expansion_order, 0)
    for m in 0:expansion_order
        # declare recursive variable for inner loop over n
        χ̂_nm1_real = multipole_expansion[1,2,i_P_m]
        χ̂_nm1_imag = multipole_expansion[2,2,i_P_m]

        # declare index
        i_n_m = i_P_m

        for n in expansion_order:-1:max(m,1)
            # transform ϕ potential
            ϕ̂_real = multipole_expansion[1,1,i_n_m]
            ϕ̂_imag = multipole_expansion[2,1,i_n_m]

            # recursively access χ̂ to reduce cost
            χ̂_real = χ̂_nm1_real
            χ̂_imag = χ̂_nm1_imag

            # ϕ̃ = ϕ̂ - im * r * m/(n+1) * χ̂
            r_m_np1 = r * m / (n+1)
            multipole_expansion[1,1,i_n_m] = ϕ̂_real + r_m_np1 * χ̂_imag
            multipole_expansion[2,1,i_n_m] = ϕ̂_imag - r_m_np1 * χ̂_real

            # recurse χ̂_nm1
            if m < n
                χ̂_nm1_real = multipole_expansion[1,2,i_n_m-n]
                χ̂_nm1_imag = multipole_expansion[2,2,i_n_m-n]
            else # χ_{n-1}^m doesn't exist
                χ̂_nm1_real = zero(eltype(multipole_expansion))
                χ̂_nm1_imag = zero(eltype(multipole_expansion))
            end

            # χ̃ = χ̂ + r/n * χ̂_nm1
            r_n = r / n
            multipole_expansion[1,2,i_n_m] = χ̂_real + r_n * χ̂_nm1_real
            multipole_expansion[2,2,i_n_m] = χ̂_imag + r_n * χ̂_nm1_imag

            # update index
            i_n_m -= n

        end
        i_P_m += 1
    end
end

function transform_lamb_helmholtz_local!(local_expansion, r, expansion_order)
    i_m_m = 1
    for m in 0:expansion_order
        # declare recursive variable for inner loop over n
        χ̂_np1_real = local_expansion[1,2,i_m_m]
        χ̂_np1_imag = local_expansion[2,2,i_m_m]

        # declare index for fast access
        i_n_m = i_m_m

        for n in m:expansion_order # skip n=0
            # recursively access χ̂ to reduce cost
            χ̂_real = χ̂_np1_real
            χ̂_imag = χ̂_np1_imag

            if n > 0
                # transform ϕ potential
                ϕ̂_real = local_expansion[1,1,i_n_m]
                ϕ̂_imag = local_expansion[2,1,i_n_m]

                # ϕ̃ = ϕ̂ + im * r * m/n * χ̂
                r_m_n = r * m / n
                local_expansion[1,1,i_n_m] = ϕ̂_real - r_m_n * χ̂_imag
                local_expansion[2,1,i_n_m] = ϕ̂_imag + r_m_n * χ̂_real
            end

            if n < expansion_order # χ_{n+1}^m = 0 when n+1>P
                # recurse χ̂_np1
                χ̂_np1_real = local_expansion[1,2,i_n_m+n+1]
                χ̂_np1_imag = local_expansion[2,2,i_n_m+n+1]

                # χ̃ = χ̂ - r/(n+1) * χ̂_nm1
                r_np1 = r / (n+1)
                local_expansion[1,2,i_n_m] = χ̂_real - r_np1 * χ̂_np1_real
                local_expansion[2,2,i_n_m] = χ̂_imag - r_np1 * χ̂_np1_imag
            end

            # update index
            i_n_m += n + 1

        end

        i_m_m += m + 2
    end
end

function multipole_to_multipole!(target_branch, source_branch, weights_tmp_1, weights_tmp_2, Ts, eimϕs, ζs_mag, Hs_π2, expansion_order, lamb_helmholtz::Val{LH}) where LH
    # extract containers
    source_weights = source_branch.multipole_expansion

    # translation vector
    Δx = target_branch.source_center - source_branch.source_center
    r, θ, ϕ = cartesian_to_spherical(Δx)

    #--- rotate coordinate system ---#

    # rotate about z axis
    rotate_z!(weights_tmp_1, source_weights, eimϕs, ϕ, expansion_order, lamb_helmholtz)

    # rotate about y axis
    # NOTE: the method used here results in an additional rotatation of π about the new z axis
    rotate_multipole_y!(weights_tmp_2, weights_tmp_1, Ts, Hs_π2, ζs_mag, θ, expansion_order, lamb_helmholtz)

    #--- translate along new z axis ---#

    translate_multipole_z!(weights_tmp_1, weights_tmp_2, r, expansion_order, lamb_helmholtz)

    #--- transform Lamb-Helmholtz decomposition for the new center ---#

    LH && transform_lamb_helmholtz_multipole!(weights_tmp_1, r, expansion_order)

    #--- back rotate coordinate system ---#

    # back rotate about y axis
    back_rotate_multipole_y!(weights_tmp_2, weights_tmp_1, Ts, ζs_mag, expansion_order, lamb_helmholtz)

    # back rotate about z axis and accumulate on target branch
    back_rotate_z!(target_branch.multipole_expansion, weights_tmp_2, eimϕs, expansion_order, lamb_helmholtz)

end

#------- MULTIPOLE TO LOCAL -------#

"""
Overwrites translated_weights
"""
function translate_multipole_to_local_z!(translated_weights, source_weights, t, expansion_order, ::Val{LH}) where LH
    one_over_t = 1/t
    n!_t_np1 = one_over_t
    i = 1
    for n in 0:expansion_order
        n_m!_t_nmp1 = n!_t_np1
        n_m = n
    	for m in 0:n
            # inner summation
            val1_real = zero(eltype(translated_weights))
            val1_imag = zero(eltype(translated_weights))

            if LH
                val2_real = zero(eltype(translated_weights))
                val2_imag = zero(eltype(translated_weights))
            end
    		n_np! = n_m!_t_nmp1
            n_np = n_m
    		for np in m:expansion_order
                i_weight = harmonic_index(np,m)
                val1_real += n_np! * source_weights[1,1,i_weight]
                val1_imag += n_np! * source_weights[2,1,i_weight]

                if LH
                    val2_real += n_np! * source_weights[1,2,i_weight]
                    val2_imag += n_np! * source_weights[2,2,i_weight]
                end

                n_np += 1
                n_np! *= n_np * one_over_t
    		end

            # set translated coefficient
            translated_weights[1,1,i] = val1_real
            translated_weights[2,1,i] = val1_imag

            if LH
                translated_weights[1,2,i] = val2_real
                translated_weights[2,2,i] = val2_imag
            end

            # increment index
            i += 1

            # recurse values
            n_m += 1
            n_m!_t_nmp1 *= n_m * one_over_t
    	end
        n!_t_np1 *= (n+1)*one_over_t
    end
end

"""
This function should receive n!_t_np1 = 1/t when n=0

Calculates ϕn0, ϕn1, χn1 for P=n
"""
function translate_multipole_to_local_z_m01_n(source_weights, t, one_over_t, ::Val{LH}, n!_t_np1, n) where LH
    # initialize values
    n_m!_t_nmp1 = n!_t_np1
    n_m = n

    #--- m = 0 ---#

    # inner summation
    ϕn0_real = zero(eltype(source_weights))

    n_np! = n_m!_t_nmp1
    n_np = n_m
    for np in 0:n
        i_weight = harmonic_index(np,0)
        ϕn0_real += n_np! * source_weights[1,1,i_weight]

        n_np += 1
        n_np! *= n_np * one_over_t
    end

    # recurse values
    n_m += 1
    n_m!_t_nmp1 *= n_m * one_over_t

    #--- m = 1 ---#

    # inner summation
    ϕn1_real = zero(eltype(source_weights))
    ϕn1_imag = zero(eltype(source_weights))
    χn1_real = zero(eltype(source_weights))
    χn1_imag = zero(eltype(source_weights))

    n_np! = n_m!_t_nmp1
    n_np = n_m
    for np in 1:n
        i_weight = harmonic_index(np,1)
        ϕn1_real += n_np! * source_weights[1,1,i_weight]
        ϕn1_imag += n_np! * source_weights[2,1,i_weight]

        if LH
            χn1_real += n_np! * source_weights[1,2,i_weight]
            χn1_imag += n_np! * source_weights[2,2,i_weight]
        end

        n_np += 1
        n_np! *= n_np * one_over_t
    end

    if LH
        t_over_n = t / n
        ϕn1_real -= t_over_n * χn1_imag
        ϕn1_imag += t_over_n * χn1_real

        # ordinarily, χn1 should also be modified here
        # (see Gumerov 2013, Eq. 48)
        # however, we're assuming, for now, that P=n,
        # so χ_{n+1}^m=0 and the correction reduces to 0
    end

    return ϕn0_real, ϕn1_real, ϕn1_imag, χn1_real, χn1_imag
end

function multipole_to_local!(target_branch::Branch{TF}, source_branch, expansion_order, lamb_helmholtz, ε) where TF
    error_check = !isnothing(ε)
    weights_tmp_1 = initialize_expansion(expansion_order + error_check, TF)
    weights_tmp_2 = initialize_expansion(expansion_order + error_check, TF)
    Ts = zeros(TF, length_Ts(expansion_order + error_check))
    eimϕs = zeros(TF, 2, expansion_order+1+error_check)

    return multipole_to_local!(target_branch, source_branch, weights_tmp_1, weights_tmp_2, Ts, eimϕs, ζs_mag, ηs_mag, Hs_π2, expansion_order, lamb_helmholtz, ε)
end

"""
Expects ζs_mag, ηs_mag, and Hs_π2 to be computed a priori.
"""
function multipole_to_local!(target_branch, source_branch, weights_tmp_1, weights_tmp_2, Ts, eimϕs, ζs_mag, ηs_mag, Hs_π2, expansion_order, lamb_helmholtz::Val{LH}, ε::Nothing) where LH
    # extract containers
    source_weights = source_branch.multipole_expansion

    # translation vector
    Δx = target_branch.target_center - source_branch.source_center
    r, θ, ϕ = cartesian_to_spherical(Δx)

    #--- rotate coordinate system ---#

    # rotate about z axis
    rotate_z!(weights_tmp_1, source_weights, eimϕs, ϕ, expansion_order, lamb_helmholtz)

    # rotate about y axis
    # NOTE: the method used here results in an additional rotatation of π about the new z axis
    rotate_multipole_y!(weights_tmp_2, weights_tmp_1, Ts, Hs_π2, ζs_mag, θ, expansion_order, lamb_helmholtz)

    #--- translate along new z axis ---#

    translate_multipole_to_local_z!(weights_tmp_1, weights_tmp_2, r, expansion_order, lamb_helmholtz)

    #--- transform Lamb-Helmholtz decomposition for the new center ---#

    LH && transform_lamb_helmholtz_local!(weights_tmp_1, r, expansion_order)

    #--- back rotate coordinate system ---#

    # back rotate about y axis
    back_rotate_local_y!(weights_tmp_2, weights_tmp_1, Ts, Hs_π2, ηs_mag, expansion_order, lamb_helmholtz)

    # back rotate about z axis and accumulate on target branch
    back_rotate_z!(target_branch.local_expansion, weights_tmp_2, eimϕs, expansion_order, lamb_helmholtz)

end

function dynamic_expansion_order!(weights_tmp_1, weights_tmp_2, Ts, eimϕs, ζs_mag, ηs_mag, source_weights, Hs_π2, expansion_order, lamb_helmholtz::Val{LH}, r, θ, ϕ, r_mp, r_l, ε_tol) where LH

    #--- initialize recursive values ---#

    rinv = one(r) / r
    r_mp_inv = one(r_mp) / r_mp
    r_mp_inv_np2 = r_mp_inv * r_mp_inv * r_mp_inv
    r_l_inv = one(r_l) / r_l
    r_l_nm1 = one(r_l)

    #--- n = 0 coefficients ---#

    # rotate about z axis
    eiϕ_real, eiϕ_imag, eimϕ_real, eimϕ_imag = rotate_z_0!(weights_tmp_1, source_weights, eimϕs, ϕ, expansion_order, lamb_helmholtz)

    # get y-axis Wigner rotation matrix
    sβ, cβ, _1_n = update_Ts_0!(Ts, Hs_π2, θ, expansion_order)

    # perform rotation
    i_ζ, i_T = _rotate_multipole_y_n!(weights_tmp_2, weights_tmp_1, Ts, ζs_mag, expansion_order, lamb_helmholtz, 0, 0, 0)

    # (can't predict error yet, as we need degree n+1 coefficients)
    # preallocate recursive values
    n!_t_np1 = rinv * rinv
    nm1!_inv = 1.0
    np1! = 2.0

    #--- n > 0 and error predictions ---#

    # preallocate here for debugging only
    ε_mp, ε_l = zero(r), zero(r)
    ϕn1_real, ϕn1_imag, ϕn0_real = zero(r), zero(r), zero(r)

    # perform computation
    n = 1
    for _ in 1:expansion_order+1

        # rotate about z axis
        eimϕ_real, eimϕ_imag = rotate_z_n!(weights_tmp_1, source_weights, eimϕs, eiϕ_real, eiϕ_imag, eimϕ_real, eimϕ_imag, lamb_helmholtz, n)

        # get y-axis Wigner rotation matrix
        _1_n = update_Ts_n!(Ts, Hs_π2, sβ, cβ, _1_n, n)

        # perform rotation about the y axis
        # NOTE: the method used here results in an additional rotatation of π about the new z axis
        i_ζ, i_T = _rotate_multipole_y_n!(weights_tmp_2, weights_tmp_1, Ts, ζs_mag, expansion_order, lamb_helmholtz, i_ζ, i_T, n)

        # check multipole error
        in0 = (n*(n+1))>>1 + 1
        ϕn0 = weights_tmp_2[1,1,in0]
        ϕn1_real = weights_tmp_2[1,1,in0+1]
        ϕn1_imag = weights_tmp_2[2,1,in0+1]
        if LH
            χn1_real = weights_tmp_2[1,2,in0+1]
            χn1_imag = weights_tmp_2[2,2,in0+1]
        end

        # calculate multipole error
        ε_mp = (sqrt(ϕn1_real * ϕn1_real + ϕn1_imag * ϕn1_imag) + abs(ϕn0)) * np1! * r_mp_inv_np2
        if LH
            ε_mp += sqrt(χn1_real * χn1_real + χn1_imag * χn1_imag) * np1! * r_mp_inv_np2 * r_mp
        end

        # check local error if multipole error passes
        if ε_mp <= ε_tol

            # extract degree n, order 0-1 coefficients for error prediction
            # this function also performs the lamb-helmholtz transformation
            ϕn0_real, ϕn1_real, ϕn1_imag, χn1_real, χn1_imag = translate_multipole_to_local_z_m01_n(weights_tmp_2, r, rinv, lamb_helmholtz, n!_t_np1, n)

            # calculate local error
            ε_l = (abs(ϕn1_real) + abs(ϕn1_imag) + abs(ϕn0_real)) * r_l_nm1 * nm1!_inv
            if LH
                ε_l += sqrt(χn1_real * χn1_real + χn1_imag * χn1_imag) * r_l_nm1 * r_l * nm1!_inv
            end

            # check total error
            if ε_mp + ε_l <= ε_tol
                return n-1
            end
        end

        # recurse
        if n < expansion_order + 1 # if statement ensures that n == 1 + desired P
                                   # when tolerance is reached (breaks the loop early)
                                   # or that n == 1 + Pmax when tolerance is not reached
            n!_t_np1 *= (n+1) * rinv
            r_l_nm1 *= r_l
            r_mp_inv_np2 *= r_mp_inv
            nm1!_inv /= n

            # increment n
            n += 1
            np1! *= n+1
        end
    end

    @warn "Error tolerance $(ε_tol * ONE_OVER_4π) not reached! Using max expansion order P=$(n-1).\n\tε_mp = $ε_mp, \n\tε_l = $ε_l"

    return n-1
end

"""
Expects ζs_mag, ηs_mag, and Hs_π2 to be computed a priori.
"""
function multipole_to_local!(target_branch, source_branch, weights_tmp_1, weights_tmp_2, Ts, eimϕs, ζs_mag, ηs_mag, Hs_π2, expansion_order, lamb_helmholtz::Val{LH}, ε_tol) where LH
    # extract containers
    source_weights = source_branch.multipole_expansion

    # temporary error variables
    ε_tol *= 4π # multiply from the RHS of the inequality to reduce computational cost

    # translation vector
    Δx = target_branch.target_center - source_branch.source_center
    r, θ, ϕ = cartesian_to_spherical(Δx)

    # multipole error location
    Δx, Δy, Δz = minimum_distance(source_branch.source_center, target_branch.target_center, target_branch.target_box)
    r_mp = sqrt(Δx * Δx + Δy * Δy + Δz * Δz)

    # local error location
    r_l = sum(target_branch.target_box) * 0.33333333333333333333 * sqrt(3)

    #------- rotate multipole coefficients (and determine P) -------#

    expansion_order = dynamic_expansion_order!(weights_tmp_1, weights_tmp_2, Ts, eimϕs, ζs_mag, ηs_mag, source_weights, Hs_π2, expansion_order, lamb_helmholtz, r, θ, ϕ, r_mp, r_l, ε_tol)

    #------- translate multipole to local expansion -------#

    #--- translate along new z axis ---#

    translate_multipole_to_local_z!(weights_tmp_1, weights_tmp_2, r, expansion_order, lamb_helmholtz)

    #--- transform Lamb-Helmholtz decomposition for the new center ---#

    LH && transform_lamb_helmholtz_local!(weights_tmp_1, r, expansion_order)

    #--- back rotate coordinate system ---#

    # back rotate about y axis
    back_rotate_local_y!(weights_tmp_2, weights_tmp_1, Ts, Hs_π2, ηs_mag, expansion_order, lamb_helmholtz)

    # back rotate about z axis and accumulate on target branch
    back_rotate_z!(target_branch.local_expansion, weights_tmp_2, eimϕs, expansion_order, lamb_helmholtz)

end

"defaults to no error prediction"
multipole_to_local!(target_branch, source_branch, weights_tmp_1, weights_tmp_2, Ts, eimϕs, ζs_mag, ηs_mag, Hs_π2, expansion_order, lamb_helmholtz) =
    multipole_to_local!(target_branch, source_branch, weights_tmp_1, weights_tmp_2, Ts, eimϕs, ζs_mag, ηs_mag, Hs_π2, expansion_order, lamb_helmholtz, nothing)

#------- LOCAL TO LOCAL -------#

"""
Overwrites translated_weights
"""
function translate_local_z!(translated_weights, source_weights, t, expansion_order, ::Val{LH}) where LH
#function translate_local_z!(translated_weights, source_weights, t, ::Val{P}, ::Val{LH}) where {P,LH}
    _t = -t
    i = 1
    for n in 0:expansion_order
    	for m in 0:n
    		# inner summation
            val1_real = zero(eltype(translated_weights))
            val1_imag = zero(eltype(translated_weights))

            if LH
                val2_real = zero(eltype(translated_weights))
                val2_imag = zero(eltype(translated_weights))
            end
    		n_np = 1
    		_t_n_n! = one(t)
    		for np in n:expansion_order
                i_weight = harmonic_index(np,m)
                val1_real += _t_n_n! * source_weights[1,1,i_weight]
                val1_imag += _t_n_n! * source_weights[2,1,i_weight]
                if LH
                    val2_real += _t_n_n! * source_weights[1,2,i_weight]
                    val2_imag += _t_n_n! * source_weights[2,2,i_weight]
                end
    			_t_n_n! *= _t / n_np
    			n_np += 1
    		end

            # set translated coefficient
            translated_weights[1,1,i] = val1_real
            translated_weights[2,1,i] = val1_imag
            if LH
                translated_weights[1,2,i] = val2_real
                translated_weights[2,2,i] = val2_imag
            end

            # increment index
            i += 1
    	end
    end
end

"""
Expects ηs_mag and Hs_π2 to be precomputed. Ts and eimϕs are computed here.
"""
function local_to_local!(target_branch, source_branch, weights_tmp_1, weights_tmp_2, Ts, eimϕs, ηs_mag, Hs_π2, expansion_order, lamb_helmholtz::Val{LH}) where LH
    # extract containers
    source_weights = source_branch.local_expansion

    # translation vector
    Δx = target_branch.target_center - source_branch.target_center
    r, θ, ϕ = cartesian_to_spherical(Δx)

    #--- rotate coordinate system ---#

    # rotate about z axis (and compute eimϕs)
    rotate_z!(weights_tmp_1, source_weights, eimϕs, ϕ, expansion_order, lamb_helmholtz)

    # rotate about y axis (and compute Ts)
    # NOTE: the method used here results in an additional rotatation of π about the new z axis
    rotate_local_y!(weights_tmp_2, weights_tmp_1, Ts, Hs_π2, ηs_mag, θ, expansion_order, lamb_helmholtz)

    #--- translate along new z axis ---#

    translate_local_z!(weights_tmp_1, weights_tmp_2, r, expansion_order, lamb_helmholtz)

    #--- Lamb-Helmholtz transformation ---#

    LH && transform_lamb_helmholtz_local!(weights_tmp_1, r, expansion_order)

    #--- back rotate coordinate system ---#

    # back rotate about y axis
    back_rotate_local_y!(weights_tmp_2, weights_tmp_1, Ts, Hs_π2, ηs_mag, expansion_order, lamb_helmholtz)

    # back rotate about z axis and accumulate result to target branch
    back_rotate_z!(target_branch.local_expansion, weights_tmp_2, eimϕs, expansion_order, lamb_helmholtz)

end

