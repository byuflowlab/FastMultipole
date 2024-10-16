#------- MULTIPOLE TO MULTIPOLE -------#

"""
Overwrites translated_weights
"""
function translate_multipole_z!(translated_weights, source_weights, t, P, ::Val{LH}) where LH
#function translate_multipole_z!(translated_weights, source_weights, t, expansion_order::Val{P}, ::Val{LH}) where {P,LH}
    _t = -t
    i = 1
    for n in 0:P
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

function transform_lamb_helmholtz_multipole!(multipole_expansion, r, P)
    i_P_m = harmonic_index(P, 0)
    for m in 0:P
        # declare recursive variable for inner loop over n
        χ̂_nm1_real = multipole_expansion[1,2,i_P_m]
        χ̂_nm1_imag = multipole_expansion[2,2,i_P_m]

        # declare index
        i_n_m = i_P_m

        for n in P:-1:max(m,1)
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

function transform_lamb_helmholtz_local!(local_expansion, r, P)
    i_m_m = 1
    for m in 0:P
        # declare recursive variable for inner loop over n
        χ̂_np1_real = local_expansion[1,2,i_m_m]
        χ̂_np1_imag = local_expansion[2,2,i_m_m]

        # declare index for fast access
        i_n_m = i_m_m

        for n in m:P # skip n=0
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

            if n < P # χ_{n+1}^m = 0 when n+1>P
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

function multipole_to_multipole!(target_branch, source_branch, weights_tmp_1, weights_tmp_2, Ts, eimϕs, ζs_mag, Hs_π2, expansion_order::Val{P}, lamb_helmholtz::Val{LH}) where {P,LH}
    # extract containers
    source_weights = source_branch.multipole_expansion

    # translation vector
    Δx = target_branch.center - source_branch.center
    r, θ, ϕ = cartesian_to_spherical(Δx)

    #--- rotate coordinate system ---#

    # rotate about z axis
    rotate_z!(weights_tmp_1, source_weights, eimϕs, ϕ, P, lamb_helmholtz)

    # rotate about y axis
    # NOTE: the method used here results in an additional rotatation of π about the new z axis
    rotate_multipole_y!(weights_tmp_2, weights_tmp_1, Ts, Hs_π2, ζs_mag, θ, P, lamb_helmholtz)

    #--- translate along new z axis ---#

    translate_multipole_z!(weights_tmp_1, weights_tmp_2, r, P, lamb_helmholtz)

    #--- transform Lamb-Helmholtz decomposition for the new center ---#

    LH && transform_lamb_helmholtz_multipole!(weights_tmp_1, r, P)

    #--- back rotate coordinate system ---#

    # back rotate about y axis
    back_rotate_multipole_y!(weights_tmp_2, weights_tmp_1, Ts, ζs_mag, P, lamb_helmholtz)

    # back rotate about z axis and accumulate on target branch
    back_rotate_z!(target_branch.multipole_expansion, weights_tmp_2, eimϕs, P, lamb_helmholtz)

end

#------- MULTIPOLE TO LOCAL -------#

"""
Overwrites translated_weights
"""
function translate_multipole_to_local_z!(translated_weights, source_weights, t, P, ::Val{LH}) where LH
#function translate_multipole_to_local_z!(translated_weights, source_weights, t, ::Val{P}, ::Val{LH}) where {P,LH}
    one_over_t = 1/t
    n!_t_np1 = one_over_t
    i = 1
    for n in 0:P
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
    		for np in m:P
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

@inline function get_dynamic_expansion_order(center_distance, source_radius, target_radius, relative_error, expansion_order::Val{Pmax}) where Pmax
    target = relative_error * (center_distance - source_radius - target_radius)
    f1 = target_radius / (center_distance-source_radius)
    f2 = source_radius / (center_distance-target_radius)
    t1 = target_radius
    t2 = source_radius
    for p in 0:Pmax-1
        (t1 + t2 < target) && (return p)
        t1 *= f1
        t2 *= f2
    end
    return Pmax
end

@inline function get_dynamic_expansion_order(center_distance, source_radius, target_radius, relative_error::Nothing, expansion_order::Val{P}) where P
    return P
end

"""
Expects ζs_mag, ηs_mag, and Hs_π2 to be computed a priori.
"""
function multipole_to_local!(target_branch, source_branch, weights_tmp_1, weights_tmp_2, Ts, eimϕs, ζs_mag, ηs_mag, Hs_π2, expansion_order::Val, lamb_helmholtz::Val{LH}, relative_error) where LH
    # extract containers
    source_weights = source_branch.multipole_expansion

    # translation vector
    Δx = target_branch.center - source_branch.center
    r, θ, ϕ = cartesian_to_spherical(Δx)

    # choose expansion order
    dynamic_expansion_order = get_dynamic_expansion_order(r, source_branch.radius, target_branch.radius, relative_error, expansion_order)

    #--- rotate coordinate system ---#

    # rotate about z axis
    rotate_z!(weights_tmp_1, source_weights, eimϕs, ϕ, dynamic_expansion_order, lamb_helmholtz)

    # rotate about y axis
    # NOTE: the method used here results in an additional rotatation of π about the new z axis
    rotate_multipole_y!(weights_tmp_2, weights_tmp_1, Ts, Hs_π2, ζs_mag, θ, dynamic_expansion_order, lamb_helmholtz)

    #--- translate along new z axis ---#

    translate_multipole_to_local_z!(weights_tmp_1, weights_tmp_2, r, dynamic_expansion_order, lamb_helmholtz)

    #--- transform Lamb-Helmholtz decomposition for the new center ---#

    LH && transform_lamb_helmholtz_local!(weights_tmp_1, r, dynamic_expansion_order)

    #--- back rotate coordinate system ---#

    # back rotate about y axis
    back_rotate_local_y!(weights_tmp_2, weights_tmp_1, Ts, Hs_π2, ηs_mag, dynamic_expansion_order, lamb_helmholtz)

    # back rotate about z axis and accumulate on target branch
    back_rotate_z!(target_branch.local_expansion, weights_tmp_2, eimϕs, dynamic_expansion_order, lamb_helmholtz)

end

#------- LOCAL TO LOCAL -------#

"""
Overwrites translated_weights
"""
function translate_local_z!(translated_weights, source_weights, t, P, ::Val{LH}) where LH
#function translate_local_z!(translated_weights, source_weights, t, ::Val{P}, ::Val{LH}) where {P,LH}
    _t = -t
    i = 1
    for n in 0:P
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
    		for np in n:P
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
function local_to_local!(target_branch, source_branch, weights_tmp_1, weights_tmp_2, Ts, eimϕs, ηs_mag, Hs_π2, expansion_order::Val{P}, lamb_helmholtz::Val{LH}) where {P,LH}
    # extract containers
    source_weights = source_branch.local_expansion

    # translation vector
    Δx = target_branch.center - source_branch.center
    r, θ, ϕ = cartesian_to_spherical(Δx)

    #--- rotate coordinate system ---#

    # rotate about z axis (and compute eimϕs)
    rotate_z!(weights_tmp_1, source_weights, eimϕs, ϕ, P, lamb_helmholtz)

    # rotate about y axis (and compute Ts)
    # NOTE: the method used here results in an additional rotatation of π about the new z axis
    rotate_local_y!(weights_tmp_2, weights_tmp_1, Ts, Hs_π2, ηs_mag, θ, P, lamb_helmholtz)

    #--- translate along new z axis ---#

    translate_local_z!(weights_tmp_1, weights_tmp_2, r, P, lamb_helmholtz)

    #--- Lamb-Helmholtz transformation ---#

    LH && transform_lamb_helmholtz_local!(weights_tmp_1, r, P)

    #--- back rotate coordinate system ---#

    # back rotate about y axis
    back_rotate_local_y!(weights_tmp_2, weights_tmp_1, Ts, Hs_π2, ηs_mag, P, lamb_helmholtz)

    # back rotate about z axis and accumulate result to target branch
    back_rotate_z!(target_branch.local_expansion, weights_tmp_2, eimϕs, P, lamb_helmholtz)

end

