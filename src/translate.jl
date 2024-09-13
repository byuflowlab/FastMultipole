#------- MULTIPOLE TO MULTIPOLE -------#

function update_multipole!(target_branch, multipole_weights, expansion_order::Val{P}) where P
    i = 1
    for n in 0:P
        for m in 0:n
            target_branch.multipole_expansion[1,1,i] += multipole_weights[1,1,i]
            target_branch.multipole_expansion[2,1,i] += multipole_weights[2,1,i]
            i += 1
        end
    end
end

"""
Overwrites translated_weights
"""
function translate_multipole_z!(translated_weights, source_weights, t, expansion_order::Val{P}) where P
    _t = -t
    i = 1
    for n in 0:P
    	for m in 0:n
    		# inner summation
            val_real = zero(eltype(translated_weights))
            val_imag = zero(eltype(translated_weights))
    		n_np! = 1.0
    		n_np = 1
    		_t_n_np = one(t)
    		for np in n:-1:m
                tmp = _t_n_np / n_np!
                val_real += tmp * source_weights[1,1,harmonic_index(np,m)]
                val_imag += tmp * source_weights[2,1,harmonic_index(np,m)]
    			_t_n_np *= _t
    			n_np! *= n_np
    			n_np += 1
    		end

            # set translated coefficient
            translated_weights[1,1,i] = val_real
            translated_weights[2,1,i] = val_imag

            # increment index
            i += 1
    	end
    end
end

function multipole_to_multipole!(target_branch, source_branch, weights_tmp_1, weights_tmp_2, Ts, eimϕs, ζs_mag, Hs_π2, expansion_order)
    # extract containers
    source_weights = source_branch.multipole_expansion

    # translation vector
    Δx = target_branch.center - source_branch.center
    r, θ, ϕ = cartesian_to_spherical(Δx)

    #--- rotate coordinate system ---#

    # rotate about z axis
    rotate_z!(weights_tmp_1, source_weights, eimϕs, ϕ, expansion_order)

    # rotate about y axis
    # NOTE: the method used here results in an additional rotatation of π about the new z axis
    rotate_multipole_y!(weights_tmp_2, weights_tmp_1, Ts, Hs_π2, ζs_mag, θ, expansion_order)

    #--- translate along new z axis ---#

    translate_multipole_z!(weights_tmp_1, weights_tmp_2, r, expansion_order)

    #--- back rotate coordinate system ---#

    # back rotate about y axis
    back_rotate_multipole_y!(weights_tmp_2, weights_tmp_1, Ts, Hs_π2, ζs_mag, expansion_order)

    # back rotate about z axis
    back_rotate_z!(weights_tmp_1, weights_tmp_2, eimϕs, expansion_order)

    #--- store coefficients ---#

    update_multipole!(target_branch, weights_tmp_1, expansion_order)

end

#------- MULTIPOLE TO LOCAL -------#

function update_local!(target_branch, local_weights, expansion_order::Val{P}) where P
    i = 1
    for n in 0:P
        for m in 0:n
            target_branch.local_expansion[1,1,i] += local_weights[1,1,i]
            target_branch.local_expansion[2,1,i] += local_weights[2,1,i]
            i += 1
        end
    end
end

"""
Overwrites translated_weights
"""
function translate_multipole_to_local_z!(translated_weights, source_weights, t, expansion_order::Val{P}) where P
    one_over_t = 1/t
    n!_t_np1 = one_over_t
    i = 1
    for n in 0:P
        n_m!_t_nmp1 = n!_t_np1
        n_m = n
    	for m in 0:n
            # inner summation
            val_real = zero(eltype(translated_weights))
            val_imag = zero(eltype(translated_weights))
    		n_np! = n_m!_t_nmp1
            n_np = n_m
    		for np in m:P
                val_real += n_np! * source_weights[1,1,harmonic_index(np,m)]
                val_imag += n_np! * source_weights[2,1,harmonic_index(np,m)]

                n_np += 1
                n_np! *= n_np * one_over_t
    		end

            # set translated coefficient
            translated_weights[1,1,i] = val_real
            translated_weights[2,1,i] = val_imag

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
Expects ζs_mag, ηs_mag, and Hs_π2 to be computed a priori.
"""
function multipole_to_local!(target_branch, source_branch, weights_tmp_1, weights_tmp_2, Ts, eimϕs, ζs_mag, ηs_mag, Hs_π2, expansion_order)
    # extract containers
    source_weights = source_branch.multipole_expansion

    # translation vector
    Δx = target_branch.center - source_branch.center
    r, θ, ϕ = cartesian_to_spherical(Δx)

    #--- rotate coordinate system ---#

    # rotate about z axis
    rotate_z!(weights_tmp_1, source_weights, eimϕs, ϕ, expansion_order)

    # rotate about y axis
    # NOTE: the method used here results in an additional rotatation of π about the new z axis
    rotate_multipole_y!(weights_tmp_2, weights_tmp_1, Ts, Hs_π2, ζs_mag, θ, expansion_order)

    #--- translate along new z axis ---#

    translate_multipole_to_local_z!(weights_tmp_1, weights_tmp_2, r, expansion_order)

    #--- back rotate coordinate system ---#

    # back rotate about y axis
    back_rotate_local_y!(weights_tmp_2, weights_tmp_1, Ts, Hs_π2, ηs_mag, expansion_order)

    # back rotate about z axis
    back_rotate_z!(weights_tmp_1, weights_tmp_2, eimϕs, expansion_order)

    #--- store coefficients ---#

    update_local!(target_branch, weights_tmp_1, expansion_order)

end

#------- LOCAL TO LOCAL -------#

"""
Overwrites translated_weights
"""
function translate_local_z!(translated_weights, source_weights, t, expansion_order)
    _t = -t
    i = 1
    for n in 0:expansion_order
    	for m in 0:n
    		# inner summation
            val_real = zero(eltype(translated_weights))
            val_imag = zero(eltype(translated_weights))
    		n_np = 1
    		_t_n_n! = one(t)
    		for np in n:expansion_order
                val_real += _t_n_n! * source_weights[1,1,harmonic_index(np,m)]
                val_imag += _t_n_n! * source_weights[2,1,harmonic_index(np,m)]
    			_t_n_n! *= _t / n_np
    			n_np += 1
    		end

            # set translated coefficient
            translated_weights[1,1,i] = val_real
            translated_weights[2,1,i] = val_imag

            # increment index
            i += 1
    	end
    end
end

"""
Expects ηs_mag and Hs_π2 to be precomputed. Ts and eimϕs are computed here.
"""
function local_to_local!(target_branch, source_branch, weights_tmp_1, weights_tmp_2, Ts, eimϕs, ηs_mag, Hs_π2, expansion_order)
    # extract containers
    source_weights = source_branch.local_expansion

    # translation vector
    Δx = target_branch.center - source_branch.center
    r, θ, ϕ = cartesian_to_spherical(Δx)

    #--- rotate coordinate system ---#

    # rotate about z axis (and compute eimϕs)
    rotate_z!(weights_tmp_1, source_weights, eimϕs, ϕ, expansion_order)

    # rotate about y axis (and compute Ts)
    # NOTE: the method used here results in an additional rotatation of π about the new z axis
    rotate_local_y!(weights_tmp_2, weights_tmp_1, Ts, Hs_π2, ηs_mag, θ, expansion_order)

    #--- translate along new z axis ---#

    translate_local_z!(weights_tmp_1, weights_tmp_2, r, expansion_order)

    #--- back rotate coordinate system ---#

    # back rotate about y axis
    back_rotate_local_y!(weights_tmp_2, weights_tmp_1, Ts, Hs_π2, ηs_mag, expansion_order)

    # back rotate about z axis
    back_rotate_z!(weights_tmp_1, weights_tmp_2, eimϕs, expansion_order)

    #--- store coefficients ---#

    update_local!(target_branch, weights_tmp_1, expansion_order)

end

