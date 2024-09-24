#------- GENERAL FUNCTIONS -------#

#--- z axis rotation matrices (diagonal) ---#

function update_eimϕs!(eimϕs, ϕ, ::Val{expansion_order}) where expansion_order
    eiϕ_imag, eiϕ_real = sincos(ϕ)
    eimϕ_real, eimϕ_imag = one(ϕ), zero(ϕ)
    for m in 0:expansion_order
        eimϕs[1,m+1] = eimϕ_real
        eimϕs[2,m+1] = eimϕ_imag
        eimϕ_real_tmp = eimϕ_real
        eimϕ_real = eimϕ_real * eiϕ_real - eimϕ_imag * eiϕ_imag
        eimϕ_imag = eimϕ_real_tmp * eiϕ_imag + eimϕ_imag * eiϕ_real
    end
end

#--- precompute y axis rotation matrices by π/2 ---#

function length_Hs(expansion_order)
    p = expansion_order
    n = p+1
    return div(n*(n+1)*(n+2), 6)
end

function length_H(n)
    return ((n+1)*(n+2))>>1
end

function get_H(Hs, n)
    if n==0
        return view(Hs,1:1)
    else
        i1 = length_Hs(n-1)+1
        i2 = length_Hs(n)
        return view(Hs,i1:i2)
    end
end

@inline function H_index(mp, m)
    return (m * (m+1))>>1 + mp + 1
end

@inline function Pnm_index(n, m)
    return (n * (n+1))>>1 + m + 1
end

function compute_Pnms(expansion_order, s, c)
    p = expansion_order
    Pnms = zeros(((p+1)*(p+2))>>1)
    Pn0, fact = 1.0, 1.0
    for m in 0:p
        Pnm = Pn0
        i = Pnm_index(m, m)
        Pnms[i] = Pnm
        tmp1 = Pnm
        Pnm = c * (2 * m + 1) * tmp1
        for n in m+1:p
            i = Pnm_index(n, m)
            Pnms[i] = Pnm
            tmp2 = tmp1
            tmp1 = Pnm
            Pnm = (c * (2 * n + 1) * tmp1 - (n + m) * tmp2) / (n - m + 1)
        end
        Pn0 = -Pn0 * fact * s
        fact += 2
    end

    return Pnms
end

function compute_Hs!(Hs, β, n_start, n_end)
    @assert n_start != 0

    sβ, cβ = sincos(β)
    Pnms = compute_Pnms(n_end+1, sβ, cβ)
    H_np1 = zeros(n_end+2)

    # b
    bnms = zeros(2(n_end+1)+1)

    # d
    dnms = zeros(n_end+1)

    # Hs
    i_start = length_Hs(n_start-1) + 1
    i_end = i_start + length_H(n_start) - 1
    for n in n_start:n_end
        #H = get_H(Hs, n)
        H = view(Hs, i_start:i_end)

        # H_n_0m
        mp = 0
        _1_m, fact = 1, 1
        for m in 0:n
            # next value
            H[H_index(mp,m)] = _1_m * fact * Pnms[Pnm_index(n, m)]

            # recurse
            _1_m *= -1
            fact /= sqrt((n+m+1)*(n-m))
        end

        # H_np1_0m
        _1_m, fact = 1, 1
        for m in 0:n+1
            # next value
            H_np1[m+1] = _1_m * fact * Pnms[Pnm_index(n+1,m)]

            # recurse
            _1_m *= -1
            fact /= sqrt((n+1-m)*(n+m+2))
        end

        # b_np1_m
        b_n = n+1
        for m in -b_n:-1
            bnms[b_n + m + 1] = -sqrt((b_n-m-1)/(2b_n-1)*(b_n-m)/(2b_n+1))
        end
        for m in 0:b_n
            bnms[b_n + m + 1] = sqrt((b_n-m-1)/(2b_n-1)*(b_n-m)/(2b_n+1))
        end

        # H_n_1m
        one_m_cβ_2 = (1.0 - cβ)/2
        one_p_cβ_2 = (1.0 + cβ)/2
        for m in 1:n
            b_np1__mm1 = bnms[b_n-m]
            b_np1_mm1 = bnms[b_n+m]
            b_np1_0 = bnms[b_n+1]
            a_n_m = sqrt((n+1+m)/(2n+1)*(n+1-m)/(2n+3))
            H[H_index(1,m)] = (b_np1__mm1 * one_m_cβ_2 * H_np1[m+2] - b_np1_mm1 * one_p_cβ_2 * H_np1[m] - a_n_m * sβ * H_np1[m+1]) / b_np1_0
        end

        # update dnms
        for m in 0:n
            dnms[m+1] = sqrt((n-m)*(n+m+1)) / 2
        end

        # H_n_mp+1_m
        for mp in 1:n-1
            for m in mp+1:n-1
                val = dnms[mp] * H[H_index(mp-1,m)] - dnms[m] * H[H_index(mp,m-1)] + dnms[m+1] * H[H_index(mp,m+1)]
                H[H_index(mp+1,m)] = val / dnms[mp+1]
            end
            m = n
            val = dnms[mp] * H[H_index(mp-1,m)] - dnms[m] * H[H_index(mp,m-1)]
            H[H_index(mp+1,m)] = val / dnms[mp+1]
        end

        # recurse indices
        i_start = i_end + 1
        i_end = i_start + length_H(n+1) - 1

    end
end

function update_Hs_π2!(Hs_π2, expansion_order::Val{P}) where P
    l = length(Hs_π2)
    l_desired = length_Hs(P)

    if l_desired > l # we need to compute more matrices
        # determine how many we already have
        n_already, l_already = 0, 0
        for n in 0:P
            l_already += (n+1)^2
            l_already == l && (n_already = n)
        end

        # calculate the rest
        β = pi/2
        resize!(Hs_π2, l_desired)
        compute_Hs!(Hs_π2, β, n_already+1, P)

    end # otherwise, we already have enough, so do nothing
end

#--- use precomputed H matrices to compute rotations for arbitrary angle ---#

function length_Ts(expansion_order)
    p = expansion_order
    return div((p+1)*(p+2)*(2p+3),6)
end

function length_T(n)
    np1 = n+1
    return np1*np1
end

function get_T(Ts, n)
    if n == 0
        return view(Ts, 1:1)
    else
        i1 = length_Ts(n-1)+1
        i2 = length_Ts(n)
        return view(Ts,i1:i2)
    end
end

@inline function T_index(mp, m)
    return m*m + mp + m + 1
end

@inline function get_scalar(m_mp)
    isodd(m_mp) && (m_mp += 1)
    if iseven(m_mp>>1)
        return 1.0
    else
        return -1.0
    end
end

function update_Ts!(Ts, Hs_π2, β::TF, ::Val{expansion_order}) where {TF,expansion_order}
    p = expansion_order

    # initialize recursive indices
    i_H_start = length_Hs(0) + 1 # 2
    i_H_end = i_H_start + 2 # n=1 matrix contains 3 unique entries
    i_T_start = 2
    i_T_end = i_T_start + 3 # n=1 T matrix contains 4 unique entries

    # other recursive values
    sβ, cβ = sincos(β)
    _1_n = -1.0

    # rotate each expansion order
    Ts[1] = one(eltype(Ts))
    @inbounds for n in 1:p
        # Wigner matrix about y axis by π/2
        H_π2 = view(Hs_π2, i_H_start:i_H_end)

        # Wigner matrix by β
        T = view(Ts, i_T_start:i_T_end)

        for m in 0:n
            H_π2_n_m_0 = H_π2[H_index(0,m)]
            _1_mp_odd = 1.0
            m_mp = m
            _1_n_mp = _1_n

            for mp in 0:m
                H_π2_n_mp_0 = H_π2[H_index(0,mp)]
                val_positive_mp = zero(eltype(Ts))
                val_negative_mp = zero(eltype(Ts))
                s_νm1_β = 0.0
                c_νm1_β = 1.0
                m_mp_even = iseven(m_mp)
                scalar = get_scalar(m_mp)
                _1_n_mp_ν = -_1_n_mp

                for ν in 1:n
                    # sin and cos (νβ)
                    c_νβ = cβ * c_νm1_β - sβ * s_νm1_β
                    s_νβ = sβ * c_νm1_β + cβ * s_νm1_β

                    # alternating sine/cosine
                    sc_νβ_m_mp = m_mp_even ? scalar * c_νβ : scalar * s_νβ

                    # accumulate
                    i, j = minmax(mp,ν)
                    H_π2_n_mp_ν = H_π2[H_index(i,j)]
                    i, j = minmax(m,ν)
                    H_π2_n_m_ν = H_π2[H_index(i,j)]

                    # leverage symmetry to compute p/m mp
                    val = H_π2_n_mp_ν * H_π2_n_m_ν * sc_νβ_m_mp #
                    val_positive_mp += val
                    val_negative_mp += val * _1_n_mp_ν * _1_mp_odd #

                    # recurse
                    c_νm1_β = c_νβ
                    s_νm1_β = s_νβ
                    _1_n_mp_ν = -_1_n_mp_ν
                end

                val_positive_mp *= 2.0
                val_negative_mp *= 2.0
                scalar = m_mp_even ? scalar : 0.0
                val = H_π2_n_m_0 * H_π2_n_mp_0 * scalar
                val_positive_mp += val
                val_negative_mp += val * _1_n_mp * _1_mp_odd #
                T[T_index(mp,m)] = val_positive_mp
                T[T_index(-mp,m)] = val_negative_mp

                # recurse
                _1_n_mp = -_1_n_mp
                _1_mp_odd = -_1_mp_odd
                m_mp += 1
            end
        end

        # recurse
        i_H_start = i_H_end + 1
        i_H_end = i_H_start + length_H(n+1) - 1
        i_T_start = i_T_end + 1
        i_T_end = i_T_start + length_T(n+1) - 1
        _1_n = -_1_n
    end
end

#--- functions to actually rotate ---#

"""
Performs a z-axis rotation of the supplied solid harmonic coefficients. Computes e^{imϕ} as well.
"""
function rotate_z!(rotated_weights, source_weights, eimϕs, ϕ, expansion_order::Val{P}, ::ExpansionSwitch{SP,VP}) where {P,SP,VP}

    update_eimϕs!(eimϕs, ϕ, expansion_order)

    i_weight = 1
    for n in 0:P
        for m in 0:n
            # get e^{imΔϕ}
            eimϕ_real, eimϕ_imag = eimϕs[1,m+1], eimϕs[2,m+1]

            # rotate coefficients
            if SP
                Xnm_real, Xnm_imag = source_weights[1,1,i_weight], source_weights[2,1,i_weight]
                rotated_weights[1,1,i_weight] = Xnm_real * eimϕ_real - Xnm_imag * eimϕ_imag
                rotated_weights[2,1,i_weight] = Xnm_real * eimϕ_imag + Xnm_imag * eimϕ_real
            end
            if VP
                Xnm_real, Xnm_imag = source_weights[1,2,i_weight], source_weights[2,2,i_weight]
                rotated_weights[1,2,i_weight] = Xnm_real * eimϕ_real - Xnm_imag * eimϕ_imag
                rotated_weights[2,2,i_weight] = Xnm_real * eimϕ_imag + Xnm_imag * eimϕ_real
                Xnm_real, Xnm_imag = source_weights[1,3,i_weight], source_weights[2,3,i_weight]
                rotated_weights[1,3,i_weight] = Xnm_real * eimϕ_real - Xnm_imag * eimϕ_imag
                rotated_weights[2,3,i_weight] = Xnm_real * eimϕ_imag + Xnm_imag * eimϕ_real
            end

            # increment index
            i_weight += 1
        end
    end
end

"""
Assumes eimϕs have already been computed. DOES NOT overwrite rotated weights (unlike other rotate functions); rather, accumulates on top of it.
"""
function back_rotate_z!(rotated_weights, source_weights, eimϕs, expansion_order::Val{P}, ::ExpansionSwitch{SP,VP}) where {P,SP,VP}
    i_weight = 1
    for n in 0:P
        for m in 0:n
            # get e^{-imΔϕ}
            eimϕ_real, eimϕ_imag = eimϕs[1,m+1], -eimϕs[2,m+1]

            # rotate coefficients
            if SP
                Xnm_real, Xnm_imag = source_weights[1,1,i_weight], source_weights[2,1,i_weight]
                rotated_weights[1,1,i_weight] += Xnm_real * eimϕ_real - Xnm_imag * eimϕ_imag
                rotated_weights[2,1,i_weight] += Xnm_real * eimϕ_imag + Xnm_imag * eimϕ_real
            end
            if VP
                Xnm_real, Xnm_imag = source_weights[1,2,i_weight], source_weights[2,2,i_weight]
                rotated_weights[1,2,i_weight] += Xnm_real * eimϕ_real - Xnm_imag * eimϕ_imag
                rotated_weights[2,2,i_weight] += Xnm_real * eimϕ_imag + Xnm_imag * eimϕ_real
                Xnm_real, Xnm_imag = source_weights[1,3,i_weight], source_weights[2,3,i_weight]
                rotated_weights[1,3,i_weight] += Xnm_real * eimϕ_real - Xnm_imag * eimϕ_imag
                rotated_weights[2,3,i_weight] += Xnm_real * eimϕ_imag + Xnm_imag * eimϕ_real
            end

            # increment index
            i_weight += 1
        end
    end
end

#------- MULTIPOLE ROTATIONS -------#

function length_ζs(expansion_order)
    p = expansion_order
    return div((p+1)*(p+2)*(2p+3),6)
end

function update_ζs_mag!(ζs_mag, expansion_order::Val{P}) where P
    l = length(ζs_mag)
    l_desired = length_ζs(P)

    if l_desired > l # we need to compute more matrices
        # determine how many we already have
        n_already, l_already = 0, 0
        for n in 0:P
            l_already += (n+1)^2
            l_already == l && (n_already = n)
        end

        # calculate the rest
        resize!(ζs_mag, l_desired)
        update_ζs_mag!(ζs_mag, n_already+1, P)

    end # otherwise, we already have enough, so do nothing
end

@inline function ζ_sign(magnitude, mp, m)
    mod = (abs(mp) - abs(m)) % 4
    mod < 0 && (mod += 4)
    if mod == 0
        return magnitude, 0.0
    elseif mod == 1
        return 0.0, magnitude
    elseif mod == 2
        return -magnitude, 0.0
    elseif mod == 3
        return 0.0, -magnitude
    end
end

function update_ζs_mag!(ζs, n_start, n_end)
    if n_start == 0
        ζs[1] = 1.0
        n_start = 1
    end
    for n in n_start:n_end
        ζ = get_ζ(ζs, n)
        i = 1
        mag_m = 1.0 # magnitude of ζ
        for m in 0:n
            mag_mp = mag_m

            for mp in 0:n
                # store magnitude (we'll compute the sign on the fly)
                ζ[i] = mag_mp

                # recurse over mp
                num = n + mp + 1
                denom = mp+1 < n ? n - mp : 1
                mag_mp *= sqrt(num/denom)

                # augment index
                i += 1
            end

            # recurse over m
            num = m < n - 1 ? n - m : 1
            denom = n+m+1
            mag_m *= sqrt(num/denom)
        end
    end
end

function get_ζ(ζs, n)
    i1 = length_ζs(n-1)+1
    i2 = length_ζs(n)
    return view(ζs, i1:i2)
end

#--- functions to actually rotate ---#

function _rotate_multipole_y!(rotated_weights, source_weights, Ts, ζs_mag, ::Val{P}, ::ExpansionSwitch{SP,VP}) where {P,SP,VP}
    # reset container
    rotated_weights .= zero(eltype(rotated_weights))

    # rotate each order
    i_ζ = 0
    i_T = 0
    #println("=== START ===")
    @inbounds for n in 0:P
        #println("== n=$n ==")
        for m in 0:n
            #println("= m=$m =")
            if SP
                val1_real = zero(eltype(Ts))
                val1_imag = zero(eltype(Ts))
            end
            if VP
                val2_real = zero(eltype(Ts))
                val2_imag = zero(eltype(Ts))
                val3_real = zero(eltype(Ts))
                val3_imag = zero(eltype(Ts))
            end
            i_ζ += n+1
            i_ζ_mp = i_ζ

            for mp in -n:-m-1 # -n <= mp < -m
                #println("\tm'=$mp")
                # T_n^{m',m} = T_n^{m,m'} = T_n^{-m,-m'}
                # @inline function T_index(mp, m)
                #     return m*m + mp + m + 1
                # end
                T_mp_m = Ts[i_T + T_index(-m,-mp)] # should always be real
                #@show T_index(-m,-mp)

                # ζ_n^{m',m} = β_n^{m'} / β_n^m
                ζ_real, ζ_imag = ζ_sign(ζs_mag[i_ζ_mp] * T_mp_m, mp, m)

                # current weights
                i_weight = harmonic_index(n, -mp)

                # sum contribution
                if SP
                    # (-1)^mp * complex conjugate for -mp
                    M_n_mp_real, M_n_mp_imag = source_weights[1,1,i_weight], -source_weights[2,1,i_weight]
                    if isodd(mp)
                        M_n_mp_real = -M_n_mp_real
                        M_n_mp_imag = -M_n_mp_imag
                    end
                    val1_real += M_n_mp_real * ζ_real - M_n_mp_imag * ζ_imag
                    val1_imag += M_n_mp_real * ζ_imag + M_n_mp_imag * ζ_real
                end
                if VP
                    # (-1)^mp * complex conjugate for -mp
                    M_n_mp_real, M_n_mp_imag = source_weights[1,2,i_weight], -source_weights[2,2,i_weight]
                    if isodd(mp)
                        M_n_mp_real = -M_n_mp_real
                        M_n_mp_imag = -M_n_mp_imag
                    end
                    val2_real += M_n_mp_real * ζ_real - M_n_mp_imag * ζ_imag
                    val2_imag += M_n_mp_real * ζ_imag + M_n_mp_imag * ζ_real

                    # (-1)^mp * complex conjugate for -mp
                    M_n_mp_real, M_n_mp_imag = source_weights[1,3,i_weight], -source_weights[2,3,i_weight]
                    if isodd(mp)
                        M_n_mp_real = -M_n_mp_real
                        M_n_mp_imag = -M_n_mp_imag
                    end
                    val3_real += M_n_mp_real * ζ_real - M_n_mp_imag * ζ_imag
                    val3_imag += M_n_mp_real * ζ_imag + M_n_mp_imag * ζ_real
                end

                # decrement index
                i_ζ_mp -= 1
            end

            for mp in -m:0 # -m <= mp <= 0
                #println("\tm'=$mp")
                i1, i2 = minmax(mp,m)
                T_mp_m = Ts[i_T + T_index(i1,i2)]
                # T_mp_m = T[T_index(mp,m)]
                #@show T_index(mp,m)

                # ζ_n^{m',m} = β_n^{m'} / β_n^m
                ζ_real, ζ_imag = ζ_sign(ζs_mag[i_ζ_mp] * T_mp_m, mp, m)

                # current weights
                i_weight = harmonic_index(n, -mp)

                # sum contribution
                if SP
                    # (-1)^mp * complex conjugate for -mp
                    M_n_mp_real, M_n_mp_imag = source_weights[1,1,i_weight], -source_weights[2,1,i_weight]
                    if isodd(mp)
                        M_n_mp_real = -M_n_mp_real
                        M_n_mp_imag = -M_n_mp_imag
                    end

                    val1_real += M_n_mp_real * ζ_real - M_n_mp_imag * ζ_imag
                    val1_imag += M_n_mp_real * ζ_imag + M_n_mp_imag * ζ_real
                end
                if VP
                    # (-1)^mp * complex conjugate for -mp
                    M_n_mp_real, M_n_mp_imag = source_weights[1,2,i_weight], -source_weights[2,2,i_weight]
                    if isodd(mp)
                        M_n_mp_real = -M_n_mp_real
                        M_n_mp_imag = -M_n_mp_imag
                    end
                    val2_real += M_n_mp_real * ζ_real - M_n_mp_imag * ζ_imag
                    val2_imag += M_n_mp_real * ζ_imag + M_n_mp_imag * ζ_real

                    # (-1)^mp * complex conjugate for -mp
                    M_n_mp_real, M_n_mp_imag = source_weights[1,3,i_weight], -source_weights[2,3,i_weight]
                    if isodd(mp)
                        M_n_mp_real = -M_n_mp_real
                        M_n_mp_imag = -M_n_mp_imag
                    end
                    val3_real += M_n_mp_real * ζ_real - M_n_mp_imag * ζ_imag
                    val3_imag += M_n_mp_real * ζ_imag + M_n_mp_imag * ζ_real
                end

                # decrement index
                i_ζ_mp -= 1
            end

            # augment index back to m'=1
            i_ζ_mp += 2

            #for mp in 1:m # 0 < mp <= m
            for mp in 1:n # mp > 0
                # T_mp_m = T[T_index(mp,m)]
                i1, i2 = minmax(mp,m)
                T_mp_m = Ts[i_T + T_index(i1,i2)]
                #@show T_index(mp,m)

                # ζ_n^{m',m} = β_n^{m'} / β_n^m
                ζ_real, ζ_imag = ζ_sign(ζs_mag[i_ζ_mp] * T_mp_m, mp, m)

                # current weights
                i_weight = harmonic_index(n, mp)

                # sum contribution
                if SP
                    M_n_mp_real, M_n_mp_imag = source_weights[1,1,i_weight], source_weights[2,1,i_weight]
                    val1_real += M_n_mp_real * ζ_real - M_n_mp_imag * ζ_imag
                    val1_imag += M_n_mp_real * ζ_imag + M_n_mp_imag * ζ_real
                end
                if VP
                    M_n_mp_real, M_n_mp_imag = source_weights[1,2,i_weight], source_weights[2,2,i_weight]
                    val2_real += M_n_mp_real * ζ_real - M_n_mp_imag * ζ_imag
                    val2_imag += M_n_mp_real * ζ_imag + M_n_mp_imag * ζ_real

                    M_n_mp_real, M_n_mp_imag = source_weights[1,3,i_weight], source_weights[2,3,i_weight]
                    val3_real += M_n_mp_real * ζ_real - M_n_mp_imag * ζ_imag
                    val3_imag += M_n_mp_real * ζ_imag + M_n_mp_imag * ζ_real
                end

                # increment index
                i_ζ_mp += 1
            end

            # update rotated_weights
            i = harmonic_index(n,m)
            if SP
                rotated_weights[1,1,i] += val1_real
                rotated_weights[2,1,i] += val1_imag
            end
            if VP
                rotated_weights[1,2,i] += val2_real
                rotated_weights[2,2,i] += val2_imag
                rotated_weights[1,3,i] += val3_real
                rotated_weights[2,3,i] += val3_imag
            end
        end

        # next T matrix
        np1 = n+1
        i_T += np1*np1
    end
end

"""
Rotate solid harmonic weights about the y axis by θ. Note that Hs_π2 and ζs_mag must be updated a priori, but Ts is updated en situ. Resets rotated_weights before computing.
"""
function rotate_multipole_y!(rotated_weights, source_weights, Ts, Hs_π2, ζs_mag, θ, expansion_order, expansion_switch)
    # get y-axis Wigner rotation matrix
    update_Ts!(Ts, Hs_π2, θ, expansion_order)

    # perform rotation
    _rotate_multipole_y!(rotated_weights, source_weights, Ts, ζs_mag, expansion_order, expansion_switch)
end

"""
Assumes Ts, Hs_π2, and ζs_mag have all been precomputed. Resets target_weights.
"""
function back_rotate_multipole_y!(target_weights, rotated_weights, Ts, ζs_mag, expansion_order, expansion_switch)
    _rotate_multipole_y!(target_weights, rotated_weights, Ts, ζs_mag, expansion_order, expansion_switch)
end

#------- LOCAL ROTATIONS -------#

function length_ηs(expansion_order)
    p = expansion_order
    return div((p+1)*(p+2)*(2p+3),6)
end

function update_ηs_mag!(ηs_mag, expansion_order::Val{P}) where P
    l = length(ηs_mag)
    l_desired = length_ηs(P)

    if l_desired > l # we need to compute more matrices
        # determine how many we already have
        n_already, l_already = 0, 0
        for n in 0:P
            l_already += (n+1)^2
            l_already == l && (n_already = n)
        end

        # calculate the rest
        resize!(ηs_mag, l_desired)
        update_ηs_mag!(ηs_mag, n_already+1, P)

    end # otherwise, we already have enough, so do nothing
end

@inline function η_sign(magnitude, mp, m)
    mod = (abs(m) - abs(mp)) % 4
    mod < 0 && (mod += 4)
    if mod == 0
        return magnitude, 0.0
    elseif mod == 1
        return 0.0, magnitude
    elseif mod == 2
        return -magnitude, 0.0
    elseif mod == 3
        return 0.0, -magnitude
    end
end

function update_ηs_mag!(ηs, n_start, n_end)
    if n_start == 0
        ηs[1] = 1.0
        n_start = 1
    end
    for n in n_start:n_end
        η = get_η(ηs, n)
        i = 1
        mag_m = 1.0 # magnitude of η
        for m in 0:n
            mag_mp = mag_m

            for mp in 0:n
                # store magnitude (we'll compute the sign on the fly)
                η[i] = mag_mp

                # recurse over mp
                num = mp+1 < n ? n - mp : 1
                denom = n + mp + 1
                mag_mp *= sqrt(num/denom)

                # augment index
                i += 1
            end

            # recurse over m
            num = n+m+1
            denom = m < n - 1 ? n - m : 1
            mag_m *= sqrt(num/denom)
        end
    end
end

function get_η(ηs, n)
    i1 = length_ηs(n-1)+1
    i2 = length_ηs(n)
    return view(ηs, i1:i2)
end

#--- functions to actually rotate ---#

function _rotate_local_y!(rotated_weights, source_weights, Ts, Hs_π2, ηs_mag, ::Val{P}, ::ExpansionSwitch{SP,VP}) where {P,SP,VP}
    # reset container
    rotated_weights .= zero(eltype(rotated_weights))

    i_η = 0
    i_T = 0
    # rotate each order
    @inbounds for n in 0:P

        for m in 0:n
            if SP
                val1_real = zero(eltype(Ts))
                val1_imag = zero(eltype(Ts))
            end
            if VP
                val2_real = zero(eltype(Ts))
                val2_imag = zero(eltype(Ts))
                val3_real = zero(eltype(Ts))
                val3_imag = zero(eltype(Ts))
            end

            i_η += n+1
            i_η_mp = i_η

            for mp in -n:-m-1 # -n <= mp < m
                # T_n^{m',m} = T_n^{m,m'} = T_n^{-m,-m'}
                T_mp_m = Ts[i_T + T_index(-m,-mp)] # should always be real

                # η_n^{m',m} = β_n^{m'} / β_n^m
                η_real, η_imag = η_sign(ηs_mag[i_η_mp] * T_mp_m, mp, m)

                # current weights
                i_weight = harmonic_index(n, -mp)

                # sum contribution
                if SP
                    # (-1)^mp * complex conjugate for -mp
                    M_n_mp_real, M_n_mp_imag = source_weights[1,1,i_weight], -source_weights[2,1,i_weight]
                    if isodd(mp)
                        M_n_mp_real = -M_n_mp_real
                        M_n_mp_imag = -M_n_mp_imag
                    end
                    val1_real += M_n_mp_real * η_real - M_n_mp_imag * η_imag
                    val1_imag += M_n_mp_real * η_imag + M_n_mp_imag * η_real
                end
                if VP
                    # (-1)^mp * complex conjugate for -mp
                    M_n_mp_real, M_n_mp_imag = source_weights[1,2,i_weight], -source_weights[2,2,i_weight]
                    if isodd(mp)
                        M_n_mp_real = -M_n_mp_real
                        M_n_mp_imag = -M_n_mp_imag
                    end
                    val2_real += M_n_mp_real * η_real - M_n_mp_imag * η_imag
                    val2_imag += M_n_mp_real * η_imag + M_n_mp_imag * η_real

                    # (-1)^mp * complex conjugate for -mp
                    M_n_mp_real, M_n_mp_imag = source_weights[1,3,i_weight], -source_weights[2,3,i_weight]
                    if isodd(mp)
                        M_n_mp_real = -M_n_mp_real
                        M_n_mp_imag = -M_n_mp_imag
                    end
                    val3_real += M_n_mp_real * η_real - M_n_mp_imag * η_imag
                    val3_imag += M_n_mp_real * η_imag + M_n_mp_imag * η_real
                end


                # decrement index
                i_η_mp -= 1
            end

            for mp in -m:0 # -m <= mp <= 0
                T_mp_m = Ts[i_T + T_index(mp,m)]

                # η_n^{m',m} = β_n^{m'} / β_n^m
                η_real, η_imag = η_sign(ηs_mag[i_η_mp] * T_mp_m, mp, m)

                # current weights
                i_weight = harmonic_index(n, -mp)

                # sum contribution
                if SP
                    # (-1)^mp * complex conjugate for -mp
                    M_n_mp_real, M_n_mp_imag = source_weights[1,1,i_weight], -source_weights[2,1,i_weight]
                    if isodd(mp)
                        M_n_mp_real = -M_n_mp_real
                        M_n_mp_imag = -M_n_mp_imag
                    end
                    val1_real += M_n_mp_real * η_real - M_n_mp_imag * η_imag
                    val1_imag += M_n_mp_real * η_imag + M_n_mp_imag * η_real
                end
                if VP
                    # (-1)^mp * complex conjugate for -mp
                    M_n_mp_real, M_n_mp_imag = source_weights[1,2,i_weight], -source_weights[2,2,i_weight]
                    if isodd(mp)
                        M_n_mp_real = -M_n_mp_real
                        M_n_mp_imag = -M_n_mp_imag
                    end
                    val2_real += M_n_mp_real * η_real - M_n_mp_imag * η_imag
                    val2_imag += M_n_mp_real * η_imag + M_n_mp_imag * η_real

                    # (-1)^mp * complex conjugate for -mp
                    M_n_mp_real, M_n_mp_imag = source_weights[1,3,i_weight], -source_weights[2,3,i_weight]
                    if isodd(mp)
                        M_n_mp_real = -M_n_mp_real
                        M_n_mp_imag = -M_n_mp_imag
                    end
                    val3_real += M_n_mp_real * η_real - M_n_mp_imag * η_imag
                    val3_imag += M_n_mp_real * η_imag + M_n_mp_imag * η_real
                end

                # decrement index
                i_η_mp -= 1
            end

            # augment index back to m'=1
            i_η_mp += 2

            for mp in 1:n # m < mp <= n
                # T_n^{m',m} = T_n^{m,m'}
                i1, i2 = minmax(mp,m)
                T_mp_m = Ts[i_T + T_index(i1,i2)]

                # η_n^{m',m} = β_n^{m'} / β_n^m
                η_real, η_imag = η_sign(ηs_mag[i_η_mp] * T_mp_m, mp, m)

                # current weights
                i_weight = harmonic_index(n, mp)

                # sum contribution
                if SP
                    M_n_mp_real, M_n_mp_imag = source_weights[1,1,i_weight], source_weights[2,1,i_weight]
                    val1_real += M_n_mp_real * η_real - M_n_mp_imag * η_imag
                    val1_imag += M_n_mp_real * η_imag + M_n_mp_imag * η_real
                end
                if VP
                    M_n_mp_real, M_n_mp_imag = source_weights[1,2,i_weight], source_weights[2,2,i_weight]
                    val2_real += M_n_mp_real * η_real - M_n_mp_imag * η_imag
                    val2_imag += M_n_mp_real * η_imag + M_n_mp_imag * η_real

                    M_n_mp_real, M_n_mp_imag = source_weights[1,3,i_weight], source_weights[2,3,i_weight]
                    val3_real += M_n_mp_real * η_real - M_n_mp_imag * η_imag
                    val3_imag += M_n_mp_real * η_imag + M_n_mp_imag * η_real
                end

                # increment index
                i_η_mp += 1
            end

            # update rotated_weights
            i = harmonic_index(n,m)
            if SP
                rotated_weights[1,1,i] += val1_real
                rotated_weights[2,1,i] += val1_imag
            end
            if VP
                rotated_weights[1,2,i] += val2_real
                rotated_weights[2,2,i] += val2_imag
                rotated_weights[1,3,i] += val3_real
                rotated_weights[2,3,i] += val3_imag
            end
        end

        # next T matrix
        np1 = n+1
        i_T += np1*np1
    end
end

function rotate_local_y!(rotated_weights, source_weights, Ts, Hs_π2, ηs_mag, θ, expansion_order, expansion_switch)
    # get y-axis Wigner rotation matrix
    update_Ts!(Ts, Hs_π2, θ, expansion_order)

    # perform rotation
    _rotate_local_y!(rotated_weights, source_weights, Ts, Hs_π2, ηs_mag, expansion_order, expansion_switch)
end

"""
Assumes Ts, Hs_π2, and ηs_mag have all been precomputed. Resets target_weights.
"""
function back_rotate_local_y!(target_weights, rotated_weights, Ts, Hs_π2, ηs_mag, expansion_order, expansion_switch)
    _rotate_local_y!(target_weights, rotated_weights, Ts, Hs_π2, ηs_mag, expansion_order, expansion_switch)
end

