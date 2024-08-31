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
    return (m*(m+1))>>1 + mp + 1
end

@inline function Pnm_index(n, m)
    return ((n)*(n+1))>>1 + m + 1
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

function update_Hs_π2!(expansion_order)
    l = length(Hs_π2)
    l_desired = length_Hs(expansion_order)

    if l_desired > l # we need to compute more matrices
        # determine how many we already have
        n_already, l_already = 0, 0
        for n in 0:expansion_order
            l_already += (n+1)^2
            l_already == l && (n_already = n)
        end

        # calculate the rest
        β = pi/2
        resize!(Hs_π2, l_desired)
        compute_Hs!(Hs_π2, β, n_already+1, expansion_order)

    end # otherwise, we already have enough, so do nothing
end

#--- use precomputed H matrices to compute rotations for arbitrary angle ---#

function length_Ts(expansion_order)
    p = expansion_order
    return div((p+1)*(p+2)*(2p+3),6)
end

function length_T(n)
    return (n+1)^2
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
    return m^2 + mp + m + 1
end

@inline function get_scalar(m_mp)
    isodd(m_mp) && (m_mp += 1)
    if iseven(m_mp>>1)
        return 1.0
    else
        return -1.0
    end
end

function update_Ts!(Ts, β, expansion_order)
    p = expansion_order

    # initialize recursive indices
    i_H_start = length_Hs(0) + 1 # 2
    i_H_end = i_H_start + 2 # n=1 matrix contains 3 unique entries
    i_T_start = 2
    i_T_end = i_T_start + 3 # n=1 T matrix contains 4 unique entries

    # other recursive values
    sβ, cβ = sincos(β)
    _1_n = -1

    # rotate each expansion order
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
                    val = H_π2_n_mp_ν * H_π2_n_m_ν * sc_νβ_m_mp
                    val_positive_mp += val
                    val_negative_mp += val * _1_n_mp_ν * _1_mp_odd

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
                val_negative_mp += val * _1_n_mp * _1_mp_odd
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

