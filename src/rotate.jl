function length_Hs(expansion_order)
    p = expansion_order
    return div((p+1)*(p+2)*(2p+3), 6)
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
    return m^2 + mp + m + 1
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
    dnms = zeros(2*n_end+1)

    # Hs
    for n in n_start:n_end
        H = get_H(Hs, n)

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
        for m in -n:-1
            dnms[n+m+1] = -sqrt((n-m)*(n+m+1)) / 2
        end
        for m in 0:n
            dnms[n+m+1] = sqrt((n-m)*(n+m+1)) / 2
        end

        # H_n_mp+1_m
        for mp in 1:n-1
            for m in mp+1:n-1
                val = dnms[n+mp] * H[H_index(mp-1,m)] - dnms[n+m] * H[H_index(mp,m-1)] + dnms[n+m+1] * H[H_index(mp,m+1)]
                H[H_index(mp+1,m)] = val / dnms[n+mp+1]
            end
            m = n
            val = dnms[n+mp] * H[H_index(mp-1,m)] - dnms[n+m] * H[H_index(mp,m-1)]
            H[H_index(mp+1,m)] = val / dnms[n+mp+1]
        end

        # H_n_mp-1_m
        for mp in 0:-1:-n+1
            for m in -mp+1:n-1
                val = dnms[n+mp+1] * H[H_index(mp+1,m)] + dnms[n+m] * H[H_index(mp,m-1)] - dnms[n+m+1] * H[H_index(mp,m+1)]
                H[H_index(mp-1,m)] = val / dnms[n+mp]
                @show mp-1
            end
            m = n
            val = dnms[n+mp+1] * H[H_index(mp+1,m)] + dnms[n+m] * H[H_index(mp,m-1)]
            H[H_index(mp-1,m)] = val / dnms[n+mp]
            @show mp-1
        end

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

#update_Hs_π2!(Hs_π2, 20)
