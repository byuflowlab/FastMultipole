function evaluate_multipole!(system, bodies_index, multipole_weights, derivatives_switch::DerivativesSwitch{PS,VPS,VS,GS}, expansion_order, expansion_center) where {PS,VPS,VS,GS}
    for i_body in bodies_index
		scalar_potential, vector_potential, velocity, gradient = evaluate_multipole(system[i_body,POSITION], expansion_center, multipole_weights, derivatives_switch, expansion_order)

        if PS
            system[i_body,SCALAR_POTENTIAL] += scalar_potential
        end
        if VPS
            vpx, vpy, vpz = system[i_body,VECTOR_POTENTIAL]
            system[i_body,VECTOR_POTENTIAL] = SVector{3}(vpx+vector_potential[1],vpy+vector_potential[2],vpz+vector_potential[3])
        end
        if VS
            vpx, vpy, vpz = system[i_body,VELOCITY]
			system[i_body,VELOCITY] = SVector{3}(velocity[1]+vpx, velocity[2]+vpy, velocity[3]+vpz)
        end
        if GS
            v1, v2, v3, v4, v5, v6, v7, v8, v9 = system[i_body,VELOCITY_GRADIENT]
            system[i_body,VELOCITY_GRADIENT] = SMatrix{3,3}(
                gradient[1] + v1,
                gradient[2] + v2,
                gradient[3] + v3,
                gradient[4] + v4,
                gradient[5] + v5,
                gradient[6] + v6,
                gradient[7] + v7,
                gradient[8] + v8,
                gradient[9] + v9
            )
        end

    end
end

function evaluate_multipole(x_target, expansion_center, multipole_weights::AbstractArray{TF}, derivatives_switch::DerivativesSwitch{PS,VPS,VS,GS}, expansion_order::Val{P}) where {TF,PS,VPS,VS,GS,P}
    # outputs
    potential = zero(TF)
    velocity = zero(SVector{3,TF})
    velocity_gradient = zero(SMatrix{3,3,TF,9})

    # distance vector
    Δx = x_target - expansion_center
    ρ, θ, ϕ = cartesian_to_spherical(Δx)

    # compute irregular solid harmonics on the fly
    y, x = sincos(θ)
    fact = 1.0 # 2m+1 (odd integers)
    pn = 1 # Legendre polynomial of degree n, order n
    one_over_ρ = 1.0 / ρ # NOTE: this should never be singular, as we only evaluate irregular harmonics far away
    ρm = one_over_ρ # (-1)^m / ρ^(n+1)
    i_ei_imag, i_ei_real = sincos(ϕ+π/2) # i e^(iϕ) = e^[i(ϕ+π/2)]
    i_eim_real, i_eim_imag = 1.0, 0.0 # i^m e^(i * m * phi) = e^[i m (ϕ+π/2)]

    for m=0:P # n=m
        p = pn
        i = harmonic_index(m, m)
        ρm_p = ρm * p

        # irregular solid harmonics
        Snm_real = ρm_p * i_eim_real
        Snm_imag = ρm_p * i_eim_imag

        # multipole weights
        Mnm_real = multipole_weights[1,1,i]
        Mnm_imag = multipole_weights[2,1,i]

        # evaluate
        if PS
            Δ_potential = Mnm_real * Snm_real - Mnm_imag * Snm_imag
            m > 0 && (Δ_potential *= 2)
            potential += Δ_potential
        end

        p1 = p
        p = x * (2 * m + 1) * p1
        ρm *= -one_over_ρ
        ρn = -ρm

        for n=m+1:P # n>m
            i = harmonic_index(n, m)
            ρn_p = ρn * p

            # irregular solid harmonics
            Snm_real = ρn_p * i_eim_real
            Snm_imag = ρn_p * i_eim_imag

            # multipole weights
            Mnm_real = multipole_weights[1,1,i]
            Mnm_imag = multipole_weights[2,1,i]

            # evaluate
            if PS
                Δ_potential = Mnm_real * Snm_real - Mnm_imag * Snm_imag
                m > 0 && (Δ_potential *= 2)
                potential += Δ_potential
            end

            p2 = p1
            p1 = p
            p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1)
            ρn *= one_over_ρ * (n - m + 1)
        end

        # recurse
        pn *= -fact * y
        fact += 2
        i_eim_real_tmp = i_eim_real
        i_eim_imag_tmp = i_eim_imag
        i_eim_real = i_eim_real_tmp * i_ei_real - i_eim_imag_tmp * i_ei_imag
        i_eim_imag = i_eim_real_tmp * i_ei_imag + i_eim_imag_tmp * i_ei_real
    end

    return potential, velocity, velocity_gradient
end
