function evaluate_local!(systems, branch::MultiBranch, derivatives_switches, expansion_order)
    for i in eachindex(systems)
        evaluate_local!(systems[i], branch.bodies_index[i], branch.local_expansion, branch.center, derivatives_switches[i], expansion_order)
    end
end

function evaluate_local!(system, branch::SingleBranch, derivatives_switch, expansion_order)
    evaluate_local!(system, branch.bodies_index, branch.local_expansion, branch.center, derivatives_switch, expansion_order)
end

function evaluate_local!(system, bodies_index, local_expansion, expansion_center, derivatives_switch::DerivativesSwitch{PS,VPS,VS,GS}, expansion_order::Val{P}) where {PS,VPS,VS,GS,P}
    for i_body in bodies_index
        scalar_potential, velocity, gradient = evaluate_local(system[i_body,POSITION] - expansion_center, local_expansion, derivatives_switch, expansion_order)
        PS && (system[i_body, SCALAR_POTENTIAL] += scalar_potential)
        if VS
            vpx, vpy, vpz = system[i_body, VELOCITY]
            system[i_body, VELOCITY] += SVector{3}(velocity[1]+vpx, velocity[2]+vpy, velocity[3]+vpz)
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

function evaluate_local(Δx, harmonics, local_expansion, derivatives_switch, expansion_order::Val{P}) where P
    # convert to spherical coordinates
    r, θ, ϕ = cartesian_to_spherical(Δx)

    # expansion basis is the regular solid harmonics
    regular_harmonics!(harmonics, r, θ, ϕ, expansion_order)

    # multiply by weights
    val = zero(eltype(local_expansion))
    i = 1
    for n in 0:P
        # m=0
        Rnm_real, Rnm_imag = harmonics[1,i], harmonics[2,i]
        Lnm_real, Lnm_imag = local_expansion[1,1,i], local_expansion[2,1,i]
        val += Rnm_real * Lnm_real - Rnm_imag * Lnm_imag
        i += 1

        # m>0
        for m in 1:n
            Rnm_real, Rnm_imag = harmonics[1,i], harmonics[2,i]
            Lnm_real, Lnm_imag = local_expansion[1,1,i], local_expansion[2,1,i]
            val += 2 * (Rnm_real * Lnm_real - Rnm_imag * Lnm_imag)
            i += 1
        end
    end

    return val, zero(SVector{3,eltype(local_expansion)}), zero(SMatrix{3,3,eltype(local_expansion),9})
end
