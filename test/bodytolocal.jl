function body_to_local_point!(::Type{Point{Source}}, multipole_coefficients, harmonics, Δx, strength, expansion_order::Val{P}) where P

    # invert the sign of the strength so v=∇ϕ instead of v=-∇ϕ
    # this will allow induced velocities to be added to those induced by vortex elements
    # we'll invert the potential again at the end
    strength = -strength

    # irregular harmonics
    ρ, θ, ϕ = FastMultipole.cartesian_to_spherical(Δx)
    FastMultipole.irregular_harmonics!(harmonics, ρ, θ, ϕ, expansion_order)

    # update coefficients
    _1_n = 1.0
    i = 1
    for n in 0:P
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

