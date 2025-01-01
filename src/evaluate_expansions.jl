function evaluate_local!(systems, branch::MultiBranch, harmonics, velocity_n_m, expansion_order, lamb_helmholtz, derivatives_switches)
    for i in eachindex(systems)
        evaluate_local!(systems[i], branch.bodies_index[i], harmonics, velocity_n_m, branch.local_expansion, branch.target_center, expansion_order, lamb_helmholtz, derivatives_switches[i])
    end
end

function evaluate_local!(system, branch, expansion_order, lamb_helmholtz, derivatives_switch)
    harmonics = branch.harmonics
    velocity_n_m = zeros(2,3,size(branch.multipole_expansion,3))
    evaluate_local!(system, branch, harmonics, velocity_n_m, expansion_order, lamb_helmholtz, derivatives_switch)
end

function evaluate_local!(system, branch::SingleBranch, harmonics, velocity_n_m, expansion_order, lamb_helmholtz, derivatives_switch)
    evaluate_local!(system, branch.bodies_index, harmonics, velocity_n_m, branch.local_expansion, branch.target_center, expansion_order, lamb_helmholtz, derivatives_switch)
end

function evaluate_local!(system, bodies_index, harmonics, velocity_n_m, local_expansion, expansion_center, expansion_order::Val{P}, lamb_helmholtz, derivatives_switch::DerivativesSwitch{PS,VS,GS}) where {P,PS,VS,GS}
    for i_body in bodies_index
        scalar_potential, velocity, gradient = evaluate_local(system[i_body,POSITION] - expansion_center, harmonics, velocity_n_m, local_expansion, expansion_order, lamb_helmholtz, derivatives_switch)
        PS && (system[i_body, SCALAR_POTENTIAL] += scalar_potential)
        if VS
            vpx, vpy, vpz = system[i_body, VELOCITY]
            system[i_body, VELOCITY] = SVector{3}(velocity[1]+vpx, velocity[2]+vpy, velocity[3]+vpz)
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

function evaluate_local(Δx, harmonics, velocity_n_m, local_expansion, expansion_order::Val{P}, ::Val{LH}, ::DerivativesSwitch{PS,VS,GS}) where {P,LH,PS,VS,GS}
    # convert to spherical coordinates
    r, θ, ϕ = cartesian_to_spherical(Δx)

    # expansion basis is the regular solid harmonics
    regular_harmonics!(harmonics, r, θ, ϕ, expansion_order)

    #--- declare/reset variables ---#

    # scalar potential
    u = zero(eltype(local_expansion))

    # velocity
    vx, vy, vz = zero(eltype(local_expansion)), zero(eltype(local_expansion)), zero(eltype(local_expansion))
    GS && ( velocity_n_m .= zero(eltype(local_expansion)) )

    # velocity gradient
    vxx, vxy, vxz = zero(eltype(local_expansion)), zero(eltype(local_expansion)), zero(eltype(local_expansion))
    vyx, vyy, vyz = zero(eltype(local_expansion)), zero(eltype(local_expansion)), zero(eltype(local_expansion))
    vzx, vzy, vzz = zero(eltype(local_expansion)), zero(eltype(local_expansion)), zero(eltype(local_expansion))

    # index
    i_n_m = 0

    #------- n = 0, m = 0 -------#

    n = 0

    # update index
    i_n_m += 1

    # get regular harmonic
    Rnm_real, Rnm_imag = harmonics[1,1,i_n_m], harmonics[2,1,i_n_m]

    # scalar potential
    if PS && !LH # scalar potential is transformed (and nonsensical) to preserve velocity when using Lamb-Helmholtz
        ϕ_n_m_real = local_expansion[1,1,i_n_m]
        ϕ_n_m_imag = local_expansion[2,1,i_n_m]
        u += Rnm_real * ϕ_n_m_real - Rnm_imag * ϕ_n_m_imag
    end

    # velocity
    if VS || GS
        #=
        vx = -Im[ϕ_{1}^{1}]
        vy = -Re[ϕ_{1}^{1}]
        vz = -ϕ_{1}^0
        =#

        # expansion coefficients
        vx_n_m_real = zero(eltype(local_expansion))
        vx_n_m_imag = zero(eltype(local_expansion))
        vy_n_m_real = zero(eltype(local_expansion))
        vy_n_m_imag = zero(eltype(local_expansion))
        vz_n_m_real = zero(eltype(local_expansion))
        vz_n_m_imag = zero(eltype(local_expansion))

        # due to ϕ
        ϕ_np1_m_real = local_expansion[1,1,i_n_m+n+1]
        ϕ_np1_m_imag = local_expansion[2,1,i_n_m+n+1]
        ϕ_np1_mp1_real = local_expansion[1,1,i_n_m+n+2]
        ϕ_np1_mp1_imag = local_expansion[2,1,i_n_m+n+2]

        vx_n_m_real -= ϕ_np1_mp1_imag
        vy_n_m_real -= ϕ_np1_mp1_real
        vz_n_m_real -= ϕ_np1_m_real
        vz_n_m_imag -= ϕ_np1_m_imag

        if VS
            # evaluate expansion
            vx += vx_n_m_real * Rnm_real
            vy += vy_n_m_real * Rnm_real
            vz += vz_n_m_real * Rnm_real - vz_n_m_imag * Rnm_imag
        end

        if GS
            # store components for computing velocity gradient
            velocity_n_m[1,1,i_n_m] = vx_n_m_real # x component
            velocity_n_m[1,2,i_n_m] = vy_n_m_real # y component
            velocity_n_m[1,3,i_n_m] = vz_n_m_real # z component
            velocity_n_m[2,3,i_n_m] = vz_n_m_imag # z component
        end

    end

    #------- n > 0 -------#

    for n in 1:P

        #--- m = 0 ---#

        # update index
        i_n_m += 1

        # get regular harmonic
        Rnm_real, Rnm_imag = harmonics[1,1,i_n_m], harmonics[2,1,i_n_m]

        # scalar potential
        if PS && !LH
            ϕ_n_m_real = local_expansion[1,1,i_n_m]
            ϕ_n_m_imag = local_expansion[2,1,i_n_m]
            u += Rnm_real * ϕ_n_m_real - Rnm_imag * ϕ_n_m_imag
        end

        # velocity
        if VS || GS
            #=
            vx = (Re[χ_1^1] - Im[ϕ_2^1]) R_n^m
            vy = -(Re[ϕ_2^1] + Im[χ_1^1]) R_n^m
            vz = -ϕ_2^0 R_n^m
            =#

            # expansion coefficients
            vx_n_m_real = zero(eltype(local_expansion))
            vx_n_m_imag = zero(eltype(local_expansion))
            vy_n_m_real = zero(eltype(local_expansion))
            vy_n_m_imag = zero(eltype(local_expansion))
            vz_n_m_real = zero(eltype(local_expansion))
            vz_n_m_imag = zero(eltype(local_expansion))

            # due to ϕ
            if n < P
                ϕ_np1_m_real = local_expansion[1,1,i_n_m+n+1]
                ϕ_np1_m_imag = local_expansion[2,1,i_n_m+n+1]
                ϕ_np1_mp1_real = local_expansion[1,1,i_n_m+n+2]
                ϕ_np1_mp1_imag = local_expansion[2,1,i_n_m+n+2]

                vx_n_m_real -= ϕ_np1_mp1_imag
                vy_n_m_real -= ϕ_np1_mp1_real
                vz_n_m_real -= ϕ_np1_m_real
                vz_n_m_imag -= ϕ_np1_m_imag
            end


            # due to χ
            if LH
                χ_n_mp1_real = local_expansion[1,2,i_n_m+1]
                χ_n_mp1_imag = local_expansion[2,2,i_n_m+1]

                vx_n_m_real += n * χ_n_mp1_real
                vy_n_m_real -= n * χ_n_mp1_imag
            end

            if VS
                # evaluate expansion
                vx += vx_n_m_real * Rnm_real
                vy += vy_n_m_real * Rnm_real
                vz += vz_n_m_real * Rnm_real - vz_n_m_imag * Rnm_imag
            end

            if GS
                # store components for computing velocity gradient
                velocity_n_m[1,1,i_n_m] = vx_n_m_real # x component
                velocity_n_m[1,2,i_n_m] = vy_n_m_real # y component
                velocity_n_m[1,3,i_n_m] = vz_n_m_real # z component
                velocity_n_m[2,3,i_n_m] = vz_n_m_imag # z component
            end

        end

        #--- m > 0 ---#

        for m in 1:n

            # update index
            i_n_m += 1

            # get regular harmonic
            Rnm_real, Rnm_imag = harmonics[1,1,i_n_m], harmonics[2,1,i_n_m]

            # scalar potential
            if PS && !LH
                ϕ_n_m_real = local_expansion[1,1,i_n_m]
                ϕ_n_m_imag = local_expansion[2,1,i_n_m]
                u += 2 * (Rnm_real * ϕ_n_m_real - Rnm_imag * ϕ_n_m_imag)
            end

            # velocity
            if VS || GS
                #=
                vx = (im * (ϕ_{n+1}^{m-1} + ϕ_{n+1}^{m+1}) + (n-m) χ_n^{m+1} - (n+m) χ_n^{m-1}) / 2
                vy = (ϕ_{n+1}^{m-1} - ϕ_{n+1}^{m+1} + im * (n-m) * χ_n^{m+1} + im * (n+m) χ_n^{m-1}) / 2
                vz = -ϕ_{n+1}^m - im * m * χ_n^m
                =#

                # expansion coefficients
                vx_n_m_real = zero(eltype(local_expansion))
                vx_n_m_imag = zero(eltype(local_expansion))
                vy_n_m_real = zero(eltype(local_expansion))
                vy_n_m_imag = zero(eltype(local_expansion))
                vz_n_m_real = zero(eltype(local_expansion))
                vz_n_m_imag = zero(eltype(local_expansion))

                if n < P
                    ϕ_np1_mm1_real = local_expansion[1,1,i_n_m+n]
                    ϕ_np1_mm1_imag = local_expansion[2,1,i_n_m+n]
                    ϕ_np1_m_real = local_expansion[1,1,i_n_m+n+1]
                    ϕ_np1_m_imag = local_expansion[2,1,i_n_m+n+1]
                    ϕ_np1_mp1_real = local_expansion[1,1,i_n_m+n+2]
                    ϕ_np1_mp1_imag = local_expansion[2,1,i_n_m+n+2]

                    vx_n_m_real -= (ϕ_np1_mm1_imag + ϕ_np1_mp1_imag) * 0.5
                    vx_n_m_imag += (ϕ_np1_mm1_real + ϕ_np1_mp1_real) * 0.5
                    vy_n_m_real += (ϕ_np1_mm1_real - ϕ_np1_mp1_real) * 0.5
                    vy_n_m_imag += (ϕ_np1_mm1_imag - ϕ_np1_mp1_imag) * 0.5
                    vz_n_m_real -= ϕ_np1_m_real
                    vz_n_m_imag -= ϕ_np1_m_imag
                end

                if LH
                    # extract expansion coefficients
                    χ_n_mm1_real = local_expansion[1,2,i_n_m-1]
                    χ_n_mm1_imag = local_expansion[2,2,i_n_m-1]

                    # form velocity coefficients
                    vx_n_m_real -= (n+m) * χ_n_mm1_real * 0.5
                    vx_n_m_imag -= (n+m) * χ_n_mm1_imag * 0.5
                    vy_n_m_real -= (n+m) * χ_n_mm1_imag * 0.5
                    vy_n_m_imag += (n+m) * χ_n_mm1_real * 0.5
                    if m < n
                        χ_n_mp1_real = local_expansion[1,2,i_n_m+1]
                        χ_n_mp1_imag = local_expansion[2,2,i_n_m+1]

                        vx_n_m_real += (n-m) * χ_n_mp1_real * 0.5
                        vx_n_m_imag += (n-m) * χ_n_mp1_imag * 0.5
                        vy_n_m_real -= (n-m) * χ_n_mp1_imag * 0.5
                        vy_n_m_imag += (n-m) * χ_n_mp1_real * 0.5
                    end

                    # extract expansion coefficients
                    χ_n_m_real = local_expansion[1,2,i_n_m]
                    χ_n_m_imag = local_expansion[2,2,i_n_m]

                    # form velocity coefficients
                    vz_n_m_real += m * χ_n_m_imag
                    vz_n_m_imag -= m * χ_n_m_real
                end

                # evaluate expansion
                vx += 2 * (vx_n_m_real * Rnm_real - vx_n_m_imag * Rnm_imag)
                vy += 2 * (vy_n_m_real * Rnm_real - vy_n_m_imag * Rnm_imag)
                vz += 2 * (vz_n_m_real * Rnm_real - vz_n_m_imag * Rnm_imag)

                if GS
                    # store components for computing velocity gradient
                    velocity_n_m[1,1,i_n_m] = vx_n_m_real # x component
                    velocity_n_m[2,1,i_n_m] = vx_n_m_imag # x component
                    velocity_n_m[1,2,i_n_m] = vy_n_m_real # y component
                    velocity_n_m[2,2,i_n_m] = vy_n_m_imag # y component
                    velocity_n_m[1,3,i_n_m] = vz_n_m_real # z component
                    velocity_n_m[2,3,i_n_m] = vz_n_m_imag # z component
                end

            end
        end
    end

    if GS

        # index
        i_n_m = 0

        for n in 0:P-1

            #--- m = 0 ---#

            # update index
            i_n_m += 1

            # get regular harmonic
            Rnm_real, Rnm_imag = harmonics[1,1,i_n_m], harmonics[2,1,i_n_m]

            vg_xx_real = -velocity_n_m[2,1,i_n_m+n+2]
            vxx += vg_xx_real * Rnm_real

            vg_yx_real = -velocity_n_m[1,1,i_n_m+n+2]
            vyx += vg_yx_real * Rnm_real

            vg_zx_real = -velocity_n_m[1,1,i_n_m+n+1]
            vg_zx_imag = -velocity_n_m[2,1,i_n_m+n+1]
            vzx += vg_zx_real * Rnm_real - vg_zx_imag * Rnm_imag

            vg_xy_real = -velocity_n_m[2,2,i_n_m+n+2]
            vxy += vg_xy_real * Rnm_real

            vg_yy_real = -velocity_n_m[1,2,i_n_m+n+2]
            vyy += vg_yy_real * Rnm_real

            vg_zy_real = -velocity_n_m[1,2,i_n_m+n+1]
            vg_zy_imag = -velocity_n_m[2,2,i_n_m+n+1]
            vzy += vg_zy_real * Rnm_real - vg_zy_imag * Rnm_imag

            vg_xz_real = -velocity_n_m[2,3,i_n_m+n+2]
            vxz += vg_xz_real * Rnm_real

            vg_yz_real = -velocity_n_m[1,3,i_n_m+n+2]
            vyz += vg_yz_real * Rnm_real

            vg_zz_real = -velocity_n_m[1,3,i_n_m+n+1]
            vg_zz_imag = -velocity_n_m[2,3,i_n_m+n+1]
            vzz += vg_zz_real * Rnm_real - vg_zz_imag * Rnm_imag


            #--- m > 0 ---#

            for m in 1:n

                # update index
                i_n_m += 1

                # get regular harmonic
                Rnm_real, Rnm_imag = harmonics[1,1,i_n_m], harmonics[2,1,i_n_m]

                vg_xx_real = -(velocity_n_m[2,1,i_n_m+n] + velocity_n_m[2,1,i_n_m+n+2]) * 0.5
                vg_xx_imag = (velocity_n_m[1,1,i_n_m+n] + velocity_n_m[1,1,i_n_m+n+2]) * 0.5
                vxx += 2 * (vg_xx_real * Rnm_real - vg_xx_imag * Rnm_imag)

                vg_yx_real = (velocity_n_m[1,1,i_n_m+n] - velocity_n_m[1,1,i_n_m+n+2]) * 0.5
                vg_yx_imag = (velocity_n_m[2,1,i_n_m+n] - velocity_n_m[2,1,i_n_m+n+2]) * 0.5
                vyx += 2 * (vg_yx_real * Rnm_real - vg_yx_imag * Rnm_imag)

                vg_zx_real = -velocity_n_m[1,1,i_n_m+n+1]
                vg_zx_imag = -velocity_n_m[2,1,i_n_m+n+1]
                vzx += 2 * (vg_zx_real * Rnm_real - vg_zx_imag * Rnm_imag)

                vg_xy_real = -(velocity_n_m[2,2,i_n_m+n] + velocity_n_m[2,2,i_n_m+n+2]) * 0.5
                vg_xy_imag = (velocity_n_m[1,2,i_n_m+n] + velocity_n_m[1,2,i_n_m+n+2]) * 0.5
                vxy += 2 * (vg_xy_real * Rnm_real - vg_xy_imag * Rnm_imag)

                vg_yy_real = (velocity_n_m[1,2,i_n_m+n] - velocity_n_m[1,2,i_n_m+n+2]) * 0.5
                vg_yy_imag = (velocity_n_m[2,2,i_n_m+n] - velocity_n_m[2,2,i_n_m+n+2]) * 0.5
                vyy += 2 * (vg_yy_real * Rnm_real - vg_yy_imag * Rnm_imag)

                vg_zy_real = -velocity_n_m[1,2,i_n_m+n+1]
                vg_zy_imag = -velocity_n_m[2,2,i_n_m+n+1]
                vzy += 2 * (vg_zy_real * Rnm_real - vg_zy_imag * Rnm_imag)

                vg_xz_real = -(velocity_n_m[2,3,i_n_m+n] + velocity_n_m[2,3,i_n_m+n+2]) * 0.5
                vg_xz_imag = (velocity_n_m[1,3,i_n_m+n] + velocity_n_m[1,3,i_n_m+n+2]) * 0.5
                vxz += 2 * (vg_xz_real * Rnm_real - vg_xz_imag * Rnm_imag)

                vg_yz_real = (velocity_n_m[1,3,i_n_m+n] - velocity_n_m[1,3,i_n_m+n+2]) * 0.5
                vg_yz_imag = (velocity_n_m[2,3,i_n_m+n] - velocity_n_m[2,3,i_n_m+n+2]) * 0.5
                vyz += 2 * (vg_yz_real * Rnm_real - vg_yz_imag * Rnm_imag)

                vg_zz_real = -velocity_n_m[1,3,i_n_m+n+1]
                vg_zz_imag = -velocity_n_m[2,3,i_n_m+n+1]
                vzz += 2 * (vg_zz_real * Rnm_real - vg_zz_imag * Rnm_imag)

            end
        end
    end

    return -u * ONE_OVER_4π, SVector{3}(vx,vy,vz) * ONE_OVER_4π, SMatrix{3,3,eltype(local_expansion),9}(vxx, vxy, vxz, vyx, vyy, vyz, vzx, vzy, vzz) * ONE_OVER_4π
end

