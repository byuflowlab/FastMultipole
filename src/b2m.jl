
#####
##### multipole generation convenience functions
#####
@inline function B2M!_vortexpoint(system, branch, bodies_index, harmonics, expansion_order::Val{P}) where P
    c_x, c_y, c_z = branch.center
    for i_body in bodies_index
        b_x, b_y, b_z = system[i_body,POSITION]
        x = b_x - c_x
        z = b_z - c_z
        y = b_y - c_y
        qx, qy, qz = system[i_body,STRENGTH]
        r, theta, phi = cartesian_2_spherical(x,y,z)
        regular_harmonic!(harmonics, r ,theta, -phi, P) # Ylm^* -> -dx[3]
        # update values
        for l in 0:P
            for m in 0:l
                i_solid_harmonic = l*l + l + m + 1
                i_compressed = 1 + (l * (l + 1)) >> 1 + m # only save half as Yl{-m} = conj(Ylm)
                #branch.multipole_expansion[2,i_compressed] += harmonics[i_solid_harmonic] * qx * ONE_OVER_4PI
                #branch.multipole_expansion[3,i_compressed] += harmonics[i_solid_harmonic] * qy * ONE_OVER_4PI
                #branch.multipole_expansion[4,i_compressed] += harmonics[i_solid_harmonic] * qz * ONE_OVER_4PI
                branch.multipole_expansion[1,2,i_compressed] += harmonics[1,i_solid_harmonic] * qx * ONE_OVER_4PI
                branch.multipole_expansion[2,2,i_compressed] += harmonics[2,i_solid_harmonic] * qx * ONE_OVER_4PI
                branch.multipole_expansion[1,3,i_compressed] += harmonics[1,i_solid_harmonic] * qy * ONE_OVER_4PI
                branch.multipole_expansion[2,3,i_compressed] += harmonics[2,i_solid_harmonic] * qy * ONE_OVER_4PI
                branch.multipole_expansion[1,4,i_compressed] += harmonics[1,i_solid_harmonic] * qz * ONE_OVER_4PI
                branch.multipole_expansion[2,4,i_compressed] += harmonics[2,i_solid_harmonic] * qz * ONE_OVER_4PI
            end
        end
    end
end

@inline function B2M!_sourcepoint(system, branch, bodies_index, harmonics, expansion_order::Val{P}) where P
    c_x, c_y, c_z = branch.center
    for i_body in bodies_index
        b_x, b_y, b_z = system[i_body,POSITION]
        x = b_x - c_x
        z = b_z - c_z
        y = b_y - c_y
        q = system[i_body,STRENGTH]
        branch.charge[] += abs(q)
        r, theta, phi = cartesian_2_spherical(x,y,z)
        regular_harmonic!(harmonics, r, theta, -phi, P) # Ylm^* -> -dx[3]
        # update values
        for l in 0:P
            for m in 0:l
                i_solid_harmonic = l*l + l + m + 1
                i_compressed = 1 + (l * (l + 1)) >> 1 + m # only save half as Yl{-m} = conj(Ylm)
                #branch.multipole_expansion[1,i_compressed] += harmonics[i_solid_harmonic] * q
                branch.multipole_expansion[1,1,i_compressed] += harmonics[1,i_solid_harmonic] * q
                branch.multipole_expansion[2,1,i_compressed] += harmonics[2,i_solid_harmonic] * q
            end
        end
    end
end

@inline function update_multipole_expansion_panel!(multipole_expansion, i_compressed, strength_J_over_4pi, inm_real, inm_imag, iam_real, iam_imag, inm_prev_mm1_real, inm_prev_mm1_imag, inm_prev_m_real, inm_prev_m_imag, inm_prev_mp1_real, inm_prev_mp1_imag, n, panel_type::UniformSourcePanel{sign}, negative_1_n) where sign
    multipole_expansion[1,1,i_compressed] += strength_J_over_4pi * complex_multiply_real(inm_real, inm_imag, iam_real, iam_imag) * sign
    multipole_expansion[2,1,i_compressed] += -strength_J_over_4pi * complex_multiply_imag(inm_real, inm_imag, iam_real, iam_imag) * sign
end

@inline function update_multipole_expansion_panel!(multipole_expansion, i_compressed, strength_J_over_4pi, inm_real, inm_imag, iam_real, iam_imag, inm_prev_mm1_real, inm_prev_mm1_imag, inm_prev_m_real, inm_prev_m_imag, inm_prev_mp1_real, inm_prev_mp1_imag, n, panel_type::UniformNormalDipolePanel{sign}, negative_1_n) where sign
    sum_iterm_real = (inm_prev_mp1_real + inm_prev_mm1_real)/2
    sum_iterm_imag = (inm_prev_mp1_imag + inm_prev_mm1_imag)/2
    x1_real, x1_imag = complex_multiply(inm_prev_mm1_real + inm_prev_mp1_real, inm_prev_mm1_imag + inm_prev_mp1_imag, 0.0, n[1]/2)
    x2_real, x2_imag = n[2] * (inm_prev_mp1_real - inm_prev_mm1_real)/2, n[2] * (inm_prev_mp1_imag - inm_prev_mm1_imag)/2
    x3_real, x3_imag = n[3] * inm_prev_m_real, n[3] * inm_prev_m_imag
    lnm_real, lnm_imag = x1_real + x2_real - x3_real, x1_imag + x2_imag - x3_imag
    # coeff_gumerov = strength_J_over_4pi * (lnm_real - im*lnm_imag) * (-1)^m / (1.0im)^m
    Mnm_real, Mnm_imag = complex_multiply(lnm_real, -lnm_imag, iam_real, iam_imag, negative_1_n*strength_J_over_4pi, 0.0)

    multipole_expansion[1,1,i_compressed] += Mnm_real * sign
    multipole_expansion[2,1,i_compressed] += Mnm_imag * sign
end

@inline function update_multipole_expansion_panel!(multipole_expansion, i_compressed, strength_J_over_4pi, inm_real, inm_imag, iam_real, iam_imag, inm_prev_mm1_real, inm_prev_mm1_imag, inm_prev_m_real, inm_prev_m_imag, inm_prev_mp1_real, inm_prev_mp1_imag, n, panel_type::UniformSourceNormalDipolePanel{sign}, negative_1_n) where {sign}
    update_multipole_expansion_panel!(multipole_expansion, i_compressed, strength_J_over_4pi[1], inm_real, inm_imag, iam_real, iam_imag, inm_prev_mm1_real, inm_prev_mm1_imag, inm_prev_m_real, inm_prev_m_imag, inm_prev_mp1_real, inm_prev_mp1_imag, n, UniformSourcePanel(sign), negative_1_n)
    update_multipole_expansion_panel!(multipole_expansion, i_compressed, strength_J_over_4pi[2], inm_real, inm_imag, iam_real, iam_imag, inm_prev_mm1_real, inm_prev_mm1_imag, inm_prev_m_real, inm_prev_m_imag, inm_prev_mp1_real, inm_prev_mp1_imag, n, UniformNormalDipolePanel(sign), negative_1_n)
end

@inline function update_multipole_expansion_panel!(multipole_expansion, i_compressed, strength_J_over_4pi, inm_real, inm_imag, panel_type::UniformSourcePanel{sign}) where sign
    multipole_expansion[1,1,i_compressed] += strength_J_over_4pi * inm_real * sign
    multipole_expansion[2,1,i_compressed] += -strength_J_over_4pi * inm_imag * sign
end

@inline function update_multipole_expansion_panel!(multipole_expansion, i_compressed, strength_J_over_4pi, inm_real, inm_imag, panel_type::UniformNormalDipolePanel)
    return nothing
end

@inline function update_multipole_expansion_panel!(multipole_expansion, i_compressed, strength_J_over_4pi, inm_real, inm_imag, panel_type::UniformSourceNormalDipolePanel{sign}) where sign
    update_multipole_expansion_panel!(multipole_expansion, i_compressed, strength_J_over_4pi[1], inm_real, inm_imag, UniformSourcePanel(sign))
    update_multipole_expansion_panel!(multipole_expansion, i_compressed, strength_J_over_4pi[2], inm_real, inm_imag, UniformNormalDipolePanel(sign))
end

@inline function _B2M!_panel(multipole_expansion, qnm_prev, jnm_prev, inm_prev, R0, Ru, Rv, strength, normal, expansion_order::Val{P}, panel_type::AbstractKernel) where P
    J = norm(cross(Ru,Rv))

    # transform into xi, eta, z
    xi0_real = R0[1]/2
    xi0_imag = R0[2]/2
    xiu_real = Ru[1]/2
    xiu_imag = Ru[2]/2
    xiv_real = Rv[1]/2
    xiv_imag = Rv[2]/2
    eta0_real = R0[1]/2
    eta0_imag = -R0[2]/2
    etau_real = Ru[1]/2
    etau_imag = -Ru[2]/2
    etav_real = Rv[1]/2
    etav_imag = -Rv[2]/2
    z0 = R0[3]
    zu = Ru[3]
    zv = Rv[3]

    # initialize recurrent values
    qnm_real = 1.0
    qnm_imag = 0.0
    qnm_m1_real = 0.0
    qnm_m1_imag = 0.0
    jnm_real = 1.0
    jnm_imag = 0.0
    jnm_m1_real = 0.0
    jnm_m1_imag = 0.0
    inm_real = 0.5
    inm_imag = 0.0
    inm_m1_real = 0.0
    inm_m1_imag = 0.0
    strength_J_over_4pi = J/4/pi * strength

    # expansion coefficient for l=0, m=0
    dim = 1
    # multipole_expansion[1,dim,1] += real(inm * strength_J_over_4pi)
    # multipole_expansion[2,dim,1] += imag(inm * strength_J_over_4pi) # always zero
    # multipole_expansion[1,dim,1] += inm_real * strength_J_over_4pi
    # multipole_expansion[2,dim,1] += inm_imag * strength_J_over_4pi # always zero

    update_multipole_expansion_panel!(multipole_expansion, 1, strength_J_over_4pi, inm_real, inm_imag, panel_type)

    qnm_prev[1,1] = qnm_real
    qnm_prev[2,1] = qnm_imag
    jnm_prev[1,1] = jnm_real
    jnm_prev[2,1] = jnm_imag
    inm_prev[1,1] = inm_real
    inm_prev[2,1] = inm_imag

    # i^(-1)
    one_over_i_real = 0
    one_over_i_imag = -1

    # recurse
    for n in 1:P
        # m=0
        iam_real = 1
        iam_imag = 0

        # qnm = (im*(xi0+xiu)*-conj(qnm_prev[2]) + im*(eta0+etau)*qnm_prev[2] - (z0+zu)*qnm_prev[1])/n
        x1_real, x1_imag = complex_multiply(xi0_real+xiu_real, xi0_imag+xiu_imag, -qnm_prev[1,2], qnm_prev[2,2], 0, 1)
        x2_real, x2_imag = complex_multiply(eta0_real+etau_real, eta0_imag+etau_imag, qnm_prev[1,2], qnm_prev[2,2], 0, 1)
        x3_real = (z0+zu)*qnm_prev[1,1]
        x3_imag = (z0+zu)*qnm_prev[2,1]
        qnm_real = (x1_real + x2_real - x3_real) / n
        qnm_imag = (x1_imag + x2_imag - x3_imag) / n

        # jnm = (im*(xi0+xiv)*-conj(jnm_prev[2]) + im*(eta0+etav)*jnm_prev[2] - (z0+zv)*jnm_prev[1] + qnm)/(n+1)
        x1_real, x1_imag = complex_multiply(xi0_real+xiv_real, xi0_imag+xiv_imag, -jnm_prev[1,2], jnm_prev[2,2], 0, 1)
        x2_real, x2_imag = complex_multiply(eta0_real+etav_real, eta0_imag+etav_imag, jnm_prev[1,2], jnm_prev[2,2], 0, 1)
        x3_real = (z0+zv) * jnm_prev[1,1]
        x3_imag = (z0+zv) * jnm_prev[2,1]
        jnm_real = (x1_real + x2_real - x3_real + qnm_real) / (n+1)
        jnm_imag = (x1_imag + x2_imag - x3_imag + qnm_imag) / (n+1)

        # inm = (im*xi0*-conj(inm_prev[2]) + im*eta0*inm_prev[2] - z0*inm_prev[1] + jnm)/(n+2)
        x1_real, x1_imag = complex_multiply(xi0_real, xi0_imag, -inm_prev[1,2], inm_prev[2,2], 0, 1)
        x2_real, x2_imag = complex_multiply(eta0_real, eta0_imag, inm_prev[1,2], inm_prev[2,2], 0, 1)
        x3_real, x3_imag = z0 * inm_prev[1,1], z0 * inm_prev[2,1]
        inm_real = (x1_real + x2_real - x3_real + jnm_real) / (n+2)
        inm_imag = (x1_imag + x2_imag - x3_imag + jnm_imag) / (n+2)

        i_compressed = 1 + (n * (n + 1)) >> 1
        # multipole_expansion[1,1,i_compressed] += real(conj(inm * strength_J_over_4pi))
        # multipole_expansion[2,1,i_compressed] += imag(conj(inm * strength_J_over_4pi))
        # multipole_expansion[1,1,i_compressed] += inm_real * strength_J_over_4pi
        # multipole_expansion[2,1,i_compressed] += -inm_imag * strength_J_over_4pi

        # update_multipole_expansion_panel!(multipole_expansion, i_compressed, strength_J_over_4pi, inm_real, inm_imag, panel_type)
        update_multipole_expansion_panel!(multipole_expansion, i_compressed, strength_J_over_4pi, inm_real, inm_imag, iam_real, iam_imag, -inm_prev[1,2], inm_prev[2,2], inm_prev[1,1], inm_prev[2,1], inm_prev[1,2], inm_prev[2,2], normal, panel_type, 1)

        qnm_m1_real = qnm_real
        qnm_m1_imag = qnm_imag
        jnm_m1_real = jnm_real
        jnm_m1_imag = jnm_imag
        inm_m1_real = inm_real
        inm_m1_imag = inm_imag

        # (-1)^m
        negative_1_m = -1

        for m in 1:n
            iam_real, iam_imag = complex_multiply(iam_real,iam_imag,one_over_i_real,one_over_i_imag)

            # qnm = (im*(xi0+xiu)*qnm_prev[m] + im*(eta0+etau)*qnm_prev[m+2] - (z0+zu)*qnm_prev[m+1])/n
            x1_real, x1_imag = complex_multiply(xi0_real+xiu_real, xi0_imag+xiu_imag, qnm_prev[1,m], qnm_prev[2,m], 0, 1)
            x2_real, x2_imag = complex_multiply(eta0_real+etau_real, eta0_imag+etau_imag, qnm_prev[1,m+2], qnm_prev[2,m+2], 0, 1)
            x3_real, x3_imag = (z0+zu)*qnm_prev[1,m+1], (z0+zu)*qnm_prev[2,m+1]
            qnm_real = (x1_real + x2_real - x3_real) / n
            qnm_imag = (x1_imag + x2_imag - x3_imag) / n

            # jnm = (im*(xi0+xiv)*jnm_prev[m] + im*(eta0+etav)*jnm_prev[m+2] - (z0+zv)*jnm_prev[m+1] + qnm)/(n+1)
            x1_real, x1_imag = complex_multiply(xi0_real+xiv_real, xi0_imag+xiv_imag, jnm_prev[1,m], jnm_prev[2,m], 0, 1)
            x2_real, x2_imag = complex_multiply(eta0_real+etav_real, eta0_imag+etav_imag, jnm_prev[1,m+2], jnm_prev[2,m+2], 0, 1)
            x3_real, x3_imag = (z0+zv)*jnm_prev[1,m+1], (z0+zv)*jnm_prev[2,m+1]
            jnm_real = (x1_real + x2_real - x3_real + qnm_real) / (n+1)
            jnm_imag = (x1_imag + x2_imag - x3_imag + qnm_imag) / (n+1)

            # inm = (im*xi0*inm_prev[m] + im*eta0*inm_prev[m+2] - z0*inm_prev[m+1] + jnm)/(n+2)
            x1_real, x1_imag = complex_multiply(xi0_real, xi0_imag, inm_prev[1,m], inm_prev[2,m], 0, 1)
            x2_real, x2_imag = complex_multiply(eta0_real, eta0_imag, inm_prev[1,m+2], inm_prev[2,m+2], 0, 1)
            x3_real, x3_imag = z0*inm_prev[1,m+1], z0*inm_prev[2,m+1]
            inm_real = (x1_real + x2_real - x3_real + jnm_real) / (n+2)
            inm_imag = (x1_imag + x2_imag - x3_imag + jnm_imag) / (n+2)

            i_compressed += 1
            # multipole_expansion[1,1,i_compressed] += real(conj(inm * strength_J_over_4pi * iam))
            # multipole_expansion[2,1,i_compressed] += imag(conj(inm * strength_J_over_4pi * iam))

            # multipole_expansion[1,1,i_compressed] += strength_J_over_4pi * complex_multiply_real(inm_real, inm_imag, iam_real, iam_imag)
            # multipole_expansion[2,1,i_compressed] += -strength_J_over_4pi * complex_multiply_imag(inm_real, inm_imag, iam_real, iam_imag)
            update_multipole_expansion_panel!(multipole_expansion, i_compressed, strength_J_over_4pi, inm_real, inm_imag, iam_real, iam_imag, inm_prev[1,m], inm_prev[2,m], inm_prev[1,m+1], inm_prev[2,m+1], inm_prev[1,m+2], inm_prev[2,m+2], normal, panel_type, negative_1_m)

            qnm_prev[1,m] = qnm_m1_real
            qnm_prev[2,m] = qnm_m1_imag
            jnm_prev[1,m] = jnm_m1_real
            jnm_prev[2,m] = jnm_m1_imag
            inm_prev[1,m] = inm_m1_real
            inm_prev[2,m] = inm_m1_imag
            qnm_m1_real = qnm_real
            qnm_m1_imag = qnm_imag
            jnm_m1_real = jnm_real
            jnm_m1_imag = jnm_imag
            inm_m1_real = inm_real
            inm_m1_imag = inm_imag

            # update (-1)^m
            negative_1_m *= -1
        end
        qnm_prev[1,n+1] = qnm_real
        qnm_prev[2,n+1] = qnm_imag
        jnm_prev[1,n+1] = jnm_real
        jnm_prev[2,n+1] = jnm_imag
        inm_prev[1,n+1] = inm_real
        inm_prev[2,n+1] = inm_imag
    end
end

function B2M!_quadpanel(system, branch, bodies_index, harmonics::AbstractArray{TF}, expansion_order::Val{P}, panel_type::AbstractKernel{sign}; compute_normal=false, invert_normal=false) where {TF,P,sign}
    if P > 2
        harmonics .= zero(TF)
        qnm_prev = view(harmonics,1:2,1:P+2)
        jnm_prev = view(harmonics,1:2,P+3:P<<1+4)
        inm_prev = view(harmonics,1:2,P<<1+5:3*P+6)
    else
        qnm_prev = zeros(TF, 2, P+2)
        jnm_prev = zeros(TF, 2, P+2)
        inm_prev = zeros(TF, 2, P+2)
    end
    for i_body in bodies_index
        strength = system[i_body, STRENGTH]
        R0_global = system[i_body, VERTEX, 1]
        R0 = R0_global - branch.center
        for i_side in 2:3
            Ru = system[i_body,VERTEX,i_side] - R0_global
            Rv = system[i_body,VERTEX,i_side + 1] - R0_global
            if compute_normal
                normal = cross(Rv, Ru)
                normal /= norm(normal) * (-1)^invert_normal
            else
                normal = system[i_body,NORMAL]
            end
            _B2M!_panel(branch.multipole_expansion, qnm_prev, jnm_prev, inm_prev, R0, Ru, Rv, strength, normal, expansion_order, panel_type)
        end
    end
end

function B2M!_tripanel(system, branch, bodies_index, harmonics::AbstractArray{TF}, expansion_order::Val{P}, panel_type::AbstractKernel) where {TF,P}
    if P > 2
        harmonics .= zero(TF)
        qnm_prev = view(harmonics,1:2,1:P+2)
        jnm_prev = view(harmonics,1:2,P+3:P<<1+4)
        inm_prev = view(harmonics,1:2,P<<1+5:3*P+6)
    else
        qnm_prev = zeros(TF, 2, P+2)
        jnm_prev = zeros(TF, 2, P+2)
        inm_prev = zeros(TF, 2, P+2)
    end
    for i_body in bodies_index
        strength = system[i_body, STRENGTH]
        R0_global = system[i_body, VERTEX, 1]
        R0 = R0_global - branch.center
        Ru = system[i_body,VERTEX,2] - R0_global
        Rv = system[i_body,VERTEX,3] - R0_global
        normal = system[i_body,NORMAL]
        _B2M!_panel(branch.multipole_expansion, qnm_prev, jnm_prev, inm_prev, R0, Ru, Rv, strength, normal, expansion_order, panel_type)
    end
end
