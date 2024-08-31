
#------- multipole generation convenience functions -------#

#--- point elements ---#

@inline function B2M!(element::Type{<:Point}, system, branch, bodies_index, harmonics, expansion_order::Val{P}) where P

    c_x, c_y, c_z = branch.center
    multipole_coefficients = branch.multipole_expansion
    for i_body in bodies_index
        b_x, b_y, b_z = system[i_body,POSITION]
        x = b_x - c_x
        y = b_y - c_y
        z = b_z - c_z
        r, theta, phi = cartesian_2_spherical(x,y,z)
        regular_harmonic!(harmonics, r, theta, -phi, P) # Ylm^* -> -dx[3]

        # update values
        B2M_point!(element, multipole_coefficients, harmonics, system, i_body, expansion_order)
    end

end

@inline function B2M_point!(element, multipole_coefficients, harmonics, system, i_body, expansion_order) where scale
    multipole_coefficients = view(multipole_coefficients,:,1,:)
    q = system[i_body, STRENGTH]
    B2M_point!(element, multipole_coefficients, harmonics, q, expansion_order)
end

@inline function B2M_point!(::Type{Point{Source{scale}}}, multipole_coefficients, harmonics, q, expansion_order::Val{P}) where {scale,P}
    multipole_coefficients = view(multipole_coefficients,:,1,:)
    q *= scale
    for l in 0:P
        for m in 0:l
            i_solid_harmonic = l*l + l + m + 1
            i_compressed = 1 + (l * (l + 1)) >> 1 + m # only save half as Yl{-m} = conj(Ylm)
            multipole_coefficients[1,i_compressed] += harmonics[1,i_solid_harmonic] * q
            multipole_coefficients[2,i_compressed] += harmonics[2,i_solid_harmonic] * q
        end
    end
end

#======= TO-DO =======#
#=
#--- filament elements ---#

@inline function B2M!(element::Filament, system, branch, bodies_index, harmonics, expansion_order::Val{P}) where P
    println("Filament")
    c_x, c_y, c_z = branch.center
    multipole_coefficients = branch.multipole_expansion
    for i_body in bodies_index
        b_x, b_y, b_z = system[i_body,POSITION]
        x = b_x - c_x
        z = b_z - c_z
        y = b_y - c_y
        r, theta, phi = cartesian_2_spherical(x,y,z)
        regular_harmonic!(harmonics, r ,theta, -phi, P) # Ylm^* -> -dx[3]

        # update values
        B2M_point!(multipole_coefficients, harmonics, system, i_body, element)
    end

end

#--- panel elements ---#

function B2M!(element::Type{Panel{NS,<:Any}}, system, branch, bodies_index, harmonics, expansion_order::Val{P}) where {NS,P}

    # identify containers
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

    # loop over panels
    for i_body in bodies_index
        strength = system[i_body, SCALAR_STRENGTH]
        R0 = system[i_body, VERTEX, 1] - branch.center
        normal = system[i_body,NORMAL]
        for i_side in 2:NS-1
            Ru = system[i_body,VERTEX,i_side] - R0_global
            Rv = system[i_body,VERTEX,i_side + 1] - R0_global
            if compute_normal
                normal = cross(Rv, Ru)
                normal /= norm(normal) * (-1)^invert_normal
            end
            B2M_panel!(element, branch.multipole_expansion, qnm_prev, jnm_prev, inm_prev, R0, Ru, Rv, normal, system, i_body, expansion_order)
        end
    end

end

@inline function B2M_panel!(element::Type{Panel{<:Any,<:Union{ConstantSource{scale}, ConstantNormalDipole{scale}, ConstantSourceNormalDipole{scale}}}}, multipole_expansion, qnm_prev, jnm_prev, inm_prev, R0, Ru, Rv, normal, system, i_body, expansion_order) where scale
    strength = system[i_body, SCALAR_STRENGTH] * scale
    coefficients = view(multipole_expansion, :, 1, :)
    _B2M!_panel(coefficients, qnm_prev, jnm_prev, inm_prev, R0, Ru, Rv, strength, normal, expansion_order, element)
end

@inline function B2M_panel!(element::Type{Panel{<:Any,Vortex{scale}}}, multipole_expansion, qnm_prev, jnm_prev, inm_prev, R0, Ru, Rv, normal, system, i_body, expansion_order) where scale
    qx, qy, qz = system[i_body, VECTOR_STRENGTH]

    # x component
    coefficients = view(multipole_expansion, :, 2, :)
    _B2M!_panel(coefficients, qnm_prev, jnm_prev, inm_prev, R0, Ru, Rv, qx * scale, normal, expansion_order, Panel{0,Source{scale}})

    # y component
    coefficients = view(multipole_expansion, :, 3, :)
    _B2M!_panel(coefficients, qnm_prev, jnm_prev, inm_prev, R0, Ru, Rv, qy * scale, normal, expansion_order, Panel{0,Source{scale}})

    # z component
    coefficients = view(multipole_expansion, :, 4, :)
    _B2M!_panel(coefficients, qnm_prev, jnm_prev, inm_prev, R0, Ru, Rv, qz * scale, normal, expansion_order, Panel{0,Source{scale}})
end

@inline function update_multipole_expansion_panel!(coefficients, i_compressed, strength_J, inm_real, inm_imag, iam_real, iam_imag, inm_prev_mm1_real, inm_prev_mm1_imag, inm_prev_m_real, inm_prev_m_imag, inm_prev_mp1_real, inm_prev_mp1_imag, n, panel_type::Type{<:Panel{<:Any,<:ConstantSource}}, negative_1_n)
    coefficients[1,i_compressed] += strength_J * complex_multiply_real(inm_real, inm_imag, iam_real, iam_imag)
    coefficients[2,i_compressed] -= strength_J * complex_multiply_imag(inm_real, inm_imag, iam_real, iam_imag)
end

@inline function update_multipole_expansion_panel!(coefficients, i_compressed, strength_J, inm_real, inm_imag, iam_real, iam_imag, inm_prev_mm1_real, inm_prev_mm1_imag, inm_prev_m_real, inm_prev_m_imag, inm_prev_mp1_real, inm_prev_mp1_imag, n, panel_type::Type{<:Panel{<:Any,<:ConstantNormalDipole}}, negative_1_n)
    sum_iterm_real = (inm_prev_mp1_real + inm_prev_mm1_real)/2
    sum_iterm_imag = (inm_prev_mp1_imag + inm_prev_mm1_imag)/2
    x1_real, x1_imag = complex_multiply(inm_prev_mm1_real + inm_prev_mp1_real, inm_prev_mm1_imag + inm_prev_mp1_imag, 0.0, n[1]/2)
    x2_real, x2_imag = n[2] * (inm_prev_mp1_real - inm_prev_mm1_real)/2, n[2] * (inm_prev_mp1_imag - inm_prev_mm1_imag)/2
    x3_real, x3_imag = n[3] * inm_prev_m_real, n[3] * inm_prev_m_imag
    lnm_real, lnm_imag = x1_real + x2_real - x3_real, x1_imag + x2_imag - x3_imag
    # coeff_gumerov = strength_J * (lnm_real - im*lnm_imag) * (-1)^m / (1.0im)^m
    Mnm_real, Mnm_imag = complex_multiply(lnm_real, -lnm_imag, iam_real, iam_imag, negative_1_n*strength_J, 0.0)

    coefficients[1,i_compressed] += Mnm_real
    coefficients[2,i_compressed] += Mnm_imag
end

@inline function update_multipole_expansion_panel!(coefficients, i_compressed, strength_J, inm_real, inm_imag, iam_real, iam_imag, inm_prev_mm1_real, inm_prev_mm1_imag, inm_prev_m_real, inm_prev_m_imag, inm_prev_mp1_real, inm_prev_mp1_imag, n, panel_type::Type{<:Panel{<:Any,<:ConstantSourceNormalDipole}}, negative_1_n)
    update_multipole_expansion_panel!(coefficients, i_compressed, strength_J[1], inm_real, inm_imag, iam_real, iam_imag, inm_prev_mm1_real, inm_prev_mm1_imag, inm_prev_m_real, inm_prev_m_imag, inm_prev_mp1_real, inm_prev_mp1_imag, n, Panel{<:Any,ConstantSource{scale}}, negative_1_n)
    update_multipole_expansion_panel!(coefficients, i_compressed, strength_J[2], inm_real, inm_imag, iam_real, iam_imag, inm_prev_mm1_real, inm_prev_mm1_imag, inm_prev_m_real, inm_prev_m_imag, inm_prev_mp1_real, inm_prev_mp1_imag, n, Panel{<:Any,ConstantNormalDipole{scale}}, negative_1_n)
end

@inline function update_multipole_expansion_panel!(coefficients, i_compressed, strength_J, inm_real, inm_imag, panel_type::Type{<:Panel{<:Any,<:ConstantSource}})
    coefficients[1,i_compressed] += strength_J * inm_real
    coefficients[2,i_compressed] -= strength_J * inm_imag
end

@inline function update_multipole_expansion_panel!(coefficients, i_compressed, strength_J, inm_real, inm_imag, panel_type::Type{<:Panel{<:Any,<:ConstantNormalDipole}})
    return nothing
end

@inline function update_multipole_expansion_panel!(coefficients, i_compressed, strength_J, inm_real, inm_imag, panel_type::Type{Panel{<:Any,ConstantSourceNormalDipole{scale}}}) where scale
    update_multipole_expansion_panel!(coefficients, i_compressed, strength_J[1], inm_real, inm_imag, Panel{0,ConstantSource{scale}})
    update_multipole_expansion_panel!(coefficients, i_compressed, strength_J[2], inm_real, inm_imag, Panel{0,ConstantNormalDipole{scale}})
end

@inline function _B2M!_panel(multipole_expansion, qnm_prev, jnm_prev, inm_prev, R0, Ru, Rv, strength, normal, expansion_order::Val{P}, panel_type::Type{<:Panel}) where P
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
    strength_J = J * strength

    # expansion coefficient for l=0, m=0
    # coefficients[1,1] += real(inm * strength_J)
    # coefficients[2,1] += imag(inm * strength_J) # always zero
    # coefficients[1,1] += inm_real * strength_J
    # coefficients[2,1] += inm_imag * strength_J # always zero

    update_multipole_expansion_panel!(multipole_expansion, 1, strength_J, inm_real, inm_imag, panel_type)

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
        # coefficients[1,1,i_compressed] += real(conj(inm * strength_J))
        # coefficients[2,1,i_compressed] += imag(conj(inm * strength_J))
        # coefficients[1,1,i_compressed] += inm_real * strength_J
        # coefficients[2,1,i_compressed] += -inm_imag * strength_J

        # update_multipole_expansion_panel!(multipole_expansion, i_compressed, strength_J, inm_real, inm_imag, panel_type)
        update_multipole_expansion_panel!(multipole_expansion, i_compressed, strength_J, inm_real, inm_imag, iam_real, iam_imag, -inm_prev[1,2], inm_prev[2,2], inm_prev[1,1], inm_prev[2,1], inm_prev[1,2], inm_prev[2,2], normal, panel_type, 1)

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
            # coefficients[1,i_compressed] += real(conj(inm * strength_J * iam))
            # coefficients[2,i_compressed] += imag(conj(inm * strength_J * iam))

            # coefficients[1,i_compressed] += strength_J * complex_multiply_real(inm_real, inm_imag, iam_real, iam_imag)
            # coefficients[2,i_compressed] -= strength_J * complex_multiply_imag(inm_real, inm_imag, iam_real, iam_imag)
            update_multipole_expansion_panel!(multipole_expansion, i_compressed, strength_J, inm_real, inm_imag, iam_real, iam_imag, inm_prev[1,m], inm_prev[2,m], inm_prev[1,m+1], inm_prev[2,m+1], inm_prev[1,m+2], inm_prev[2,m+2], normal, panel_type, negative_1_m)

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

# function B2M!_quadpanel(system, branch, bodies_index, harmonics::AbstractArray{TF}, expansion_order::Val{P}, panel_type::AbstractKernel{scale}; compute_normal=false, invert_normal=false) where {TF,P,scale}
#     if P > 2
#         harmonics .= zero(TF)
#         qnm_prev = view(harmonics,1:2,1:P+2)
#         jnm_prev = view(harmonics,1:2,P+3:P<<1+4)
#         inm_prev = view(harmonics,1:2,P<<1+5:3*P+6)
#     else
#         qnm_prev = zeros(TF, 2, P+2)
#         jnm_prev = zeros(TF, 2, P+2)
#         inm_prev = zeros(TF, 2, P+2)
#     end
#     for i_body in bodies_index
#         strength = system[i_body, SCALAR_STRENGTH]
#         R0_global = system[i_body, VERTEX, 1]
#         R0 = R0_global - branch.center
#         for i_side in 2:3
#             Ru = system[i_body,VERTEX,i_side] - R0_global
#             Rv = system[i_body,VERTEX,i_side + 1] - R0_global
#             if compute_normal
#                 normal = cross(Rv, Ru)
#                 normal /= norm(normal) * (-1)^invert_normal
#             else
#                 normal = system[i_body,NORMAL]
#             end
#             _B2M!_panel(branch.multipole_expansion, qnm_prev, jnm_prev, inm_prev, R0, Ru, Rv, strength, normal, expansion_order, panel_type)
#         end
#     end
# end
#
# function B2M!_tripanel(system, branch, bodies_index, harmonics::AbstractArray{TF}, expansion_order::Val{P}, panel_type::AbstractKernel) where {TF,P}
#     if P > 2
#         harmonics .= zero(TF)
#         qnm_prev = view(harmonics,1:2,1:P+2)
#         jnm_prev = view(harmonics,1:2,P+3:P<<1+4)
#         inm_prev = view(harmonics,1:2,P<<1+5:3*P+6)
#     else
#         qnm_prev = zeros(TF, 2, P+2)
#         jnm_prev = zeros(TF, 2, P+2)
#         inm_prev = zeros(TF, 2, P+2)
#     end
#     for i_body in bodies_index
#         strength = system[i_body, SCALAR_STRENGTH]
#         R0_global = system[i_body, VERTEX, 1]
#         R0 = R0_global - branch.center
#         Ru = system[i_body,VERTEX,2] - R0_global
#         Rv = system[i_body,VERTEX,3] - R0_global
#         normal = system[i_body,NORMAL]
#         _B2M!_panel(branch.multipole_expansion, qnm_prev, jnm_prev, inm_prev, R0, Ru, Rv, strength, normal, expansion_order, panel_type)
#     end
# end
=#
