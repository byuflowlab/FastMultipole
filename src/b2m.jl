
#####
##### multipole generation convenience functions
#####
@inline function B2M!_vortexpoint(system, branch, bodies_index, harmonics, expansion_order)
    c_x, c_y, c_z = branch.center
    for i_body in bodies_index
        b_x, b_y, b_z = system[i_body,POSITION]
        x = b_x - c_x
        z = b_z - c_z
        y = b_y - c_y
        qx, qy, qz = system[i_body,VECTOR_STRENGTH]
        r, theta, phi = cartesian_2_spherical(x,y,z)
        regular_harmonic!(harmonics, r ,theta, -phi, expansion_order) # Ylm^* -> -dx[3]
        # update values
        for l in 0:expansion_order
            for m in 0:l
                i_solid_harmonic = l*l + l + m + 1
                i_compressed = 1 + (l * (l + 1)) >> 1 + m # only save half as Yl{-m} = conj(Ylm)
                branch.multipole_expansion[2,i_compressed] += harmonics[i_solid_harmonic] * qx * ONE_OVER_4PI
                branch.multipole_expansion[3,i_compressed] += harmonics[i_solid_harmonic] * qy * ONE_OVER_4PI
                branch.multipole_expansion[4,i_compressed] += harmonics[i_solid_harmonic] * qz * ONE_OVER_4PI
            end
        end
    end
end

@inline function B2M!_sourcepoint(system, branch, bodies_index, harmonics, expansion_order)
    c_x, c_y, c_z = branch.center
    for i_body in bodies_index
        b_x, b_y, b_z = system[i_body,POSITION]
        x = b_x - c_x
        z = b_z - c_z
        y = b_y - c_y
        q = system[i_body,SCALAR_STRENGTH]
        r, theta, phi = cartesian_2_spherical(x,y,z)
        regular_harmonic!(harmonics, r, theta, -phi, expansion_order) # Ylm^* -> -dx[3]
        # update values
        for l in 0:expansion_order
            for m in 0:l
                i_solid_harmonic = l*l + l + m + 1
                i_compressed = 1 + (l * (l + 1)) >> 1 + m # only save half as Yl{-m} = conj(Ylm)
                branch.multipole_expansion[1,i_compressed] += harmonics[i_solid_harmonic] * q
            end
        end
    end
end

# function B2M!_sourcepanel(system, branch, bodies_index, harmonics, expansion_order)
#     qnm_prev = zeros(Complex{Float64}, expansion_order+2)
#     jnm_prev = zeros(Complex{Float64}, expansion_order+2)
#     inm_prev = zeros(Complex{Float64}, expansion_order+2)
#     for i_body in bodies_index
#         strength = system[i_body, SCALAR_STRENGTH]
#         R0_global = system[i_body, POSITION]
#         R0 = R0_global - branch.center
#         n_sides = length()
#         for i_side in 1:4
#             i_side_p1 = i_side == 4 ? 1 : i_side + 1
#             Ru = system[i_body,VERTEX,i_side] - R0_global
#             Rv = system[i_body,VERTEX,i_side_p1] - R0_global
#             J = norm(cross(Ru,Rv))

#             # transform into xi, eta, z
#             xi0 = (R0[1] + im*R0[2])/2
#             xiu = (Ru[1] + im*Ru[2])/2
#             xiv = (Rv[1] + im*Rv[2])/2
#             eta0 = (R0[1] - im*R0[2])/2
#             etau = (Ru[1] - im*Ru[2])/2
#             etav = (Rv[1] - im*Rv[2])/2
#             z0 = R0[3]
#             zu = Ru[3]
#             zv = Rv[3]

#             # initialize recurrent values
#             qnm = 1.0
#             qnm_m1 = 0.0
#             jnm = 1.0
#             jnm_m1 = 0.0
#             inm = 0.5
#             inm_m1 = 0.0
#             strength_J_over_4pi = J/4/pi * strength

#             # expansion coefficient for l=0, m=0
#             dim = 1
#             branch.multipole_expansion[dim,1] += inm * strength_J_over_4pi
#             qnm_prev[1] = qnm
#             jnm_prev[1] = jnm
#             inm_prev[1] = inm

#             # i^(-1)
#             one_over_i = -im

#             # recurse
#             for n in 1:expansion_order
#                 # m=0
#                 iam = Complex(1)
#                 qnm = (im*(xi0+xiu)*-conj(qnm_prev[2]) + im*(eta0+etau)*qnm_prev[2] - (z0+zu)*qnm_prev[1])/n
#                 jnm = (im*(xi0+xiv)*-conj(jnm_prev[2]) + im*(eta0+etav)*jnm_prev[2] - (z0+zv)*jnm_prev[1] + qnm)/(n+1)
#                 inm = (im*xi0*-conj(inm_prev[2]) + im*eta0*inm_prev[2] - z0*inm_prev[1] + jnm)/(n+2)
#                 i_compressed = 1 + (n * (n + 1)) >> 1
#                 branch.multipole_expansion[1,i_compressed] += conj(inm * strength_J_over_4pi)
#                 qnm_m1 = qnm
#                 jnm_m1 = jnm
#                 inm_m1 = inm
#                 for m in 1:n
#                     iam *= one_over_i
#                     qnm = (im*(xi0+xiu)*qnm_prev[m] + im*(eta0+etau)*qnm_prev[m+2] - (z0+zu)*qnm_prev[m+1])/n
#                     jnm = (im*(xi0+xiv)*jnm_prev[m] + im*(eta0+etav)*jnm_prev[m+2] - (z0+zv)*jnm_prev[m+1] + qnm)/(n+1)
#                     inm = (im*xi0*inm_prev[m] + im*eta0*inm_prev[m+2] - z0*inm_prev[m+1] + jnm)/(n+2)
#                     i_compressed += 1
#                     branch.multipole_expansion[1,i_compressed] += conj(inm * strength_J_over_4pi * iam)
#                     qnm_prev[m] = qnm_m1
#                     jnm_prev[m] = jnm_m1
#                     inm_prev[m] = inm_m1
#                     qnm_m1 = qnm
#                     jnm_m1 = jnm
#                     inm_m1 = inm
#                 end
#                 qnm_prev[n+1] = qnm
#                 jnm_prev[n+1] = jnm
#                 inm_prev[n+1] = inm
#             end
#         end
#     end
# end

@inline function _B2M!_panel(branch, qnm_prev, jnm_prev, inm_prev, R0, Ru, Rv, strength, expansion_order::Val{P}) where P
    J = norm(cross(Ru,Rv))

    # transform into xi, eta, z
    xi0 = (R0[1] + im*R0[2])/2
    xiu = (Ru[1] + im*Ru[2])/2
    xiv = (Rv[1] + im*Rv[2])/2
    eta0 = (R0[1] - im*R0[2])/2
    etau = (Ru[1] - im*Ru[2])/2
    etav = (Rv[1] - im*Rv[2])/2
    z0 = R0[3]
    zu = Ru[3]
    zv = Rv[3]

    # initialize recurrent values
    qnm = 1.0
    qnm_m1 = 0.0
    jnm = 1.0
    jnm_m1 = 0.0
    inm = 0.5
    inm_m1 = 0.0
    strength_J_over_4pi = J/4/pi * strength

    # expansion coefficient for l=0, m=0
    dim = 1
    branch.multipole_expansion[dim,1] += inm * strength_J_over_4pi
    qnm_prev[1] = qnm
    jnm_prev[1] = jnm
    inm_prev[1] = inm

    # i^(-1)
    one_over_i = -im

    # recurse
    for n in 1:P
        # m=0
        iam = Complex(1)
        qnm = (im*(xi0+xiu)*-conj(qnm_prev[2]) + im*(eta0+etau)*qnm_prev[2] - (z0+zu)*qnm_prev[1])/n
        jnm = (im*(xi0+xiv)*-conj(jnm_prev[2]) + im*(eta0+etav)*jnm_prev[2] - (z0+zv)*jnm_prev[1] + qnm)/(n+1)
        inm = (im*xi0*-conj(inm_prev[2]) + im*eta0*inm_prev[2] - z0*inm_prev[1] + jnm)/(n+2)
        i_compressed = 1 + (n * (n + 1)) >> 1
        branch.multipole_expansion[1,i_compressed] += conj(inm * strength_J_over_4pi)
        qnm_m1 = qnm
        jnm_m1 = jnm
        inm_m1 = inm
        for m in 1:n
            iam *= one_over_i
            qnm = (im*(xi0+xiu)*qnm_prev[m] + im*(eta0+etau)*qnm_prev[m+2] - (z0+zu)*qnm_prev[m+1])/n
            jnm = (im*(xi0+xiv)*jnm_prev[m] + im*(eta0+etav)*jnm_prev[m+2] - (z0+zv)*jnm_prev[m+1] + qnm)/(n+1)
            inm = (im*xi0*inm_prev[m] + im*eta0*inm_prev[m+2] - z0*inm_prev[m+1] + jnm)/(n+2)
            i_compressed += 1
            branch.multipole_expansion[1,i_compressed] += conj(inm * strength_J_over_4pi * iam)
            qnm_prev[m] = qnm_m1
            jnm_prev[m] = jnm_m1
            inm_prev[m] = inm_m1
            qnm_m1 = qnm
            jnm_m1 = jnm
            inm_m1 = inm
        end
        qnm_prev[n+1] = qnm
        jnm_prev[n+1] = jnm
        inm_prev[n+1] = inm
    end
end

function B2M!_sourcequadpanel(system, branch, bodies_index, harmonics, expansion_order::Val{P}) where P
    if P > 2
        harmonics .= zero(eltype(harmonics))
        qnm_prev = view(harmonics,1:P+2)
        jnm_prev = view(harmonics,P+3:P<<1+4)
        inm_prev = view(harmonics,P<<1+5:3*P+6)
    else
        qnm_prev = zeros(Complex{Float64}, expansion_order+2)
        jnm_prev = zeros(Complex{Float64}, expansion_order+2)
        inm_prev = zeros(Complex{Float64}, expansion_order+2)
    end
    for i_body in bodies_index
        strength = system[i_body, SCALAR_STRENGTH]
        R0_global = system[i_body, VERTEX, 1]
        R0 = R0_global - branch.center
        for i_side in 2:3
            Ru = system[i_body,VERTEX,i_side] - R0_global
            Rv = system[i_body,VERTEX,i_side + 1] - R0_global
            _B2M!_panel(branch, qnm_prev, jnm_prev, inm_prev, R0, Ru, Rv, strength, expansion_order)
        end
    end
end

function B2M!_sourcetripanel(system, branch, bodies_index, harmonics, expansion_order::Val{P}) where P
    if P > 2
        harmonics .= zero(eltype(harmonics))
        qnm_prev = view(harmonics,1:P+2)
        jnm_prev = view(harmonics,P+3:P<<1+4)
        inm_prev = view(harmonics,P<<1+5:3*P+6)
    else
        qnm_prev = zeros(Complex{Float64}, expansion_order+2)
        jnm_prev = zeros(Complex{Float64}, expansion_order+2)
        inm_prev = zeros(Complex{Float64}, expansion_order+2)
    end
    for i_body in bodies_index
        strength = system[i_body, SCALAR_STRENGTH]
        R0_global = system[i_body, VERTEX, 1]
        R0 = R0_global - branch.center
        Ru = system[i_body,VERTEX,2] - R0_global
        Rv = system[i_body,VERTEX,3] - R0_global
        _B2M!_panel(branch, qnm_prev, jnm_prev, inm_prev, R0, Ru, Rv, strength, expansion_order)
    end
end
