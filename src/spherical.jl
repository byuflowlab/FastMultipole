function cartesian_2_spherical!(dx; EPSILON=1e-10)
    r = sqrt(dx' * dx)
    theta = r < EPSILON ? 0.0 : acos(dx[3] / r)
    phi = r < EPSILON ? 0.0 : atan(dx[2], dx[1])
    dx .= r, theta, phi
end

function cartesian_2_spherical(dx; EPSILON=1e-10)
    r = sqrt(dx' * dx)
    theta = r < EPSILON ? 0.0 : acos(dx[3] / r)
    phi = r < EPSILON ? 0.0 : atan(dx[2], dx[1])
    dx = SVector{3}(r, theta, phi)
end

function cartesian_2_spherical(x, y, z; EPSILON=1e-10)
    r = sqrt(x*x + y*y + z*z)
    theta = r < EPSILON ? 0.0 : acos(z / r)
    phi = r < EPSILON ? 0.0 : atan(y, x)
    return r, theta, phi
end

spherical_2_cartesian!(potential_jacobian, potential_hessian, workspace, rho, theta, phi, derivatives_switch::DerivativesSwitch{PS,false,false}) where PS = nothing

@inline odd_or_even(n::Int) = (n & 1) == 1 ? -1 : 1

@inline ipow2l(n::Int) = n >= 0 ? 1 : odd_or_even(n);

function regular_harmonic!(harmonics, rho, theta, phi::TF, P) where TF
    y,x = sincos(theta)
    fact = 1.0
    pl = 1.0
    rhom = 1.0 # rho^l / (l+m)! * (-1)^l
    #ei = exp(im * phi)
    ei = SVector{2,TF}(cos(phi), sin(phi))
    #eim = 1.0 # e^(i * m * phi)
    eim = SVector{2,TF}(1.0,0.0) # e^(i * m * phi)
    for m=0:P # l=m up here
        p = pl
        lpl = m * m + 2 * m + 1
        lml = m * m + 1
        #harmonics[lpl] = rhom * p * eim
        harmonics[1,lpl] = rhom * p * eim[1]
        harmonics[2,lpl] = rhom * p * eim[2]
        #harmonics[lml] = conj(harmonics[lpl])
        harmonics[1,lml] = harmonics[1,lpl]
        harmonics[2,lml] = -harmonics[2,lpl]
        p1 = p
        p = x * (2 * m + 1) * p1
        rhom *= rho
        rhol = rhom
        for l=m+1:P # l>m in here
            lpm = l * l + l + m + 1
            lmm = l * l + l - m + 1
            rhol /= -(l + m)
            #harmonics[lpm] = rhol * p * eim
            harmonics[1,lpm] = rhol * p * eim[1]
            harmonics[2,lpm] = rhol * p * eim[2]
            #harmonics[lmm] = conj(harmonics[lpm])
            harmonics[1,lmm] = harmonics[1,lpm]
            harmonics[2,lmm] = -harmonics[2,lpm]
            p2 = p1
            p1 = p
            p = (x * (2 * l + 1) * p1 - (l + m) * p2) / (l - m + 1)
            rhol *= rho
        end
        rhom /= -(2 * m + 2) * (2 * m + 1)
        pl = -pl * fact * y
        fact += 2
        #eim *= ei
        eim = SVector{2,TF}(eim[1]*ei[1] - eim[2]*ei[2], eim[1]*ei[2] + eim[2]*ei[1])
    end
end

function irregular_harmonic!(harmonics, rho, theta, phi::TF, P) where TF
    y, x = sincos(theta)
    fact = 1
    pl = 1
    invR = -1.0 / rho
    rhom = -invR
    #ei = exp(im * phi)
    ei = SVector{2,TF}(cos(phi), sin(phi))
    #eim = 1.0 # e^(i * m * phi)
    eim = SVector{2,TF}(1.0,0.0) # e^(i * m * phi)
    for m=0:P
        p = pl
        npl = m * m + 2 * m + 1
        nml = m * m + 1
        #harmonics[npl] = rhom * p * eim
        harmonics[1,npl] = rhom * p * eim[1]
        harmonics[2,npl] = rhom * p * eim[2]
        #harmonics[nml] = conj(harmonics[npl])
        harmonics[1,nml] = harmonics[1,npl]
        harmonics[2,nml] = -harmonics[2,npl]
        p1 = p
        p = x * (2 * m + 1) * p1
        rhom *= invR
        rhon = rhom
        for l=m+1:P
            npm = l * l + l + m + 1
            nmm = l * l + l - m + 1
            # npm_max = P^2 + P + P + 1 = (P+1)^2
            #harmonics[npm] = rhon * p * eim
            harmonics[1,npm] = rhon * p * eim[1]
            harmonics[2,npm] = rhon * p * eim[2]
            #harmonics[nmm] = conj(harmonics[npm])
            harmonics[1,nmm] = harmonics[1,npm]
            harmonics[2,nmm] = -harmonics[2,npm]
            p2 = p1
            p1 = p
            p = (x * (2 * l + 1) * p1 - (l + m) * p2) / (l - m + 1)
            rhon *= invR * (l - m + 1)
        end
        pl = -pl * fact * y
        fact += 2
        #eim *= ei
        eim = SVector{2,TF}(eim[1]*ei[1] - eim[2]*ei[2], eim[1]*ei[2] + eim[2]*ei[1])
    end
end

@inline function B2M!(branch, system, harmonics, expansion_order)
    B2M!(system, branch, branch.bodies_index, harmonics, expansion_order)
end

function B2M!(branch, systems::Tuple, harmonics, expansion_order)
    # iterate over systems
    for (system,bodies_index) in zip(systems, branch.bodies_index)
        B2M!(system, branch, bodies_index, harmonics, expansion_order)
    end
end

function M2B!(target_potential, target, center, irregular_harmonics, multipole_expansion, expansion_order)
    dx = target[1:3] - center
    r, theta, phi = cartesian_2_spherical(dx)
    irregular_harmonic!(irregular_harmonics, r, theta, phi, expansion_order)
    d_potential = zeros(4)
    for l in 0:expansion_order
        for m in 0:l
            ip = l^2 + l + m + 1
            i_compressed = 1 + (l * (l + 1)) >> 1 + m # only save half as Yl{-m} = conj(Ylm)
            factor = m > 0 ? 2.0 : 1.0
            for dim in 1:4
                #d_potential[dim] += factor*real(branch.multipole_expansion[dim,i_compressed] * irregular_harmonics[ip])
                d_potential[dim] += factor*(multipole_expansion[1,dim,i_compressed] * irregular_harmonics[1,ip] - multipole_expansion[2,dim,i_compressed] * irregular_harmonics[2,ip])
            end
        end
    end
    target_potential .+= d_potential
end

function M2B!(target_potential, target, center, irregular_harmonics, multipole_expansion::Matrix{Complex{Float64}}, expansion_order)
    dx = target[1:3] - center
    r, theta, phi = cartesian_2_spherical(dx)
    irregular_harmonic!(irregular_harmonics, r, theta, phi, expansion_order)
    d_potential = zeros(4)
    for l in 0:expansion_order
        for m in 0:l
            ip = l^2 + l + m + 1
            i_compressed = 1 + (l * (l + 1)) >> 1 + m # only save half as Yl{-m} = conj(Ylm)
            factor = m > 0 ? 2.0 : 1.0
            for dim in 1:4
                #d_potential[dim] += factor*real(branch.multipole_expansion[dim,i_compressed] * irregular_harmonics[ip])
                d_potential[dim] += factor*(real(multipole_expansion[dim,i_compressed]) * irregular_harmonics[1,ip] - imag(multipole_expansion[dim,i_compressed]) * irregular_harmonics[2,ip])
            end
        end
    end
    target_potential .+= d_potential
end

function M2B!(target_potential, target, i_branch, tree::Tree{<:Any,P}) where P
    branch = tree.branches[i_branch]
    irregular_harmonics = Matrix{eltype(branch.multipole_expansion[1])}(undef, 2, (P+1)^2)
    M2B!(target_potential, target, branch.center, irregular_harmonics, branch.multipole_expansion, P)
end

function M2M!(branch, child, harmonics, M, expansion_order::Val{P}) where P
    # get distance vector
    dx, dy, dz = branch.center - child.center
    r, theta, phi = cartesian_2_spherical(dx, dy, dz)
    regular_harmonic!(harmonics, r, theta, phi, P)

    for j in 0:P # iterate over new Multipole coefficients B_j^k
        for k in 0:j
            i_jk = ((j * (j+1)) >> 1) + k + 1 # current index
            M .= zero(eltype(M))
            for l in 0:j
                for m in max(-l,-j+k+l):min(k-1,l)
                    jlkms = (((j-l) * (j-l+1)) >> 1) + k - m + 1
                    lm = l * l + l - m + 1
                    ipow = ipow2l(m)
                    oddeven = odd_or_even(l)
                    #C_tmp = harmonics[lm] * ipow * oddeven
                    C_tmp1 = harmonics[1,lm] * ipow * oddeven
                    C_tmp2 = harmonics[2,lm] * ipow * oddeven
                    for dim in 1:4
                        #@inbounds M[dim] += child.multipole_expansion[dim,jlkms] * C_tmp
                        @inbounds M[1,dim] += child.multipole_expansion[1,dim,jlkms] * C_tmp1 - child.multipole_expansion[2,dim,jlkms] * C_tmp2
                        @inbounds M[2,dim] += child.multipole_expansion[1,dim,jlkms] * C_tmp2 + child.multipole_expansion[2,dim,jlkms] * C_tmp1
                    end
                end
                for m in k:min(l,j+k-l)
                    jlkms = (((j-l) * (j-l+1)) >> 1) - k + m + 1
                    lm = l * l + l - m + 1
                    oddeven = odd_or_even(k + l + m)
                    #C_tmp = harmonics[lm] * oddeven
                    C_tmp1 = harmonics[1,lm] * oddeven
                    C_tmp2 = harmonics[2,lm] * oddeven
                    for dim in 1:4
                        #@inbounds M[dim] += conj(child.multipole_expansion[dim,jlkms]) * C_tmp
                        @inbounds M[1,dim] += child.multipole_expansion[1,dim,jlkms] * C_tmp1 + child.multipole_expansion[2,dim,jlkms] * C_tmp2
                        @inbounds M[2,dim] += child.multipole_expansion[1,dim,jlkms] * C_tmp2 - child.multipole_expansion[2,dim,jlkms] * C_tmp1
                    end
                end
            end
            for dim in 1:4
                #@inbounds branch.multipole_expansion[dim,i_jk] += M[dim]
                @inbounds branch.multipole_expansion[1,dim,i_jk] += M[1,dim]
                @inbounds branch.multipole_expansion[2,dim,i_jk] += M[2,dim]
            end
        end
    end
end

function M2L_loop!(local_expansion, L, multipole_expansion, harmonics, expansion_order::Val{P}) where P
    for j in 0:P
        Cnm = odd_or_even(j)
        for k in 0:j
            jks = (j * (j + 1)) >> 1 + k + 1
            L .= zero(eltype(L))
            for n in 0:P
                for m in -n:-1
                    nms = (n * (n+1)) >> 1 - m + 1
                    jnkm = (j + n)^2 + j + n + m - k + 1
                    # jnkm_max = (P + P)^2 + P + P + -1 - 0 + 1 = (2P)^2 + 2P = 2P(2P+1)
                    #Cnm_tmp = Cnm * harmonics[jnkm]
                    Cnm_tmp1 = Cnm * harmonics[1,jnkm]
                    Cnm_tmp2 = Cnm * harmonics[2,jnkm]
                    for dim in 1:4
                        #@inbounds L[dim] += conj(multipole_expansion[dim,nms]) * Cnm_tmp
                        @inbounds L[1,dim] += multipole_expansion[1,dim,nms] * Cnm_tmp1 + multipole_expansion[2,dim,nms] * Cnm_tmp2
                        @inbounds L[2,dim] += multipole_expansion[1,dim,nms] * Cnm_tmp2 - multipole_expansion[2,dim,nms] * Cnm_tmp1
                    end
                end
                for m in 0:n
                    nms = (n * (n+1)) >> 1 + m + 1
                    jnkm = (j + n) * (j + n) + j + n + m - k + 1
                    # jnkm_max = 2P * 2P + 2P + P + P - 0 + 1 = (2P)^2 + 2P + 2P + 1 = 4P^2 + 4P + 1 = (2P + 1)^2
                    #Cnm_tmp = Cnm * odd_or_even((k-m) * (1 >> (k>=m)) + m) * harmonics[jnkm]
                    Cnm_tmp1 = Cnm * odd_or_even((k-m) * (1 >> (k>=m)) + m) * harmonics[1,jnkm]
                    Cnm_tmp2 = Cnm * odd_or_even((k-m) * (1 >> (k>=m)) + m) * harmonics[2,jnkm]
                    for dim in 1:4
                        #@inbounds L[dim] += multipole_expansion[dim,nms] * Cnm_tmp
                        @inbounds L[1,dim] += multipole_expansion[1,dim,nms] * Cnm_tmp1 - multipole_expansion[2,dim,nms] * Cnm_tmp2
                        @inbounds L[2,dim] += multipole_expansion[1,dim,nms] * Cnm_tmp2 + multipole_expansion[2,dim,nms] * Cnm_tmp1
                    end
                end
            end
            view(local_expansion,:,:,jks) .+= L
        end
    end
end

function M2L!(target_branch, source_branch, harmonics, L, expansion_order::Val{P}) where P
    dx, dy, dz = target_branch.center - source_branch.center
    r, theta, phi = cartesian_2_spherical(dx, dy, dz)
    irregular_harmonic!(harmonics, r, theta, phi, P<<1)
    M2L_loop!(target_branch.local_expansion, L, source_branch.multipole_expansion, harmonics, expansion_order)
end

function M2L!(target_branch, source_branch, expansion_order::Val{P}) where P
    dx, dy, dz = target_branch.center - source_branch.center
    r, theta, phi = cartesian_2_spherical(dx, dy, dz)
    irregular_harmonic!(target_branch.harmonics, r, theta, phi, P<<1)
    M2L_loop!(target_branch.local_expansion, target_branch.ML, source_branch.multipole_expansion, target_branch.harmonics, expansion_order)
end

function B2L!(tree::Tree{<:Any,P}, i_branch, source_position, source_strength) where P
    branch = tree.branches[i_branch]
    irregular_harmonics = Matrix{eltype(branch.multipole_expansion)}(undef, 2, (P+1)^2)
    r, theta, phi = cartesian_2_spherical(source_position - branch.center)
    irregular_harmonic!(irregular_harmonics, r, theta, -phi, P)
    for l in 0:P
        for m in 0:l
            i_abb = (l * (l+1)) >> 1 + m + 1
            i_exp = l^2 + l + m + 1
            for dim in 1:4
                #branch.local_expansion[dim,i_abb] = irregular_harmonics[i_exp] * source_strength[dim]
                branch.local_expansion[1,dim,i_abb] = irregular_harmonics[1,i_exp] * source_strength[dim]
                branch.local_expansion[2,dim,i_abb] = irregular_harmonics[2,i_exp] * source_strength[dim]
            end
        end
    end
end

function L2L!(branch, child, regular_harmonics, L, expansion_order::Val{P}) where P
    dx, dy, dz = child.center - branch.center
    r, theta, phi = cartesian_2_spherical(dx, dy, dz)
    regular_harmonic!(regular_harmonics, r, theta, phi, P)
    for j in 0:P
        for k in 0:j
            jks = (j * (j + 1)) >> 1 + k + 1
            L .= zero(eltype(branch.local_expansion[1]))
            for n in j:P
                for m in j+k-n:-1
                    jnkm = (n-j) * (n-j) + n - j + m - k + 1
                    nms = (n * (n + 1)) >> 1 - m + 1
                    oddeven = odd_or_even(k)
                    C_tmp1 = regular_harmonics[1,jnkm] * oddeven
                    C_tmp2 = regular_harmonics[2,jnkm] * oddeven
                    for dim in 1:4
                        #@inbounds L += conj(branch.local_expansion[dim,nms]) * C_tmp
                        @inbounds L[1,dim] += branch.local_expansion[1,dim,nms] * C_tmp1 + branch.local_expansion[2,dim,nms] * C_tmp2
                        @inbounds L[2,dim] += branch.local_expansion[1,dim,nms] * C_tmp2 - branch.local_expansion[2,dim,nms] * C_tmp1
                    end
                end
                for m in 0:n
                    if n-j >= abs(m-k)
                        jnkm = (n - j) * (n - j) + n - j + m - k + 1
                        nms = (n * (n + 1)) >> 1 + m + 1
                        oddeven = odd_or_even((m-k) * (1 >> (m >= k)))
                        C_tmp1 = regular_harmonics[1,jnkm] * oddeven
                        C_tmp2 = regular_harmonics[2,jnkm] * oddeven
                        for dim in 1:4
                            #@inbounds L += conj(branch.local_expansion[dim,nms]) * C_tmp
                            @inbounds L[1,dim] += branch.local_expansion[1,dim,nms] * C_tmp1 - branch.local_expansion[2,dim,nms] * C_tmp2
                            @inbounds L[2,dim] += branch.local_expansion[1,dim,nms] * C_tmp2 + branch.local_expansion[2,dim,nms] * C_tmp1
                        end
                    end
                end
            end
            view(child.local_expansion,:,:,jks) .+= L
        end
    end
end

function L2B!(systems, branch::MultiBranch, derivatives_switches, expansion_order)
    for i in eachindex(systems)
		L2B!(systems[i], branch.bodies_index[i], branch.local_expansion, derivatives_switches[i], expansion_order, branch.center)
    end
end

function L2B!(system, branch::SingleBranch, derivatives_switch, expansion_order)
	L2B!(system, branch.bodies_index, branch.local_expansion, derivatives_switch, expansion_order, branch.center)
end

function L2B!(system, bodies_index, local_expansion, derivatives_switch::DerivativesSwitch{PS,VPS,VS,GS}, expansion_order, expansion_center) where {PS,VPS,VS,GS}
    for i_body in bodies_index
		scalar_potential, vector_potential, velocity, gradient = L2B(system[i_body,POSITION], expansion_center, local_expansion, derivatives_switch, expansion_order)
        if PS
            system[i_body,SCALAR_POTENTIAL] += scalar_potential
        end
        if VPS
            # note: system[i,VECTOR_POTENTIAL], system[i,VELOCITY], and system[i,VELOCITY_GRADIENT] must be mutable
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

@inline function L2B_scalar_potential(Lnm_real, Lnm_imag, eimp_real, eimp_imag, rn, Pnm, C_n_m)
    return complex_multiply_real(Lnm_real, Lnm_imag, eimp_real, eimp_imag) * rn * Pnm * C_n_m
end

@inline function L2B_vector_potential(Lnm_x_real::TF1, Lnm_x_imag::TF1, Lnm_y_real::TF1, Lnm_y_imag::TF1, Lnm_z_real::TF1, Lnm_z_imag::TF1, eimp_real::TF2, eimp_imag::TF2, rn, Pnm, C_n_m) where {TF1,TF2}
	TF = promote_type(TF1,TF2)
	return SVector{3,TF}(
        complex_multiply_real(Lnm_x_real, Lnm_x_imag, eimp_real, eimp_imag),
        complex_multiply_real(Lnm_y_real, Lnm_y_imag, eimp_real, eimp_imag),
        complex_multiply_real(Lnm_z_real, Lnm_z_imag, eimp_real, eimp_imag)
    ) * rn * Pnm * C_n_m
end

"Does not include the rotation matrix R yet. Eq 28 in S&L"
@inline function L2B_velocity(Lnm_real::TF1, Lnm_imag, Lnm_x_real, Lnm_x_imag, Lnm_y_real, Lnm_y_imag, Lnm_z_real, Lnm_z_imag, eimp_real::TF2, eimp_imag, rnm1, Pnm, dPdt, beta, C_n_m, R, n, m, derivatives_switch::DerivativesSwitch{PS,VPS,<:Any,<:Any}) where {TF1,TF2,PS,VPS}
	TF = promote_type(TF1,TF2)

    # initialize
    velocity = zero(SVector{3,TF})

    # intermediate quantities
    ux_real = n * Pnm
    uy_real = dPdt
    uz_imag = m * beta

    # rotate to cartesian
	u_real_cartesian, u_imag_cartesian = complex_multiply(R, ux_real, 0, uy_real, 0, 0, uz_imag)

    # velocity due to scalar potential

    Lnm_eimp_real, Lnm_eimp_imag = complex_multiply(eimp_real, eimp_imag, Lnm_real, Lnm_imag)

    velocity -= SVector{3}(
		complex_multiply_real(u_real_cartesian[1], u_imag_cartesian[1], Lnm_eimp_real, Lnm_eimp_imag),
		complex_multiply_real(u_real_cartesian[2], u_imag_cartesian[2], Lnm_eimp_real, Lnm_eimp_imag),
		complex_multiply_real(u_real_cartesian[3], u_imag_cartesian[3], Lnm_eimp_real, Lnm_eimp_imag)
    ) * rnm1 * C_n_m

    # velocity due to vector potential

    Lnm_x_eimp_real, Lnm_x_eimp_imag = complex_multiply(Lnm_x_real, Lnm_x_imag, eimp_real, eimp_imag)
	Lnm_y_eimp_real, Lnm_y_eimp_imag = complex_multiply(Lnm_y_real, Lnm_y_imag, eimp_real, eimp_imag)
    Lnm_z_eimp_real, Lnm_z_eimp_imag = complex_multiply(Lnm_z_real, Lnm_z_imag, eimp_real, eimp_imag)

	velocity += rnm1 * C_n_m * complex_cross_real(u_real_cartesian[1], u_imag_cartesian[1], u_real_cartesian[2], u_imag_cartesian[2], u_real_cartesian[3], u_imag_cartesian[3], Lnm_x_eimp_real, Lnm_x_eimp_imag, Lnm_y_eimp_real, Lnm_y_eimp_imag, Lnm_z_eimp_real, Lnm_z_eimp_imag)

    return velocity
end

"Eq 35-37 S&L"
@inline function L2B_velocity_gradient(Lnm_real::TF, Lnm_imag, Lnm_x_real, Lnm_x_imag, Lnm_y_real, Lnm_y_imag, Lnm_z_real, Lnm_z_imag, eimp_real, eimp_imag, rnm2, Pnm, dPdt, d2Pdt2, ddt_Pnm_st, alpha, beta, st, ct, sp, cp, C_n_m, R, n, m, derivatives_switch::DerivativesSwitch{PS,VPS,<:Any,<:Any}) where {TF,PS,VPS}
	ncheck = 3
	mcheck = 3
    # initialize outputs
    dudr = zero(SVector{3,TF})
    dudt_r = zero(SVector{3,TF})
    dudp_r_st = zero(SVector{3,TF})

    # intermediate values
    rnm2Cnm = rnm2 * C_n_m
    nm1 = n - 1

    #####
    ##### du/dr
    #####

    # intermediate quantities
    u1x_real = n * Pnm
    u1y_real = dPdt
    u1z_imag = m * beta

    vx_real, vx_imag = eimp_real * u1x_real, eimp_imag * u1x_real
    vy_real, vy_imag = eimp_real * u1y_real, eimp_imag * u1y_real
    vz_real, vz_imag = -eimp_imag * u1z_imag, eimp_real * u1z_imag

	# transform to cartesian
	u_real_cartesian, u_imag_cartesian = complex_multiply(R, vx_real, vx_imag, vy_real, vy_imag, vz_real, vz_imag)

    # due to scalar potential
    # if PS
        dudr -= nm1 * rnm2Cnm * SVector{3}(
			complex_multiply_real(u_real_cartesian[1], u_imag_cartesian[1], Lnm_real, Lnm_imag),
			complex_multiply_real(u_real_cartesian[2], u_imag_cartesian[2], Lnm_real, Lnm_imag),
			complex_multiply_real(u_real_cartesian[3], u_imag_cartesian[3], Lnm_real, Lnm_imag)
        )
    # end

    # due to vector potential
    # if VPS
		dudr += nm1 * rnm2Cnm * complex_cross_real(u_real_cartesian[1], u_imag_cartesian[1], u_real_cartesian[2], u_imag_cartesian[2], u_real_cartesian[3], u_imag_cartesian[3], Lnm_x_real, Lnm_x_imag, Lnm_y_real, Lnm_y_imag, Lnm_z_real, Lnm_z_imag)
    # end

    #####
    ##### du/dphi / r
    #####

    # intermediate quantities
    u2x_real = nm1 * dPdt
    u2y_real = d2Pdt2 + n * Pnm
    u2z_imag = m * ddt_Pnm_st

    vx_real, vx_imag = eimp_real * u2x_real, eimp_imag * u2x_real
    vy_real, vy_imag = eimp_real * u2y_real, eimp_imag * u2y_real
    vz_real, vz_imag = -eimp_imag * u2z_imag, eimp_real * u2z_imag

	u_real_cartesian, u_imag_cartesian = complex_multiply(R, vx_real, vx_imag, vy_real, vy_imag, vz_real, vz_imag)

    # due to scalar potential
    # if PS
        dudt_r -= rnm2Cnm * SVector{3}(
			complex_multiply_real(u_real_cartesian[1], u_imag_cartesian[1], Lnm_real, Lnm_imag),
			complex_multiply_real(u_real_cartesian[2], u_imag_cartesian[2], Lnm_real, Lnm_imag),
			complex_multiply_real(u_real_cartesian[3], u_imag_cartesian[3], Lnm_real, Lnm_imag)
        )
    # end

    # due to vector potential
    # if VPS
		dudt_r += rnm2Cnm * complex_cross_real(u_real_cartesian[1], u_imag_cartesian[1], u_real_cartesian[2], u_imag_cartesian[2], u_real_cartesian[3], u_imag_cartesian[3], Lnm_x_real, Lnm_x_imag, Lnm_y_real, Lnm_y_imag, Lnm_z_real, Lnm_z_imag)
    # end

    #####
    ##### du/dphi / r / sin(theta)
    #####

    # intermediate quantities
    Fnm = n * n * Pnm + d2Pdt2
    u3x_real = sp * Fnm
    u3x_imag = m * cp * (alpha * (m * m - 1) - Fnm)
    u3y_real = -cp * Fnm
    u3y_imag = m * sp * (alpha * (m * m - 1) - Fnm)
    u3z_real = zero(TF)
    u3z_imag = m * (n * ct * beta - dPdt)

    vx_real, vx_imag = complex_multiply(eimp_real, eimp_imag, u3x_real, u3x_imag)
    vy_real, vy_imag = complex_multiply(eimp_real, eimp_imag, u3y_real, u3y_imag)
    vz_real, vz_imag = -eimp_imag * u3z_imag, eimp_real * u3z_imag

    # due to scalar potential
    # if PS
        dudp_r_st -= SVector{3}(
            complex_multiply_real(vx_real, vx_imag, Lnm_real, Lnm_imag),
            complex_multiply_real(vy_real, vy_imag, Lnm_real, Lnm_imag),
			complex_multiply_real(vz_real, vz_imag, Lnm_real, Lnm_imag)
        ) * rnm2Cnm
    # end

    # due to vector potential
    # if VPS
        dudp_r_st += rnm2Cnm * complex_cross_real(vx_real, vx_imag, vy_real, vy_imag, vz_real, vz_imag, Lnm_x_real, Lnm_x_imag, Lnm_y_real, Lnm_y_imag, Lnm_z_real, Lnm_z_imag)
    # end

	return dudr, dudt_r, dudp_r_st
end

function L2B(body_position, expansion_center::SVector{3,TF}, local_expansion, derivatives_switch::DerivativesSwitch{PS,VPS,VS,GS}, expansion_order::Val{P}) where {TF,PS,VPS,VS,GS,P}
    #=
    R_n^m(\rho, \theta, \phi) &= (-1)^n \frac{\rho^n}{(n+\|m\|)!} P_n^{\|m\|}(\cos \theta) e^{i m \phi}\\
    I_n^m(\rho, \theta, \phi) &= (-1)^n \frac{(n-\|m\|)!}{\rho^{n+1}} P_n^{\|m\|}(\cos \theta) e^{i m \phi}\\
    phi = \sum \limits_{n=0}^P \sum \limits_{m=-n}^n L_n^m R_n^m
    =#
    dx, dy, dz = body_position - expansion_center
    r, theta, phi = cartesian_2_spherical(dx, dy, dz)

	#--- the following is inspired by Salloum and Lakkis (2020) with some improvements ---#

    # initialization
    scalar_potential = zero(TF)
    vector_potential = zero(SVector{3,TF})
    velocity = zero(SVector{3,TF})
    gradient = zero(SMatrix{3,3,TF,9})
    dudr = zero(SVector{3,TF})
    dudt_r = zero(SVector{3,TF})
    dudp_r_st = zero(SVector{3,TF})

    st, ct = sincos(theta)
    sp, cp = sincos(phi)

	if VS || GS
		R = SMatrix{3,3}(st*cp,st*sp,ct,ct*cp,ct*sp,-st,-sp,cp,0)
	end

    # note that by definition beta_n^m=0 for m<1 and alpha_n^m=0 for m<2

    # m = 0, n = 0
    n = m = 0
    rn = 1.0            # r^n
    rnm1 = 0.0          # r^(n-1)
    rnm2 = 0.0          # r^(n-2)
    C_n_m = 1.0           # (-1)^n / (n + |m|)!
    eimp_real = TF(1.0) # real(e^(i m phi))
    eimp_imag = TF(0.0) # imag(e^(i m phi))
    Pnm = 1.0           # associated legendre polynomial of degree n, order m
    alpha_n_m = 0.0     # alpha = Pnm / sin^2(theta)
    beta_n_m = 0.0      # beta = Pnm / sin(theta)
    dPdt_n_m = 0.0      # derivative w.r.t. theta of P_n_m
    d2Pdt2_n_m = 0.0    # second derivative w.r.t. theta of P_n_m
    index = 1           # index of the local expansion corresponding to n and m

    if PS
		scalar_potential += L2B_scalar_potential(local_expansion[1,1,index], local_expansion[2,1,index], eimp_real, eimp_imag, rn, Pnm, C_n_m)
		# scalar_potential += local_expansion[1,1,index]
    end

    if VPS
        vector_potential += L2B_vector_potential(local_expansion[1,2,index], local_expansion[2,2,index], local_expansion[1,3,index], local_expansion[2,3,index], local_expansion[1,4,index], local_expansion[2,4,index], eimp_real, eimp_imag, rn, Pnm, C_n_m)
        # vector_potential += SVector{3}(
        #     local_expansion[1,2,1],
        #     local_expansion[1,3,1],
        #     local_expansion[1,4,1]
        # )
    end

    #--- m = 0, n = 1 ---#

    if P > 0
        n = 1
        m = 0

        #####
        ##### recurse
        #####

        # r
        rnm2 = rnm1
        rnm1 = rn
        rn *= r

        # Pnm
        Pnm = ct
        dPdt_n_m = -st
        d2Pdt2_n_m = -ct

        # C_n_m
        C_n_m *= -1

        # eimp is unchanged
        # alpha = 0
        beta_n_m = zero(TF)

        # index
        index += 1

        #####
        ##### evaluate
        #####

        if PS
            scalar_potential += L2B_scalar_potential(local_expansion[1,1,index], local_expansion[2,1,index], eimp_real, eimp_imag, rn, Pnm, C_n_m)
        end

        if VPS
            vector_potential += L2B_vector_potential(local_expansion[1,2,index], local_expansion[2,2,index], local_expansion[1,3,index], local_expansion[2,3,index], local_expansion[1,4,index], local_expansion[2,4,index], eimp_real, eimp_imag, rn, Pnm, C_n_m)
        end

        if VS
            velocity += L2B_velocity(local_expansion[1,1,index], local_expansion[2,1,index], local_expansion[1,2,index], local_expansion[2,2,index], local_expansion[1,3,index], local_expansion[2,3,index], local_expansion[1,4,index], local_expansion[2,4,index], eimp_real, eimp_imag, rnm1, Pnm, dPdt_n_m, beta_n_m, C_n_m, R, n, m, derivatives_switch)
        end

        #--- m = 1, n = 1 ---#

        #####
        ##### recurse
        #####

        n = m = 1

        # beta
        beta_n_m = -1.0
        beta_nm1_m = 0.0
        beta_nm2_m = 0.0

        # Pnm
        Pnm = -st
        dPdt_n_m = -ct
        d2Pdt2_n_m = st

        # C_n_m
        n_plus_m = n + m
        C_n_m /= n_plus_m

        # eimp
        eimp_real, eimp_imag = complex_multiply(eimp_real, eimp_imag, cp, sp)

        # index
        index += 1

        # r is unchanged

        #####
        ##### evaluate
        #####

        if PS
            scalar_potential += 2 * L2B_scalar_potential(local_expansion[1,1,index], local_expansion[2,1,index], eimp_real, eimp_imag, rn, Pnm, C_n_m)
        end

        if VPS
            vector_potential += 2 * L2B_vector_potential(local_expansion[1,2,index], local_expansion[2,2,index], local_expansion[1,3,index], local_expansion[2,3,index], local_expansion[1,4,index], local_expansion[2,4,index], eimp_real, eimp_imag, rn, Pnm, C_n_m)
        end

        if VS
            velocity += 2 * L2B_velocity(local_expansion[1,1,index], local_expansion[2,1,index], local_expansion[1,2,index], local_expansion[2,2,index], local_expansion[1,3,index], local_expansion[2,3,index], local_expansion[1,4,index], local_expansion[2,4,index], eimp_real, eimp_imag, rnm1, Pnm, dPdt_n_m, beta_n_m, C_n_m, R, n, m, derivatives_switch)
        end

    end

    #--- m = 0, m = 1, n > 1 ---#

    if P > 1

        m = 1
        alpha_n_0 = zero(TF)        # n>1, m<2
        beta_n_0 = zero(TF)         # n>1, m<1
        P_nm2_0 = 0.0               # n=-1,m=0
        P_nm1_0 = 1.0               # n=0, m=0
        P_n_0 = ct                  # n=1, m=0
        dPdt_nm1_0 = zero(TF)       # n=0, m=0
        dPdt_n_0 = -st              # n=1, m=0
        d2Pdt2_nm1_0 = 0.0          # n=0, m=0
        d2Pdt2_n_0 = dPdt_n_m       # n=1, m=0
        d3Pdt3_nm2_0 = 0.0          # n=-1,m=0
        d3Pdt3_nm1_0 = 0.0          # n=0, m=0
        d3Pdt3_n_0 = st             # n=1, m=0
        dPdct_n_0 = 1.0             # n=1, m=0
        d2Pdct2_n_0 = 0.0           # n=1, m=0
        dPdct_nm1_0 = 0.0           # n=0, m=0
        d2Pdct2_nm1_0 = 0.0         # n=0, m=0
        P_nm1_m = 0.0               # n=0, m=1

        for n in 2:P

            #####
            ##### recurse
            #####

            # rn
            rnm2 = rnm1
            rnm1 = rn
            rn *= r

            # C_n_m
            n_plus_m += 1
            C_n_0 = -C_n_m
            C_n_m /= -n_plus_m

            # beta
            beta_nm2_m = beta_nm1_m
            beta_nm1_m = beta_n_m
            beta_n_m = ((2*n-1) * ct * beta_nm1_m - n * beta_nm2_m) / (n-1)

            # Pnm, m=1
            P_nm1_m = Pnm
            Pnm = beta_n_m * st

            # Pnm, m=0
            P_nm2_0 = P_nm1_0
            P_nm1_0 = P_n_0
            P_n_0 = 1/n*((2*n-1) * ct * P_nm1_0 - (n-1) * P_nm2_0)

            # first derivatives, m=1
            dPdt_n_m = n * ct * beta_n_m - (n + m) * beta_nm1_m

            # first derivatives, m=0
            dPdt_nm1_0 = dPdt_n_0
            dPdt_n_0 = ct * dPdt_nm1_0 - n * st * P_nm1_0

            # second derivatives
            d2Pdt2_nm1_0 = d2Pdt2_n_0
            d2Pdt2_n_0 = dPdt_n_m
            d3Pdt3_nm2_0 = d3Pdt3_nm1_0
            d3Pdt3_nm1_0 = d3Pdt3_n_0
            d3Pdt3_n_0 = d3Pdt3_nm2_0 - (2*n-1) * (st*d2Pdt2_nm1_0 + 2*ct*P_nm1_m - st*P_nm1_0)
            d2Pdt2_n_m = d3Pdt3_n_0

            # second derivatives over sin(theta)
            dPdct_nm1_0 = dPdct_n_0
            d2Pdct2_nm1_0 = d2Pdct2_n_0
            dPdct_n_0 = n * P_nm1_0 + ct * dPdct_nm1_0
            d2Pdct2_n_0 = (n+1) * dPdct_nm1_0 + ct * d2Pdct2_nm1_0
            ddt_Pnm_st = st * d2Pdct2_n_0

            # index
            index = (n * (n + 1)) >> 1 + m + 1
			index_m0 = (n * (n + 1)) >> 1 + 0 + 1

            # eimp is unchanged

            #####
            ##### evaluate
            #####

            if PS
				scalar_potential += L2B_scalar_potential(local_expansion[1,1,index_m0], local_expansion[2,1,index_m0], one(TF), zero(TF), rn, P_n_0, C_n_0)
                scalar_potential += 2 * L2B_scalar_potential(local_expansion[1,1,index], local_expansion[2,1,index], eimp_real, eimp_imag, rn, Pnm, C_n_m)
            end

            if VPS
				vector_potential += L2B_vector_potential(local_expansion[1,2,index_m0], local_expansion[2,2,index_m0], local_expansion[1,3,index_m0], local_expansion[2,3,index_m0], local_expansion[1,4,index_m0], local_expansion[2,4,index_m0], one(TF), zero(TF), rn, P_n_0, C_n_0)
                vector_potential += 2 * L2B_vector_potential(local_expansion[1,2,index], local_expansion[2,2,index], local_expansion[1,3,index], local_expansion[2,3,index], local_expansion[1,4,index], local_expansion[2,4,index], eimp_real, eimp_imag, rn, Pnm, C_n_m)
            end

            if VS
                velocity += L2B_velocity(local_expansion[1,1,index_m0], local_expansion[2,1,index_m0], local_expansion[1,2,index_m0], local_expansion[2,2,index_m0], local_expansion[1,3,index_m0], local_expansion[2,3,index_m0], local_expansion[1,4,index_m0], local_expansion[2,4,index_m0], one(TF), zero(TF), rnm1, P_n_0, dPdt_n_0, beta_n_0, C_n_0, R, n, 0, derivatives_switch)
                velocity += 2 * L2B_velocity(local_expansion[1,1,index], local_expansion[2,1,index], local_expansion[1,2,index], local_expansion[2,2,index], local_expansion[1,3,index], local_expansion[2,3,index], local_expansion[1,4,index], local_expansion[2,4,index], eimp_real, eimp_imag, rnm1, Pnm, dPdt_n_m, beta_n_m, C_n_m, R, n, m, derivatives_switch)
            end

            if GS
                # m = 0
                ddt_Pn0_st = zero(TF)
                this_dudr, this_dudt_r, this_dudp_r_st = L2B_velocity_gradient(local_expansion[1,1,index_m0], local_expansion[2,1,index_m0], local_expansion[1,2,index_m0], local_expansion[2,2,index_m0], local_expansion[1,3,index_m0], local_expansion[2,3,index_m0], local_expansion[1,4,index_m0], local_expansion[2,4,index_m0], one(TF), zero(TF), rnm2, P_n_0, dPdt_n_0, d2Pdt2_n_0, ddt_Pn0_st, alpha_n_0, beta_n_0, st, ct, sp, cp, C_n_0, R, n, 0, derivatives_switch)

                dudr += this_dudr
                dudt_r += this_dudt_r
                dudp_r_st += this_dudp_r_st

                # m = 1

                this_dudr, this_dudt_r, this_dudp_r_st = L2B_velocity_gradient(local_expansion[1,1,index], local_expansion[2,1,index], local_expansion[1,2,index], local_expansion[2,2,index], local_expansion[1,3,index], local_expansion[2,3,index], local_expansion[1,4,index], local_expansion[2,4,index], eimp_real, eimp_imag, rnm2, Pnm, dPdt_n_m, d2Pdt2_n_m, ddt_Pnm_st, alpha_n_0, beta_n_m, st, ct, sp, cp, C_n_m, R, n, m, derivatives_switch)

				dudr += 2 * this_dudr
                dudt_r += 2 * this_dudt_r
                dudp_r_st += 2 * this_dudp_r_st
            end

        end

        #--- m > 1, n > 1 ---#

        # double factorial
        # double_factorial = 3

        # sin^(n-2)(theta)
        # st_nm2 = one(TF)

		# alpha_n_n
		alpha_n_n = TF(3.0) # sin(theta)^(n-2) (2*n-1)!! for n=2

        # rn
        rnm2 = zero(TF) # n=-1...
        rnm1 = one(TF)  # n=0
        rn = r          # n=1

        # C_n_m: n = m = 1
        n_plus_m = 2
        C_n_m = -0.5 # (-1)^1 / (1 + |1|)!

        for m in 2:P

            #####
            ##### recurse
            #####

            n = m

            # rn
            rnm2 = rnm1
            rnm1 = rn
            rn *= r

            # C_n_m: increment n
            n_plus_m += 1
            C_n_m /= -n_plus_m

            # increment m
            n_plus_m += 1
            C_n_m /= n_plus_m

            # alpha_n_n = (-1)^n * st^(n-2) * double_factorial(2*n-1)
			alpha_n_m = alpha_n_n
            alpha_nm1_m = 0.0
            alpha_nm2_m = 0.0

            # beta
            beta_n_m = alpha_n_m * st
            beta_nm1_m = 0.0

            # P_n_m
            P_n_m = beta_n_m * st

            # first derivative
            dPdt_n_m = n * ct * beta_n_m - (n + m) * beta_nm1_m

            # second derivative
            d2Pdt2_n_m = (n + m) * ct * alpha_nm1_m - (n * (n * st^2 + 1) - m^2) * alpha_n_m

            # derivative of Pnm/sin(theta)
            ddt_Pnm_st = (n-1) * ct * alpha_n_m - (n+m) * alpha_nm1_m

            # eimp
            eimp_real, eimp_imag = complex_multiply(eimp_real, eimp_imag, cp, sp)

            # index
            index = (n * (n + 1)) >> 1 + m + 1

            #####
            ##### evaluate
            #####

            if PS
				scalar_potential += 2 * L2B_scalar_potential(local_expansion[1,1,index], local_expansion[2,1,index], eimp_real, eimp_imag, rn, P_n_m, C_n_m)
			end

            if VPS
                vector_potential += 2 * L2B_vector_potential(local_expansion[1,2,index], local_expansion[2,2,index], local_expansion[1,3,index], local_expansion[2,3,index], local_expansion[1,4,index], local_expansion[2,4,index], eimp_real, eimp_imag, rn, P_n_m, C_n_m)
            end

            if VS
                velocity += 2 * L2B_velocity(local_expansion[1,1,index], local_expansion[2,1,index], local_expansion[1,2,index], local_expansion[2,2,index], local_expansion[1,3,index], local_expansion[2,3,index], local_expansion[1,4,index], local_expansion[2,4,index], eimp_real, eimp_imag, rnm1, P_n_m, dPdt_n_m, beta_n_m, C_n_m, R, n, m, derivatives_switch)
            end

            if GS
                this_dudr, this_dudt_r, this_dudp_r_st = L2B_velocity_gradient(local_expansion[1,1,index], local_expansion[2,1,index], local_expansion[1,2,index], local_expansion[2,2,index], local_expansion[1,3,index], local_expansion[2,3,index], local_expansion[1,4,index], local_expansion[2,4,index], eimp_real, eimp_imag, rnm2, P_n_m, dPdt_n_m, d2Pdt2_n_m, ddt_Pnm_st, alpha_n_m, beta_n_m, st, ct, sp, cp, C_n_m, R, n, m, derivatives_switch)

				dudr += 2 * this_dudr
                dudt_r += 2 * this_dudt_r
                dudp_r_st += 2 * this_dudp_r_st
            end

			# prepare to recurse C_n_m
			C_n_m_loop = C_n_m
			n_plus_m_loop = n_plus_m
			rn_loop = rn
			rnm1_loop = rnm1
			rnm2_loop = rnm2

            for n in m+1:P

                #####
                ##### recurse
                #####

				# rn
				rnm2_loop = rnm1_loop
            	rnm1_loop = rn_loop
            	rn_loop *= r

                # alpha
                alpha_nm2_m = alpha_nm1_m
                alpha_nm1_m = alpha_n_m
                alpha_n_m = 1/(n-m) * ( (2*n-1)*ct*alpha_nm1_m - (n+m-1)*alpha_nm2_m )

                # beta
                beta_nm1_m = beta_n_m
                beta_n_m = alpha_n_m * st

                # P_n_m
                P_n_m = beta_n_m * st

                # first derivative
                dPdt_n_m = n * ct * beta_n_m - (n + m) * beta_nm1_m

                # second derivative
                d2Pdt2_n_m = (n + m) * ct * alpha_nm1_m - (n * (n * st^2 + 1) - m^2) * alpha_n_m

                # derivative of Pnm/sin(theta)
                ddt_Pnm_st = (n-1) * ct * alpha_n_m - (n+m) * alpha_nm1_m

				# C_n_m
				n_plus_m_loop += 1
				C_n_m_loop /= -n_plus_m_loop

				# index
            	index = (n * (n + 1)) >> 1 + m + 1

                #####
                ##### evaluate
                #####
                if PS
                    scalar_potential += 2 * L2B_scalar_potential(local_expansion[1,1,index], local_expansion[2,1,index], eimp_real, eimp_imag, rn_loop, P_n_m, C_n_m_loop)
                end

                if VPS
                    vector_potential += 2 * L2B_vector_potential(local_expansion[1,2,index], local_expansion[2,2,index], local_expansion[1,3,index], local_expansion[2,3,index], local_expansion[1,4,index], local_expansion[2,4,index], eimp_real, eimp_imag, rn_loop, P_n_m, C_n_m_loop)
                end

                if VS
                    velocity += 2 * L2B_velocity(local_expansion[1,1,index], local_expansion[2,1,index], local_expansion[1,2,index], local_expansion[2,2,index], local_expansion[1,3,index], local_expansion[2,3,index], local_expansion[1,4,index], local_expansion[2,4,index], eimp_real, eimp_imag, rnm1_loop, P_n_m, dPdt_n_m, beta_n_m, C_n_m_loop, R, n, m, derivatives_switch)
                end

                if GS
                    this_dudr, this_dudt_r, this_dudp_r_st = L2B_velocity_gradient(local_expansion[1,1,index], local_expansion[2,1,index], local_expansion[1,2,index], local_expansion[2,2,index], local_expansion[1,3,index], local_expansion[2,3,index], local_expansion[1,4,index], local_expansion[2,4,index], eimp_real, eimp_imag, rnm2_loop, P_n_m, dPdt_n_m, d2Pdt2_n_m, ddt_Pnm_st, alpha_n_m, beta_n_m, st, ct, sp, cp, C_n_m_loop, R, n, m, derivatives_switch)

                    dudr += 2 * this_dudr
                    dudt_r += 2 * this_dudt_r
                    dudp_r_st += 2 * this_dudp_r_st
                end
            end

            # tail recursion
			alpha_n_n *= -st * (2*n+1)
			# double_factorial *= 2*n+1
            # st_nm2 *= -st

        end

    end

    # rotate to cartesian
    # R = SMatrix{3,3}(st*cp,st*sp,ct,ct*cp,ct*sp,-st,-sp,cp,0)
    duidrk = hcat(dudr, dudt_r, dudp_r_st)
    gradient = duidrk * R'

	return scalar_potential, vector_potential, velocity, gradient
end
