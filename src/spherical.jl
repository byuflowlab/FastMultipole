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

# Defining these outside for readability and for making it easier to pass ForwardDiff through.
function get_drjdxi(rho,s_theta,s_phi,c_theta,c_phi)
    return @SMatrix [s_theta*c_phi    c_theta*c_phi/rho      -s_phi/rho/s_theta
            s_theta * s_phi  c_theta * s_phi / rho  c_phi / rho / s_theta
            c_theta         -s_theta / rho          0                    ]
end

# replaced repeated sin/cos calls with precomputed values. Reduces ReverseDiff allocations due to fewer functions recorded to the tape. Also increases readability.
function get_drkdxidxj(rho,s_theta,s_phi,c_theta,c_phi,k_coord)

    if k_coord == 1 # r coordinate
        return @SMatrix [
            (1-c_phi^2 * s_theta^2)/rho -s_theta^2*c_phi*s_phi/rho -s_theta*c_phi*c_theta/rho;
            (-s_theta^2*c_phi*s_phi)/rho (1-s_theta^2*s_phi^2)/rho -s_theta*s_phi*c_theta/rho;
            -s_theta*c_phi*c_theta/rho -s_theta*s_phi*c_theta/rho s_theta^2/rho
        ]
    elseif k_coord == 2 # theta coordinate
        return @SMatrix [
            c_theta/s_theta*(1-c_phi^2*(1+2*s_theta^2))/rho^2 -c_theta/s_theta*s_phi*c_phi*(1+2*s_theta^2)/rho^2 c_phi*(1-2*c_theta^2)/rho^2;
            -c_theta/s_theta*s_phi*c_phi*(1+2*s_theta^2)/rho^2 c_theta/s_theta*(1-s_phi^2*(1+2*s_theta^2))/rho^2 (2*s_theta^2-1)/rho^2*s_phi;
            c_phi*(1-2*c_theta^2)/rho^2 (2*s_theta^2-1)/rho^2*s_phi 2*s_theta*c_theta/rho^2
        ]
    else # phi coordinate
        return @SMatrix [
            2*c_phi*s_phi/rho^2/s_theta^2 (2*s_phi^2-1)/rho^2/s_theta^2 0;
            (2*s_phi^2-1)/rho^2/s_theta^2 -2*s_phi*c_phi/rho^2/s_theta^2 0;
            0 0 0
        ]
    end

end

function s2c_jac!(potential_jacobian, workspace, rho, theta, phi)

    s_theta, c_theta = sincos(theta)
    s_phi, c_phi = sincos(phi)

    drjdxi = get_drjdxi(rho,s_theta,s_phi,c_theta,c_phi)

    workspace .= potential_jacobian

    mul!(potential_jacobian, drjdxi, workspace)

    return potential_jacobian

end

#=function s2c_hess!(potential_jacobian, potential_hessian, workspace, rho, theta, phi)

    s_theta, c_theta = sincos(theta)
    s_phi, c_phi = sincos(phi)

    R = get_drjdxi(rho,s_theta,s_phi,c_theta,c_phi)
    RR = get_RR(rho,theta,phi)
    workspace3x3 = view(workspace,:,1:3)
    for ind in 1:4
        workspace3x3 .= potential_hessian[:,:,ind]
        potential_hessian[:,:,ind] .= R * workspace3x3 * transpose(R)
        for ν in 1:3
            potential_hessian[:,:,ind] .+= RR[:,:,ν] * potential_jacobian[ν,ind]
        end
    end
    return potential_hessian

end=#

# I split these so that I could add custom reverse diff rules to each part separately. Doing it with one function is an enormous headache; it should be possible but there is no documentation on how to handle functions that return tuples.
function spherical_2_cartesian!(potential_jacobian, potential_hessian, workspace, rho, theta, phi)
    # the order of these two matters, since the hessian needs to use the original jacobian.
    potential_hessian = s2c_hess!(potential_jacobian, potential_hessian, workspace, rho, theta, phi)
    potential_jacobian = s2c_jac!(potential_jacobian, workspace, rho, theta, phi)

    return potential_jacobian,potential_hessian
end

# will need to make use of this for more efficient computations, but that will wait until I tackle rrules for functions with multiple returned objects.
function s2c_hess!(potential_jacobian, potential_hessian, workspace, rho, theta, phi)
    # get partial derivatives of the coordinates
    s_theta, c_theta = sincos(theta)
    s_phi, c_phi = sincos(phi)
    
    drjdxi = @SMatrix [
        s_theta*c_phi c_theta*c_phi/rho -s_phi/rho/s_theta
        s_theta * s_phi c_theta * s_phi / rho c_phi / rho / s_theta
        c_theta -s_theta / rho 0.0
    ]
    # convert Hessian to cartesian coordinates
    workspace3x3 = view(workspace,:,1:3)
    for ind in 1:4
        # potential_hessian[:,:,ind] .= drjdxi * potential_hessian[:,:,ind] * transpose(drjdxi)
        mul!(workspace3x3, drjdxi, view(potential_hessian,:,:,ind))
        mul!(view(potential_hessian,:,:,ind), workspace3x3, transpose(drjdxi))
    end
    for k_coord in 1:3 # loop over r, theta, and phi to save me some space on the stack
        if k_coord == 1 # r coordinate
            drkdxidxj = @SMatrix [
                (1-c_phi^2 * s_theta^2)/rho -s_theta^2*c_phi*s_phi/rho -s_theta*c_phi*c_theta/rho;
                (-s_theta^2*c_phi*s_phi)/rho (1-s_theta^2*s_phi^2)/rho -s_theta*s_phi*c_theta/rho;
                -s_theta*c_phi*c_theta/rho -s_theta*s_phi*c_theta/rho s_theta^2/rho
            ]
        elseif k_coord == 2 # theta coordinate
            drkdxidxj = @SMatrix [
                c_theta/s_theta*(1-c_phi^2*(1+2*s_theta^2))/rho^2 -c_theta/s_theta*s_phi*c_phi*(1+2*s_theta^2)/rho^2 c_phi*(1-2*c_theta^2)/rho^2;
                -c_theta/s_theta*s_phi*c_phi*(1+2*s_theta^2)/rho^2 c_theta/s_theta*(1-s_phi^2*(1+2*s_theta^2))/rho^2 (2*s_theta^2-1)/rho^2*s_phi;
                c_phi*(1-2*c_theta^2)/rho^2 (2*s_theta^2-1)/rho^2*s_phi 2*s_theta*c_theta/rho^2
            ]
        else # phi coordinate
            drkdxidxj = @SMatrix [
                2*c_phi*s_phi/rho^2/s_theta^2 (2*s_phi^2-1)/rho^2/s_theta^2 0;
                (2*s_phi^2-1)/rho^2/s_theta^2 -2*s_phi*c_phi/rho^2/s_theta^2 0;
                0 0 0
            ]
        end
        for ind in 1:4
            view(potential_hessian,:,:,ind) .+= drkdxidxj * potential_jacobian[k_coord,ind]
        end
    end
    
    return potential_hessian
end

function flatten_jacobian!(jacobian)
    jacobian[1,1] = -jacobian[1,1] + jacobian[2,4] - jacobian[3,3]
    jacobian[2,1] = -jacobian[2,1] + jacobian[3,2] - jacobian[1,4]
    jacobian[3,1] = -jacobian[3,1] + jacobian[1,3] - jacobian[2,2]
    return jacobian
end

function flatten_hessian!(hessian)
    hessian[1,1,1] = -hessian[1,1,1]+hessian[2,1,4]-hessian[3,1,3]
    hessian[2,1,1] = -hessian[2,1,1]+hessian[3,1,2]-hessian[1,1,4]
    hessian[3,1,1] = -hessian[3,1,1]+hessian[1,1,3]-hessian[2,1,2]
    hessian[1,2,1] = -hessian[1,2,1]+hessian[2,2,4]-hessian[3,2,3]
    hessian[2,2,1] = -hessian[2,2,1]+hessian[3,2,2]-hessian[1,2,4]
    hessian[3,2,1] = -hessian[3,2,1]+hessian[1,2,3]-hessian[2,2,2]
    hessian[1,3,1] = -hessian[1,3,1]+hessian[2,3,4]-hessian[3,3,3]
    hessian[2,3,1] = -hessian[2,3,1]+hessian[3,3,2]-hessian[1,3,4]
    hessian[3,3,1] = -hessian[3,3,1]+hessian[1,3,3]-hessian[2,3,2]
    return hessian
end

function flatten_derivatives!(jacobian, hessian)
    # velocity
    jacobian .= flatten_jacobian!(jacobian)
    # velocity gradient
    hessian .= flatten_hessian!(hessian)
end

@inline function update_scalar_potential(scalar_potential,LE,h,P)

    for n in 0:P
        # nm = n * n + n + 1 # m = 0
        nms = (n * (n+1)) >> 1 + 1 # m = 0
        scalar_potential += real(LE[1,nms] * h[nms])
        for m in 1:n # m > 0
            # nm = n * n + n + m + 1
            nms = (n * (n + 1)) >> 1 + m + 1
            scalar_potential += 2 * real(LE[1,nms] * h[nms])
        end
    end
    return scalar_potential

end

@inline function update_vector_potential!(vector_potential,LE,h,P)

    for n in 0:P
        # nm = n * n + n + 1 # m = 0
        nms = (n * (n+1)) >> 1 + 1 # m = 0
        vector_potential[1] += real(LE[2,nms] * h[nms])
        vector_potential[2] += real(LE[3,nms] * h[nms])
        vector_potential[3] += real(LE[4,nms] * h[nms])
        for m in 1:n # m > 0
            # nm = n * n + n + m + 1
            nms = (n * (n + 1)) >> 1 + m + 1
            vector_potential[1] += 2 * real(LE[2,nms] * h[nms])
            vector_potential[2] += 2 * real(LE[3,nms] * h[nms])
            vector_potential[3] += 2 * real(LE[4,nms] * h[nms])
        end
    end
    return vector_potential
end

@inline function update_potential_jacobian!(potential_jacobian,LE,h,ht,P,r)

    for n in 0:P
        # nm = n * n + n + 1 # m = 0
        nms = (n * (n+1)) >> 1 + 1 # m = 0
        for ind in 1:4
            # store derivatives of the potential in spherical coordinates here
            potential_jacobian[1,ind] += n/r * real(LE[ind,nms] * h[nms]) # dPsi/dr
            potential_jacobian[2,ind] += real(LE[ind,nms] * ht[nms]) # dPsi/dtheta
            # dJ_potential[3,ind] += 0 # dPsi/dphi
        end
        for m in 1:n # m > 0
            # nm = n * n + n + m + 1
            nms = (n * (n + 1)) >> 1 + m + 1
            for ind in 1:4
                # store derivatives of the potential in spherical harmonics here
                potential_jacobian[1,ind] += 2 * n/r * real(LE[ind,nms] * h[nms]) # dPsi/dr
                potential_jacobian[2,ind] += 2 * real(LE[ind,nms] * ht[nms]) # dPsi/dtheta
                potential_jacobian[3,ind] += 2 * m * real(im * LE[ind,nms] * h[nms]) # dPsi/dphi
            end
        end
    end
    return potential_jacobian
end

@inline function update_potential_hessian!(potential_hessian,LE,h,ht,ht2,P,r)

    for n in 0:P
        # nm = n * n + n + 1 # m = 0
        nms = (n * (n+1)) >> 1 + 1 # m = 0
        for ind in 1:4
            # store derivatives of the potential in spherical coordinates here
            potential_hessian[1,1,ind] += n * (n-1) / r^2 * real(LE[ind,nms] * h[nms]) # d2Psi/dr2
            potential_hessian[2,1,ind] += n/r * real(LE[ind,nms] * ht[nms]) # d2Psi/dtheta dr
            # potential_hessian[3,1,ind] += 0 # d2Psi/dphi dr
            potential_hessian[1,2,ind] += n/r * real(LE[ind,nms] * ht[nms]) # d2Psi/dr dtheta
            potential_hessian[2,2,ind] += real(LE[ind,nms] * ht2[nms]) # d2Psi/dtheta2
            # potential_hessian[3,2,ind] += 0 # d2Psi/dphi dtheta
            # potential_hessian[1,3,ind] += 0 # d2Psi/dr dphi
            # potential_hessian[2,3,ind] += 0 # d2Psi/dtheta dphi
            # potential_hessian[3,3,ind] += 0 # d2Psi/dphi2
        end
        for m in 1:n # m > 0
            # nm = n * n + n + m + 1
            nms = (n * (n + 1)) >> 1 + m + 1
            for ind in 1:4
                # store derivatives of the potential in spherical harmonics here
                potential_hessian[1,1,ind] += 2 * n * (n-1) / r^2 * real(LE[ind,nms] * h[nms]) # d2Psi/dr2
                potential_hessian[2,1,ind] += 2 * n/r * real(LE[ind,nms] * ht[nms]) # d2Psi/dtheta dr
                potential_hessian[3,1,ind] += 2 * n * m / r * real(im * LE[ind,nms] * h[nms]) # d2Psi/dphi dr
                potential_hessian[1,2,ind] += 2 * n/r * real(LE[ind,nms] * ht[nms]) # d2Psi/dr dtheta
                potential_hessian[2,2,ind] += 2 * real(LE[ind,nms] * ht2[nms]) # d2Psi/dtheta2
                potential_hessian[3,2,ind] += 2 * m * real(im * LE[ind,nms] * ht[nms]) # d2Psi/dphi dtheta
                potential_hessian[1,3,ind] += 2 * n * m / r * real(im * LE[ind,nms] * h[nms]) # d2Psi/dr dphi
                potential_hessian[2,3,ind] += 2 * m * real(im * LE[ind,nms] * ht[nms]) # d2Psi/dtheta dphi
                potential_hessian[3,3,ind] += 2 * -m^2 * real(LE[ind,nms] * h[nms]) # d2Psi/dphi2
            end
        end
    end
    return potential_hessian
end

@inline odd_or_even(n::Int) = (n & 1) == 1 ? -1 : 1

@inline ipow2l(n::Int) = n >= 0 ? 1 : odd_or_even(n);

function regular_harmonic_derivs!(harmonics, rho, theta, phi, P)
    y,x = sincos(theta)
    invY = y == 0 ? 0 : 1 / y
    fact = 1.0
    pl = 1.0
    rhom = 1.0
    ei = exp(im * phi)
    eim = 1.0+0.0im
    for m=0:P
        p = pl
        lpl = (m * (m + 1)) >> 1 + m + 1
        # lpl = m * m + 2 * m + 1
        # lml = m * m + 1
        harmonics[1,lpl] = rhom * p * eim
        # harmonics[lml] = conj(harmonics[lpl])
        p1 = p
        p = x * (2 * m + 1) * p1
        harmonics[2,lpl] = rhom * (p - (m + 1) * x * p1) * invY * eim
        harmonics[3,lpl] = rhom * (-x * p + (-m + (m+1)^2 * x^2) * p1) * invY^2 * eim

        rhom *= rho
        rhol = rhom
        for l=m+1:P
            lpm = (l * (l + 1)) >> 1 + m + 1
            # lpm = l * l + l + m + 1
            # lmm = l * l + l - m + 1
            rhol /= -(l + m)
            harmonics[1,lpm] = rhol * p * eim
            # harmonics[lmm] = conj(harmonics[lpm])
            p2 = p1
            p1 = p
            p = (x * (2 * l + 1) * p1 - (l + m) * p2) / (l - m + 1)
            harmonics[2,lpm] = rhol * ((l - m + 1) * p - (l + 1) * x * p1) * invY * eim
            harmonics[3,lpm] = rhol * ((m-l-1) * x * p + (m^2 - l*(l+1) + (l+1)^2 * x^2) * p1) * invY^2 * eim
            rhol *= rho
        end
        rhom /= -(2 * m + 2) * (2 * m + 1)
        pl = -pl * fact * y
        fact += 2
        eim *= ei

    end
    return harmonics
end

function regular_harmonic!(harmonics, rho, theta, phi, P)
    harmonics .= zero(eltype(harmonics)) # temp
    y,x = sincos(theta)
    fact = 1.0
    pl = 1.0
    rhom = 1.0 # rho^l / (l+m)! * (-1)^l
    ei = exp(im * phi)
    eim = 1.0
    for m=0:P # l=m up here
        p = pl
        lpl = m * m + 2 * m + 1
        lml = m * m + 1
        harmonics[lpl] = rhom * p * eim
        harmonics[lml] = conj(harmonics[lpl])
        p1 = p
        p = x * (2 * m + 1) * p1
        rhom *= rho
        rhol = rhom
        for l=m+1:P # l>m in here
            lpm = l * l + l + m + 1
            lmm = l * l + l - m + 1
            rhol /= -(l + m)
            harmonics[lpm] = rhol * p * eim
            harmonics[lmm] = conj(harmonics[lpm])
            p2 = p1
            p1 = p
            p = (x * (2 * l + 1) * p1 - (l + m) * p2) / (l - m + 1)
            rhol *= rho
        end
        rhom /= -(2 * m + 2) * (2 * m + 1)
        pl = -pl * fact * y
        fact += 2
        eim *= ei
    end
    return harmonics
end

function irregular_harmonic!(harmonics, rho, theta, phi, P)
    y, x = sincos(theta)
    fact = 1
    pl = 1
    invR = -1.0 / rho
    rhom = -invR
    ei = exp(im * phi)
    eim = 1.0
    for m=0:P
        p = pl
        npl = m * m + 2 * m + 1
        nml = m * m + 1
        harmonics[npl] = rhom * p * eim
        harmonics[nml] = conj(harmonics[npl])
        p1 = p
        p = x * (2 * m + 1) * p1
        rhom *= invR
        rhon = rhom
        for l=m+1:P
            npm = l * l + l + m + 1
            nmm = l * l + l - m + 1
            # npm_max = P^2 + P + P + 1 = (P+1)^2
            harmonics[npm] = rhon * p * eim
            harmonics[nmm] = conj(harmonics[npm])
            p2 = p1
            p1 = p
            p = (x * (2 * l + 1) * p1 - (l + m) * p2) / (l - m + 1)
            rhon *= invR * (l - m + 1)
        end
        pl = -pl * fact * y
        fact += 2
        eim *= ei
    end
    return harmonics
end

@inline function B2M!(branch::SingleBranch, system, harmonics, expansion_order)
    B2M!(system, branch, branch.bodies_index, harmonics, expansion_order)
end

function B2M!(branch::MultiBranch, systems, harmonics, expansion_order)
    # iterate over systems
    for (i,system) in enumerate(systems)
        B2M!(system, branch, branch.bodies_index[i], harmonics, expansion_order)
    end
end

function M2B!(target_potential, target, i_branch, tree)
    branch = tree.branches[i_branch]
    T = eltype(branch.multipole_expansion[1])
    irregular_harmonics = Vector{eltype(branch.multipole_expansion[1])}(undef, (tree.expansion_order+1)^2)
    irregular_harmonics .= zero(eltype(irregular_harmonics))
    dx = target[1:3] - branch.center
    r, theta, phi = cartesian_2_spherical(dx)
    irregular_harmonics = irregular_harmonic!(irregular_harmonics, r, theta, phi, tree.expansion_order)
    d_potential = zeros(4)
    for l in 0:tree.expansion_order
        for m in 0:l
            ip = l^2 + l + m + 1
            i_compressed = 1 + (l * (l + 1)) >> 1 + m # only save half as Yl{-m} = conj(Ylm)
            factor = m > 0 ? 2.0 : 1.0
            for dim in 1:4
                d_potential[dim] += factor*real(branch.multipole_expansion[dim,i_compressed] * irregular_harmonics[ip])
            end
        end
    end
    target_potential .+= d_potential
end

function M2M!(branch, child, harmonics, M, expansion_order)
    # get distance vector
    dx, dy, dz = branch.center - child.center
    r, theta, phi = cartesian_2_spherical(dx, dy, dz)
    harmonics .= regular_harmonic!(harmonics, r, theta, phi, expansion_order)
    branch.multipole_expansion .= M2M_loop!(branch.multipole_expansion,child.multipole_expansion,harmonics,expansion_order)
end

function M2M_loop!(BM,CM,h,P)

    M = zeros(eltype(BM), 4)
    for j in 0:P # iterate over new multipole coefficients B_j^k
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
                    for dim in 1:4
                        @inbounds M[dim] += CM[dim,jlkms] * h[lm] * ipow * oddeven
                    end
                end
                for m in k:min(l,j+k-l)
                    jlkms = (((j-l) * (j-l+1)) >> 1) - k + m + 1
                    lm = l * l + l - m + 1
                    oddeven = odd_or_even(k + l + m)
                    C_tmp = h[lm] * oddeven
                    for dim in 1:4
                        @inbounds M[dim] += conj(CM[dim,jlkms]) * C_tmp
                    end
                end
            end
            for dim in 1:4
                @inbounds BM[dim,i_jk] += M[dim]
            end
        end
    end
    return BM
end

function M2L_loop!(local_expansion, L, multipole_expansion, harmonics, expansion_order)
    for j in 0:expansion_order
        Cnm = odd_or_even(j)
        for k in 0:j
            jks = (j * (j + 1)) >> 1 + k + 1
            L .= zero(eltype(L))
            for n in 0:expansion_order
                for m in -n:-1
                    nms = (n * (n+1)) >> 1 - m + 1
                    jnkm = (j + n)^2 + j + n + m - k + 1
                    # jnkm_max = (P + P)^2 + P + P + -1 - 0 + 1 = (2P)^2 + 2P = 2P(2P+1)
                    Cnm_tmp = Cnm * harmonics[jnkm]
                    for dim in 1:4
                        @inbounds L[dim] += conj(multipole_expansion[dim,nms]) * Cnm_tmp
                    end
                end
                for m in 0:n
                    nms = (n * (n+1)) >> 1 + m + 1
                    jnkm = (j + n) * (j + n) + j + n + m - k + 1
                    # jnkm_max = 2P * 2P + 2P + P + P - 0 + 1 = (2P)^2 + 2P + 2P + 1 = 4P^2 + 4P + 1 = (2P + 1)^2
                    Cnm_tmp = Cnm * odd_or_even((k-m) * (1 >> (k>=m)) + m) * harmonics[jnkm]
                    for dim in 1:4
                        @inbounds L[dim] += multipole_expansion[dim,nms] * Cnm_tmp
                    end
                end
            end
            view(local_expansion,:,jks) .+= L
        end
    end
    return local_expansion
end

function M2L!(target_branch, source_branch, harmonics, L, expansion_order)
    twice_expansion_order = expansion_order << 1
    #l = length(ReverseDiff.tape(target_branch.center))
    #harmonics = Vector{eltype(target_branch.multipole_expansion)}(undef, (twice_expansion_order + 1)*(twice_expansion_order + 1))
    harmonics = zeros(eltype(target_branch.multipole_expansion),(twice_expansion_order + 1)*(twice_expansion_order + 1))
    dx, dy, dz = target_branch.center - source_branch.center
    #@show length(ReverseDiff.tape(target_branch.center)) - l
    #l = length(ReverseDiff.tape(target_branch.center))
    r, theta, phi = cartesian_2_spherical(dx, dy, dz)
    irregular_harmonic!(harmonics, r, theta, phi, twice_expansion_order)
    M2L_loop!(target_branch.local_expansion, L, source_branch.multipole_expansion, harmonics, expansion_order)
end

function M2L!(target_branch, source_branch, expansion_order)
    twice_expansion_order = expansion_order << 1
    dx, dy, dz = target_branch.center - source_branch.center
    r, theta, phi = cartesian_2_spherical(dx, dy, dz)
    irregular_harmonic!(source_branch.harmonics, r, theta, phi, twice_expansion_order)
    M2L_loop!(target_branch.local_expansion, source_branch.ML, source_branch.multipole_expansion, source_branch.harmonics, expansion_order)
end

function B2L!(tree, i_branch, source_position, source_strength)
    branch = tree.branches[i_branch]
    irregular_harmonics = Vector{eltype(branch.multipole_expansion)}(undef, (tree.expansion_order+1)^2)
    r, theta, phi = cartesian_2_spherical(source_position - branch.center)
    irregular_harmonic!(irregular_harmonics, r, theta, -phi, tree.expansion_order)
    for l in 0:tree.expansion_order
        for m in 0:l
            i_abb = (l * (l+1)) >> 1 + m + 1
            i_exp = l^2 + l + m + 1
            for dim in 1:4
                branch.local_expansion[dim,i_abb] = irregular_harmonics[i_exp] * source_strength[dim]
            end
        end
    end
end

function L2L!(branch, child, harmonics, L, expansion_order) # 15-50k tape entries from this, also only in the loop.
    dx, dy, dz = child.center - branch.center
    r, theta, phi = cartesian_2_spherical(dx, dy, dz)
    r0 = 0.0*r + 0.0*im*r
    if !(ReverseDiff.hastape(child.local_expansion))
        child.local_expansion .+= r0
    end
    harmonics = regular_harmonic!(harmonics, r, theta, phi, expansion_order)
    #child.local_expansion .= L2L_loop!(child.local_expansion,branch.local_expansion,harmonics,L,expansion_order)
    #@show sum(child.local_expansion)
    #child.local_expansion = L2L_loop!(child.local_expansion,branch.local_expansion,harmonics,L,expansion_order)
    L2L_loop!(child.local_expansion,branch.local_expansion,harmonics,L,expansion_order)
    #@show sum(child.local_expansion)
    return nothing
end
@inline function L2L_loop!(CLE,BLE,h,L,P)
    for j in 0:P
        for k in 0:j
            jks = (j * (j + 1)) >> 1 + k + 1
            L .= zero(eltype(L))
            for n in j:P
                for m in j+k-n:-1
                    jnkm = (n-j) * (n-j) + n - j + m - k + 1
                    nms = (n * (n + 1)) >> 1 - m + 1
                    oddeven = odd_or_even(k)
                    #C_tmp = regular_harmonics[jnkm] * oddeven
                    for dim in 1:4
                        @inbounds L[dim] += conj(BLE[dim,nms]) * h[jnkm] * oddeven
                    end
                end
                for m in 0:n
                    if n-j >= abs(m-k)
                        jnkm = (n - j) * (n - j) + n - j + m - k + 1
                        nms = (n * (n + 1)) >> 1 + m + 1
                        oddeven = odd_or_even((m-k) * (1 >> (m >= k)))
                        C_tmp = h[jnkm] * oddeven
                        L .+= view(BLE,:,nms) .* C_tmp
                    end
                end
            end
            CLE[:,jks] .+= L
        end
    end
    return CLE
end

function L2B!(systems, branch::MultiBranch, expansion_order, vector_potential, potential_jacobian, potential_hessian, harmonics, workspace)
    for (i_system, system) in enumerate(systems)
        L2B!(system, branch.bodies_index[i_system], branch.local_expansion, expansion_order, branch.center, vector_potential, potential_jacobian, potential_hessian, harmonics, workspace)
    end
end

function L2B!(system, branch::SingleBranch, expansion_order, vector_potential, potential_jacobian, potential_hessian, harmonics, workspace)
    L2B!(system, branch.bodies_index, branch.local_expansion, expansion_order, branch.center, vector_potential, potential_jacobian, potential_hessian, harmonics, workspace)
end

function L2B!(system, bodies_index, local_expansion, expansion_order, expansion_center, vector_potential, potential_jacobian, potential_hessian, harmonics, workspace)
    for i_body in bodies_index
        #l = length(ReverseDiff.tape(system.bodies[1].position[2]))
        vector_potential = zeros(eltype(vector_potential), size(vector_potential))
        potential_jacobian = zeros(eltype(potential_jacobian), size(potential_jacobian))
        potential_hessian = zeros(eltype(potential_hessian), size(potential_hessian))
        #vector_potential .= zero(eltype(vector_potential))
        #potential_jacobian .= zero(eltype(potential_jacobian))
        #potential_hessian .= zero(eltype(potential_hessian))
        body_position = system[i_body,POSITION]
        #println("tape entries before L2B loop: $(length(ReverseDiff.tape(system.bodies[1].position[2])) - l)")
        #l = length(ReverseDiff.tape(system.bodies[1].position[2]))
        scalar_potential, vector_potential, potential_jacobian, potential_hessian = L2B_loop!(vector_potential, potential_jacobian, potential_hessian, body_position, expansion_center, local_expansion, harmonics, expansion_order, workspace)
        #println("tape entries from L2B loop: $(length(ReverseDiff.tape(system.bodies[1].position[2])) - l)")
        #l = length(ReverseDiff.tape(system.bodies[1].position[2]))
        system[i_body,SCALAR_POTENTIAL] += scalar_potential
        # note: system[i,VECTOR_POTENTIAL], system[i,VELOCITY], and system[i,VELOCITY_GRADIENT] must be mutable
        system[i_body,VECTOR_POTENTIAL] .+= vector_potential
        system[i_body,VELOCITY] .+= view(potential_jacobian,:,1)
        system[i_body,VELOCITY_GRADIENT] .+= view(potential_hessian,:,:,1)
        #println("tape entries after L2B loop: $(length(ReverseDiff.tape(system.bodies[1].position[2])) - l)")
    end
end

function L2B_loop!(vector_potential, potential_jacobian, potential_hessian, body_position, expansion_center, local_expansion, harmonics, expansion_order, workspace)
    #l = length(ReverseDiff.tape(body_position[2]))
    dx, dy, dz = body_position - expansion_center
    #println("tape entries from first line: $(length(ReverseDiff.tape(body_position[2])) - l)")
    #l = length(ReverseDiff.tape(body_position[2]))
    r, theta, phi = cartesian_2_spherical(dx, dy, dz)
    #println("tape entries from cartesian_2_spherical: $(length(ReverseDiff.tape(body_position[2])) - l)")
    #l = length(ReverseDiff.tape(body_position[2]))
    #harmonics = regular_harmonic_derivs!(harmonics, r, theta, phi, expansion_order)
    #@show typeof(harmonics)
    #harmonics .= regular_harmonic_derivs!(harmonics, r, theta, phi, expansion_order)
    harmonics = regular_harmonic_derivs!(harmonics, r, theta, phi, expansion_order) # 447423 tape entries vs. 368273 tape entries vs. 487023 tape entries
    #@show sum(harmonics)
    #@show typeof(harmonics)
    #println("tape entries from regular_harmonic_derivs!: $(length(ReverseDiff.tape(body_position[2])) - l)")
    #l = length(ReverseDiff.tape(body_position[2]))
    scalar_potential = update_scalar_potential(0.0,local_expansion,harmonics[1,:],expansion_order)
    #println("tape entries from scalar update: $(length(ReverseDiff.tape(body_position[2])) - l)")
    #l = length(ReverseDiff.tape(body_position[2]))
    #vector_potential .= update_vector_potential!(vector_potential,local_expansion,harmonics[1,:],expansion_order)
    vector_potential = update_vector_potential!(vector_potential,local_expansion,harmonics[1,:],expansion_order)
    #update_vector_potential!(vector_potential,local_expansion,harmonics[1,:],expansion_order)
    #println("tape entries from vector update: $(length(ReverseDiff.tape(body_position[2])) - l)")
    #l = length(ReverseDiff.tape(body_position[2]))
    potential_jacobian = update_potential_jacobian!(potential_jacobian,local_expansion,harmonics[1,:],harmonics[2,:],expansion_order,r)
    #println("tape entries from jacobian update: $(length(ReverseDiff.tape(body_position[2])) - l)")
    #l = length(ReverseDiff.tape(body_position[2]))
    potential_hessian = update_potential_hessian!(potential_hessian,local_expansion,harmonics[1,:],harmonics[2,:],harmonics[3,:],expansion_order,r)
    #println("tape entries from hessian updates: $(length(ReverseDiff.tape(body_position[2])) - l)")
    #l = length(ReverseDiff.tape(body_position[2]))
    #@show size(potential_hessian) size(local_expansion) size(harmonics) size(harmonics_theta) size(harmonics_theta_2)
    #spherical_2_cartesian!(potential_jacobian, potential_hessian, workspace, r, theta, phi)
    potential_hessian = s2c_hess!(potential_jacobian, potential_hessian, workspace, r, theta, phi)
    potential_jacobian = s2c_jac!(potential_jacobian, workspace, r, theta, phi)
    #println("tape entries from potential s2c: $(length(ReverseDiff.tape(body_position[2])) - l)")
    #l = length(ReverseDiff.tape(body_position[2]))
    potential_jacobian = flatten_jacobian!(potential_jacobian)
    potential_hessian = flatten_hessian!(potential_hessian)
    #non_elementify(potential_hessian,flatten_hessian!,potential_hessian)
    #println("tape entries from flattening jacobian and hessian: $(length(ReverseDiff.tape(body_position[2])) - l)")
    #l = length(ReverseDiff.tape(body_position[2]))
    #flatten_derivatives!(potential_jacobian, potential_hessian) # compute velocity and velocity gradient
    return scalar_potential, vector_potential, potential_jacobian, potential_hessian
end
