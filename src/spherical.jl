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
function get_drjdxi(s_theta,s_phi,c_theta,c_phi,rho)
    return @SMatrix [s_theta*c_phi    c_theta*c_phi/rho      -s_phi/rho/s_theta
            s_theta * s_phi  c_theta * s_phi / rho  c_phi / rho / s_theta
            c_theta         -s_theta / rho          0                    ]
end

# replaced repeated sin/cos calls with precomputed values. Reduces ReverseDiff allocations due to fewer functions recorded to the tape. Also increases readability.
function get_drkdxidxj(s_theta,s_phi,c_theta,c_phi,rho,k_coord)

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

    drjdxi = get_drjdxi(s_theta,s_phi,c_theta,c_phi,rho)

    workspace .= potential_jacobian

    mul!(potential_jacobian, drjdxi, workspace)

    return potential_jacobian

end

function s2c_hess!(potential_jacobian, potential_hessian, workspace, rho, theta, phi)

    s_theta, c_theta = sincos(theta)
    s_phi, c_phi = sincos(phi)

    drjdxi = get_drjdxi(s_theta,s_phi,c_theta,c_phi,rho)
    # we only need to evaluate these functions once instead of 4 times
    drkdxidxj_r = get_drkdxidxj(s_theta,s_phi,c_theta,c_phi,rho,1)
    drkdxidxj_theta = get_drkdxidxj(s_theta,s_phi,c_theta,c_phi,rho,2)
    drkdxidxj_phi = get_drkdxidxj(s_theta,s_phi,c_theta,c_phi,rho,3)

    # convert Hessian to cartesian coordinates
    workspace3x3 = view(workspace,:,1:3)
    for ind in 1:4
        workspace3x3 .= potential_hessian[:,:,ind]
        potential_hessian[:,:,ind] .= drjdxi * workspace3x3 * transpose(drjdxi)
        potential_hessian[:,:,ind] .+= drkdxidxj_r * potential_jacobian[1,ind]
        potential_hessian[:,:,ind] .+= drkdxidxj_theta * potential_jacobian[2,ind]
        potential_hessian[:,:,ind] .+= drkdxidxj_phi * potential_jacobian[3,ind]
    end
    return potential_hessian

end

# I split these so that I could add custom reverse diff rules to each part separately. Doing it with one function is an enormous headache; it should be possible but there is no documentation on how to handle functions that return tuples.
function spherical_2_cartesian!(potential_jacobian, potential_hessian, workspace, rho, theta, phi)
    # the order of these two matters, since the hessian needs to use the original jacobian.
    potential_hessian = s2c_hess!(potential_jacobian, potential_hessian, workspace, rho, theta, phi)
    potential_jacobian = s2c_jac!(potential_jacobian, workspace, rho, theta, phi)

    return potential_jacobian,potential_hessian
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
    jacobian = flatten_jacobian!(jacobian)
    # velocity gradient
    hessian = flatten_hessian!(hessian)
end

@inline odd_or_even(n::Int) = (n & 1) == 1 ? -1 : 1

@inline ipow2l(n::Int) = n >= 0 ? 1 : odd_or_even(n);

# assume harmonics (and theta derivatives) are real nx2 matrices rather than complex nx1 vectors
function regular_harmonic!(harmonics, harmonics_theta, harmonics_theta_2, rho, theta, phi, P)
    y,x = sincos(theta)
    invY = y == 0 ? 0 : 1 / y
    fact = 1.0
    pl = 1.0
    rhom = 1.0
    #ei = exp(im * phi) # split real and complex parts
    ei = [cos(phi) sin(phi)]
    #eim = 1.0 # split real and complex parts
    eim = eltype(harmonics)[1.0 0.0]
    for m=0:P
        p = pl
        lpl = (m * (m + 1)) >> 1 + m + 1
        # lpl = m * m + 2 * m + 1
        # lml = m * m + 1
        #harmonics[lpl] = rhom * p * eim # split real and complex parts; do complex multiplication manually
        harmonics[lpl,1] = p * rhom * eim[1]
        harmonics[lpl,2] = p * rhom * eim[2]
        # harmonics[lml] = conj(harmonics[lpl])
        p1 = p
        p = x * (2 * m + 1) * p1
        #harmonics_theta[lpl] = rhom * (p - (m + 1) * x * p1) * invY * eim # more manual complex multiplication
        harmonics_theta[lpl,1] = (p - (m + 1) * x * p1) * invY * rhom * eim[1]
        harmonics_theta[lpl,2] = (p - (m + 1) * x * p1) * invY * rhom * eim[2]
        #harmonics_theta_2[lpl] = rhom * (-x * p + (-m + (m+1)^2 * x^2) * p1) * invY^2 * eim # more manual complex multiplication
        harmonics_theta_2[lpl,1] = rhom * (-x * p + (-m + (m+1)^2 * x^2) * p1) * invY^2 * eim[1]
        harmonics_theta_2[lpl,2] = rhom * (-x * p + (-m + (m+1)^2 * x^2) * p1) * invY^2 * eim[2]

        rhom *= rho
        rhol = rhom
        for l=m+1:P
            lpm = (l * (l + 1)) >> 1 + m + 1
            # lpm = l * l + l + m + 1
            # lmm = l * l + l - m + 1
            rhol /= -(l + m)
            #harmonics[lpm] = rhol * p * eim # more manual complex multiplication
            harmonics[lpm,1] = rhol * p * eim[1]
            harmonics[lpm,2] = rhol * p * eim[2]
            # harmonics[lmm] = conj(harmonics[lpm])
            p2 = p1
            p1 = p
            p = (x * (2 * l + 1) * p1 - (l + m) * p2) / (l - m + 1)
            #harmonics_theta[lpm] = rhol * ((l - m + 1) * p - (l + 1) * x * p1) * invY * eim # more manual complex multiplication
            harmonics_theta[lpm,1] = rhol * ((l - m + 1) * p - (l + 1) * x * p1) * invY * eim[1]
            harmonics_theta[lpm,2] = rhol * ((l - m + 1) * p - (l + 1) * x * p1) * invY * eim[2]
            #harmonics_theta_2[lpm] = rhol * ((m-l-1) * x * p + (m^2 - l*(l+1) + (l+1)^2 * x^2) * p1) * invY^2 * eim # more manual complex multiplication
            harmonics_theta_2[lpm,1] = rhol * ((m-l-1) * x * p + (m^2 - l*(l+1) + (l+1)^2 * x^2) * p1) * invY^2 * eim[1]
            harmonics_theta_2[lpm,2] = rhol * ((m-l-1) * x * p + (m^2 - l*(l+1) + (l+1)^2 * x^2) * p1) * invY^2 * eim[2]
            rhol *= rho
        end
        rhom /= -(2 * m + 2) * (2 * m + 1)
        pl = -pl * fact * y
        fact += 2
        #eim *= ei # more manual complex multiplication
        eim[1] = eim[1]*ei[1] - eim[2]*ei[2]
        eim[2] = eim[1]*ei[2] + eim[2]*ei[1]

    end
end

# same real-to-complex changes as previous function.
function regular_harmonic!(harmonics, rho, theta, phi, P)
    y,x = sincos(theta)
    fact = 1.0
    pl = 1.0
    rhom = 1.0 # rho^l / (l+m)! * (-1)^l
    #ei = exp(im * phi)
    ei = [cos(phi) sin(phi)]
    #eim = 1.0
    eim = eltype(harmonics)[1.0 0.0] # e^(i * m * phi)
    for m=0:P # l=m up here
        p = pl
        lpl = m * m + 2 * m + 1
        lml = m * m + 1
        #harmonics[lpl] = rhom * p * eim
        harmonics[lpl,1] = rhom * p * eim[1]
        harmonics[lpl,2] = rhom * p * eim[2]
        #harmonics[lml] = conj(harmonics[lpl])
        harmonics[lml,1] = harmonics[lpl,1]
        harmonics[lml,2] = -harmonics[lpl,2]
        p1 = p
        p = x * (2 * m + 1) * p1
        rhom *= rho
        rhol = rhom
        for l=m+1:P # l>m in here
            lpm = l * l + l + m + 1
            lmm = l * l + l - m + 1
            rhol /= -(l + m)
            #harmonics[lpm] = rhol * p * eim
            harmonics[lpm,1] = rhol * p * eim[1]
            harmonics[lpm,2] = rhol * p * eim[2]
            #harmonics[lmm] = conj(harmonics[lpm])
            harmonics[lmm,1] = harmonics[lpm,1]
            harmonics[lmm,2] = -harmonics[lpm,2]
            p2 = p1
            p1 = p
            p = (x * (2 * l + 1) * p1 - (l + m) * p2) / (l - m + 1)
            rhol *= rho
        end
        rhom /= -(2 * m + 2) * (2 * m + 1)
        pl = -pl * fact * y
        fact += 2
        #eim *= ei
        eim[1] = eim[1]*ei[1] - eim[2]*ei[2]
        eim[2] = eim[1]*ei[2] + eim[2]*ei[1]
    end
    return harmonics
end

# same changes as before
function irregular_harmonic!(harmonics, rho, theta, phi, P)
    y, x = sincos(theta)
    fact = 1
    pl = 1
    invR = -1.0 / rho
    rhom = -invR
    #ei = exp(im * phi)
    ei = [cos(phi) sin[phi]]
    eim = eltype(harmonics)[1.0 0.0]
    for m=0:P
        p = pl
        npl = m * m + 2 * m + 1
        nml = m * m + 1
        #harmonics[npl] = rhom * p * eim
        harmonics[npl,1] = rhom * p * eim[1]
        harmonics[npl,2] = rhom * p * eim[2]
        #harmonics[nml] = conj(harmonics[npl])
        harmonics[nml,1] = harmonics[npl,1]
        harmonics[nml,2] = -harmonics[npl,2]
        p1 = p
        p = x * (2 * m + 1) * p1
        rhom *= invR
        rhon = rhom
        for l=m+1:P
            npm = l * l + l + m + 1
            nmm = l * l + l - m + 1
            # npm_max = P^2 + P + P + 1 = (P+1)^2
            #harmonics[npm] = rhon * p * eim
            harmonics[npm,1] = rhon * p * eim[1]
            harmonics[npm,2] = rhon * p * eim[2]
            #harmonics[nmm] = conj(harmonics[npm])
            harmonics[nmm,1] = harmonics[npm,1]
            harmonics[nmm,2] = -harmonics[npm,2]
            p2 = p1
            p1 = p
            p = (x * (2 * l + 1) * p1 - (l + m) * p2) / (l - m + 1)
            rhon *= invR * (l - m + 1)
        end
        pl = -pl * fact * y
        fact += 2
        #eim *= ei
        eim[1] = eim[1]*ei[1] - eim[2]*ei[2]
        eim[2] = eim[1]*ei[2] + eim[2]*ei[1]
    end
    return harmonics
end

function B2M!(tree, systems, i_branch, sources_index)
    
    branch = tree.branches[i_branch]
    T = eltype(branch.multipole_expansion[1])
    # harmonics = Vector{eltype(branch.multipole_expansion[1])}(undef, (tree.expansion_order+1)^2) # this is faster than MVector
    harmonics = Matrix{T}(undef, (tree.expansion_order+1)^2,2)
    harmonics .= zero(T)

    # iterate over elements
    for (i_iter, system) in enumerate(systems[sources_index])
        i_type = sources_index[i_iter]
        bodies_index = branch.first_body[i_type]:branch.first_body[i_type] + branch.n_bodies[i_type] - 1
        B2M!(branch, system, bodies_index, harmonics, tree.expansion_order)
    end
end

function M2B!(target_potential, target, i_branch, tree)
    branch = tree.branches[i_branch]
    T = eltype(branch.multipole_expansion[1])
    #irregular_harmonics = Vector{eltype(branch.multipole_expansion[1])}(undef, (tree.expansion_order+1)^2)
    irregular_harmonics = Array{T}(undef, (tree.expansion_order+1)^2,2) # split real and complex components
    irregular_harmonics .= zero(T)
    dx = target[1:3] - branch.center
    r, theta, phi = cartesian_2_spherical(dx)
    irregular_harmonics = irregular_harmonic!(irregular_harmonics, r, theta, phi, tree.expansion_order)
    d_potential = zeros(4)
    for l in 0:tree.expansion_order
        for m in 0:l
            ip = l^2 + l + m + 1
            i_compressed = 1 + (l * (l + 1)) >> 1 + m # only save half as Yl{-m} = conj(Ylm)
            for dim in 1:4
                d_potential[dim] = real(branch.multipole_expansion[dim][i_compressed] * irregular_harmonics[ip])
            end
            m > 0 && (d_potential .*= 2)
            target_potential .+= d_potential
        end
    end
end

function M2M!(tree, branch, child, harmonics) # 17.8k tape entries, all in the loop.
    # get distance vector
    dx, dy, dz = branch.center - child.center
    r, theta, phi = cartesian_2_spherical(dx, dy, dz)
    harmonics = regular_harmonic!(harmonics, r, theta, phi, tree.expansion_order)
    # this next line should be redundant
    #harmonics = regular_harmonic_real!(harmonics, r, theta, phi, tree.expansion_order) .+ regular_harmonic_imag!(harmonics, r, theta, phi, tree.expansion_order).*im

    #M = zeros(eltype(branch.multipole_expansion[1]), 4)
    M = zeros(eltype(branch.multipole_expansion[1]), 4, 2) # split real and complex components
    #@show length(branch.multipole_expansion[1][1,1].tape)
    for j in 0:tree.expansion_order # iterate over new Multipole coefficients B_j^k
        for k in 0:j
            i_jk = ((j * (j+1)) >> 1) + k + 1 # current index
            M .*= 0.0
            for l in 0:j
                for m in max(-l,-j+k+l):min(k-1,l)
                    jlkms = (((j-l) * (j-l+1)) >> 1) + k - m + 1
                    lm = l * l + l - m + 1
                    ipow = ipow2l(m)
                    oddeven = odd_or_even(l)
                    for dim in 1:4
                        #M[dim] += child.multipole_expansion[dim][jlkms] * harmonics[lm] * ipow * oddeven
                        M[dim,1] += (child.multipole_expansion[dim][jlkms,1] * harmonics[lm,1] - child.multipole_expansion[dim][jlkms,2] * harmonics[lm,2]) * ipow * oddeven
                        M[dim,2] += (child.multipole_expansion[dim][jlkms,1] * harmonics[lm,2] + child.multipole_expansion[dim][jlkms,2] * harmonics[lm,1]) * ipow * oddeven
                    end
                end
                for m in k:min(l,j+k-l)
                    jlkms = (((j-l) * (j-l+1)) >> 1) - k + m + 1
                    lm = l * l + l - m + 1
                    oddeven = odd_or_even(k + l + m)
                    for dim in 1:4
                        #M[dim] += conj(child.multipole_expansion[dim][jlkms]) * harmonics[lm] * oddeven
                        M[dim,1] += (child.multipole_expansion[dim][jlkms,1] * harmonics[lm,1] + child.multipole_expansion[dim][jlkms,2] * harmonics[lm,2]) * oddeven
                        M[dim,2] += (child.multipole_expansion[dim][jlkms,1] * harmonics[lm,2] - child.multipole_expansion[dim][jlkms,2] * harmonics[lm,1]) * oddeven
                    end
                end
            end
            for dim in 1:4
                #branch.multipole_expansion[dim][i_jk] += M[dim]
                branch.multipole_expansion[dim][i_jk,1] += M[dim,1]
                branch.multipole_expansion[dim][i_jk,2] += M[dim,2]
            end
        end
    end
    #@show length(branch.multipole_expansion[1][1,1].tape)
end

function M2M!(tree, i_branch)
    # expose objects
    branch = tree.branches[i_branch]
    T = eltype(branch.multipole_expansion[1])
    #initialize memory TODO: do this beforehand?
    harmonics = Matrix{T}(undef, (tree.expansion_order+1)^2,2)
    harmonics .= zero(T)

    # iterate over children
    for i_child in branch.first_branch:branch.first_branch + branch.n_branches - 1
        child = tree.branches[i_child]
        M2M!(tree, branch, child, harmonics)
    end
end

function M2L_loop!(local_expansion, L, multipole_expansion, harmonics, expansion_order)
    for j in 0:expansion_order
        Cnm = odd_or_even(j)
        for k in 0:j
            jks = (j * (j + 1)) >> 1 + k + 1
            L .*= 0.0
            for n in 0:expansion_order
                for m in -n:-1
                    nms = (n * (n+1)) >> 1 - m + 1
                    jnkm = (j + n)^2 + j + n + m - k + 1
                    # jnkm_max = (P + P)^2 + P + P + -1 - 0 + 1 = (2P)^2 + 2P = 2P(2P+1)
                    for dim in 1:4
                        #L[dim] += conj(multipole_expansion[dim][nms]) * Cnm * harmonics[jnkm]
                        L[dim,1] += (multipole_expansion[dim][nms,1] * harmonics[jnkm,1] + multipole_expansion[dim][nms,2] * harmonics[jnkm,2]) * Cnm
                        L[dim,2] += (multipole_expansion[dim][nms,1] * harmonics[jnkm,2] - multipole_expansion[dim][nms,2] * harmonics[jnkm,1]) * Cnm
                    end
                end
                for m in 0:n
                    nms = (n * (n+1)) >> 1 + m + 1
                    jnkm = (j + n) * (j + n) + j + n + m - k + 1
                    # jnkm_max = 2P * 2P + 2P + P + P - 0 + 1 = (2P)^2 + 2P + 2P + 1 = 4P^2 + 4P + 1 = (2P + 1)^2
                    Cnm2 = Cnm * odd_or_even((k-m) * (1 >> (k>=m)) + m)
                    for dim in 1:4
                        #L[dim] += multipole_expansion[dim][nms] * Cnm2 * harmonics[jnkm]
                        L[dim,1] += (multipole_expansion[dim][nms,1] * harmonics[jnkm,1] - multipole_expansion[dim][nms,2] * harmonics[jnkm,2]) * Cnm2
                        L[dim,2] += (multipole_expansion[dim][nms,1] * harmonics[jnkm,2] + multipole_expansion[dim][nms,2] * harmonics[jnkm,1]) * Cnm2
                    end
                end
            end
            for dim in 1:4
                #local_expansion[dim][jks] += L[dim]
                local_expansion[dim][jks,1] += L[dim,1]
                local_expansion[dim][jks,2] += L[dim,2]
            end
        end
    end
end

function M2L!(tree, i_local, j_multipole)
    local_branch = tree.branches[i_local]
    multipole_branch = tree.branches[j_multipole]
    # preallocate
    T = eltype(local_branch.multipole_expansion[1])
    #harmonics = Vector{eltype(local_branch.multipole_expansion[1])}(undef, (2*tree.expansion_order + 1)^2)
    harmonics = Matrix{T}(undef, (2*tree.expansion_order + 1)^2, 2)
    harmonics .= zero(T)
    # get separation vector
    dx, dy, dz = local_branch.center - multipole_branch.center
    r, theta, phi = cartesian_2_spherical(dx, dy, dz)
    irregular_harmonic!(harmonics, r, theta, phi, 2*tree.expansion_order)
    #L = zeros(eltype(local_branch.local_expansion[1]), 4)
    L = zeros(eltype(local_branch.local_expansion[1]), 4, 2)
    M2L_loop!(local_branch.local_expansion, L, multipole_branch.multipole_expansion, harmonics, tree.expansion_order)
end

# I'm assuming that source_strength is complex but I'm not sure
function B2L!(tree, i_branch, source_position, source_strength)
    branch = tree.branches[i_branch]
    T = eltype(branch.multipole_expansion[1])
    #irregular_harmonics = Vector{eltype(branch.multipole_expansion[1])}(undef, (tree.expansion_order+1)^2)
    irregular_harmonics = Matrix{T}(undef, (tree.expansion_order+1)^2, 2)
    irregular_harmonics .= zero(T)
    r, theta, phi = cartesian_2_spherical(source_position - branch.center)
    irregular_harmonic!(irregular_harmonics, r, theta, -phi, tree.expansion_order)
    for l in 0:tree.expansion_order
        for m in 0:l
            i_abb = (l * (l+1)) >> 1 + m + 1
            i_exp = l^2 + l + m + 1
            for dim in 1:4
                #branch.local_expansion[dim][i_abb] = irregular_harmonics[i_exp] * source_strength[dim]
                branch.local_expansion[dim][i_abb,1] = irregular_harmonics[i_exp,1] * source_strength[dim,1] - irregular_harmonics[i_exp,2] * source_strength[dim,2]
                branch.local_expansion[dim][i_abb,2] = irregular_harmonics[i_exp,1] * source_strength[dim,2] + irregular_harmonics[i_exp,2] * source_strength[dim,1]
            end
        end
    end
end

function L2L!(tree, branch, child, harmonics) # 15-50k tape entries from this, also only in the loop.
    dx, dy, dz = child.center - branch.center
    r, theta, phi = cartesian_2_spherical(dx, dy, dz)
    #@show length(child.local_expansion[1][1,1].tape)
    harmonics = regular_harmonic!(harmonics, r, theta, phi, tree.expansion_order)
    #@show length(child.local_expansion[1][1,1].tape)
    #regular_harmonic_mod!(harmonics, r, theta, phi, tree.expansion_order)
    for j in 0:tree.expansion_order
        for k in 0:j
            jks = (j * (j + 1)) >> 1 + k + 1
            #L = zeros(eltype(branch.local_expansion[1]), 4)
            L = zeros(eltype(branch.local_expansion[1]), 4, 2)
            for n in j:tree.expansion_order
                for m in j+k-n:-1
                    jnkm = (n-j) * (n-j) + n - j + m - k + 1
                    nms = (n * (n + 1)) >> 1 - m + 1
                    oddeven = odd_or_even(k)
                    for dim in 1:4
                        #L[dim] += conj(branch.local_expansion[dim][nms]) * harmonics[jnkm] * oddeven
                        L[dim,1] += (branch.local_expansion[dim][nms,1] * harmonics[jnkm,1] + branch.local_expansion[dim][nms,2] * harmonics[jnkm,2]) * oddeven
                        L[dim,2] += (branch.local_expansion[dim][nms,1] * harmonics[jnkm,2] - branch.local_expansion[dim][nms,2] * harmonics[jnkm,1]) * oddeven
                    end
                end
                for m in 0:n
                    if n-j >= abs(m-k)
                        jnkm = (n - j) * (n - j) + n - j + m - k + 1
                        nms = (n * (n + 1)) >> 1 + m + 1
                        oddeven = odd_or_even((m-k) * (1 >> (m >= k)))
                        for dim in 1:4
                            #L[dim] += branch.local_expansion[dim][nms] * harmonics[jnkm] * oddeven
                            L[dim,1] += (branch.local_expansion[dim][nms,1] * harmonics[jnkm,1] - branch.local_expansion[dim][nms,2] * harmonics[jnkm,2]) * oddeven
                            L[dim,2] += (branch.local_expansion[dim][nms,1] * harmonics[jnkm,2] + branch.local_expansion[dim][nms,2] * harmonics[jnkm,1]) * oddeven
                        end
                    end
                end
            end
            for dim in 1:4
                #child.local_expansion[dim][jks] += L[dim]
                child.local_expansion[dim][jks,1] += L[dim,1]
                child.local_expansion[dim][jks,2] += L[dim,2]
            end
        end
    end
    #@show length(child.local_expansion[1][1,1].tape)
end

function L2L!(tree, j_source)
    # expose branch
    branch = tree.branches[j_source]

    #initialize memory TODO: do this beforehand?
    T = eltype(branch.multipole_expansion[1])
    #harmonics = Vector{eltype(branch.multipole_expansion[1])}(undef, (tree.expansion_order+1)^2)
    harmonics = Matrix{T}(undef, (tree.expansion_order+1)^2, 2)
    harmonics .= zero(T)

    # iterate over children
    for i_child in branch.first_branch:branch.first_branch + branch.n_branches - 1
        child = tree.branches[i_child]
        L2L!(tree, branch, child, harmonics)
    end
end

"Calculates the potential at all child elements of a branch."
function L2B!(tree, systems, i_branch, targets_index)
    branch = tree.branches[i_branch]
    T = eltype(branch.multipole_expansion[1])
    #harmonics = Vector{eltype(branch.multipole_expansion[1])}(undef, ((tree.expansion_order+1) * (tree.expansion_order+2)) >> 1)
    harmonics = zeros(T, ((tree.expansion_order+1) * (tree.expansion_order+2)) >> 1, 2)
    #harmonics_theta = zeros(eltype(branch.multipole_expansion[1]), ((tree.expansion_order+1) * (tree.expansion_order+2)) >> 1)
    harmonics_theta = zeros(T, ((tree.expansion_order+1) * (tree.expansion_order+2)) >> 1, 2)
    #harmonics_theta_2 = zeros(eltype(branch.multipole_expansion[1]), ((tree.expansion_order+1) * (tree.expansion_order+2)) >> 1)
    harmonics_theta_2 = zeros(T, ((tree.expansion_order+1) * (tree.expansion_order+2)) >> 1, 2)
    workspace = zeros(eltype(tree),3,4) # this is unused?
    spherical_potential = zeros(eltype(tree),52) # purely real?
    for (i_target, system) in enumerate(systems[targets_index])
        i_type = targets_index[i_target]
        for i_body in branch.first_body[i_type]:branch.first_body[i_type] + branch.n_bodies[i_type] - 1
            L2B!(system, i_body, harmonics, harmonics_theta, harmonics_theta_2, workspace, spherical_potential, tree, branch)
            spherical_potential .*= 0
        end
    end
end

function update_potential!(potential,LE,h,P)
    @show typeof(potential) typeof(LE) typeof(h) typeof(P)
    #@show size(LE)
    #@show typeof(LE) typeof(reduce(hcat,LE))
    for n in 0:P
        nms = (n * (n+1)) >> 1 + 1
        for ind in 1:4
            potential[ind] += LE[ind][nms,1] * h[nms,1] - LE[ind][nms,2] * h[nms,2]
        end
        for m in 1:n
            nms = (n * (n + 1)) >> 1 + m + 1
            for ind in 1:4
                potential[ind] += 2 * (LE[ind][nms,1] * h[nms,1] - LE[ind][nms,2] * h[nms,2])
            end
        end
    end
    return potential
end

function update_potential_jacobian!(potential_jacobian,LE,h,ht,r,P)
    for n in 0:P
        nms = (n * (n+1)) >> 1 + 1
        for ind in 1:4
            potential_jacobian[1,ind] += n/r * (LE[ind][nms,1] * h[nms,1] - LE[ind][nms,2] * h[nms,2])
            potential_jacobian[2,ind] += LE[ind][nms,1] * ht[nms,1] - LE[ind][nms,2] * ht[nms,2]
        end
        for m in 1:n
            nms = (n * (n + 1)) >> 1 + m + 1
            for ind in 1:4
                potential_jacobian[1,ind] += 2 * n/r * (LE[ind][nms,1] * h[nms,1] - LE[ind][nms,2] * h[nms,2]) # dPsi/dr
                potential_jacobian[2,ind] += 2 * (LE[ind][nms,1] * ht[nms,1] - LE[ind][nms,2] * ht[nms,2]) # dPsi/dtheta
                potential_jacobian[3,ind] += 2 * m * (LE[ind][nms,1] * h[nms,2] + LE[ind][nms,2] * h[nms,1]) # dPsi/dphi
            end
        end
    end
    return potential_jacobian
end

function update_potential_hessian!(potential_hessian,LE,h,ht,ht2,r,P)

    for n in 0:P
        nms = (n * (n+1)) >> 1 + 1
        for ind in 1:4
            potential_hessian[1,1,ind] += n * (n-1) / r^2 * (LE[ind][nms,1] * h[nms,1] - LE[ind][nms,2] * h[nms,2])
            potential_hessian[2,1,ind] += n/r * (LE[ind][nms,1] * ht[nms,1] - LE[ind][nms,2] * ht[nms,2])
            potential_hessian[1,2,ind] += n/r * (LE[ind][nms,1] * ht[nms,1] - LE[ind][nms,2] * ht[nms,2])
            potential_hessian[2,2,ind] += LE[ind][nms,1] * ht2[nms,1] - LE[ind][nms,2] * ht2[nms,2]
        end
        for m in 1:n
            nms = (n * (n + 1)) >> 1 + m + 1
            for ind in 1:4
                potential_hessian[1,1,ind] += 2 * n * (n-1) / r^2 * (LE[ind][nms,1] * h[nms,1] - LE[ind][nms,2] * h[nms,2]) # d2Psi/dr2
                potential_hessian[2,1,ind] += 2 * n/r * (LE[ind][nms,1] * ht[nms,1] - LE[ind][nms,2] * ht[nms,2]) # d2Psi/dtheta dr
                potential_hessian[3,1,ind] += 2 * n * m / r * (LE[ind][nms,1] * h[nms,2] + LE[ind][nms,2] * h[nms,1]) # d2Psi/dphi dr
                potential_hessian[1,2,ind] += 2 * n/r * (LE[ind][nms,1] * ht[nms,1] - LE[ind][nms,2] * ht[nms,2]) # d2Psi/dr dtheta
                potential_hessian[2,2,ind] += 2 * (LE[ind][nms,1] * ht2[nms,1] - LE[ind][nms,2] * ht2[nms,2]) # d2Psi/dtheta2
                potential_hessian[3,2,ind] += 2 * m * (LE[ind][nms,1] * ht[nms,2] + LE[ind][nms,2] * ht[nms,1]) # d2Psi/dphi dtheta
                potential_hessian[1,3,ind] += 2 * n * m / r * (LE[ind][nms,1] * h[nms,2] + LE[ind][nms,2] * h[nms,1]) # d2Psi/dr dphi
                potential_hessian[2,3,ind] += 2 * m * (LE[ind][nms,1] * ht[nms,2] + LE[ind][nms,2] * ht[nms,1]) # d2Psi/dtheta dphi
                potential_hessian[3,3,ind] += 2 * -m^2 * (LE[ind][nms,1] * h[nms,1] - LE[ind][nms,2] * h[nms,2]) # d2Psi/dphi2
            end
        end
    end
    return potential_hessian

end

@inline function L2B!(system, i_body, harmonics, harmonics_theta, harmonics_theta_2, workspace, spherical_potential, tree, branch) # contains tape entries ~70k to the end. i.e. 110k tape entries.
    println("starting L2B!")
    @show length(system[i_body,POTENTIAL][1].tape)
    potential = view(spherical_potential,1:4)
    potential_jacobian = reshape(view(spherical_potential, 5:16),3,4)
    potential_hessian = reshape(view(spherical_potential, 17:52),3,3,4)
    #local_expansion = cat(branch.local_expansion[1:4]...;dims=3) # this definitely allocates... there should be a way to use a view to do this better.
    #@show size(local_expansion) typeof(local_expansion)
    body_position = system[i_body,POSITION]
    dx, dy, dz = body_position - branch.center
    r, theta, phi = cartesian_2_spherical(dx, dy, dz)
    @show length(system[i_body,POTENTIAL][1].tape) # ~ + 12 tape entries
    regular_harmonic!(harmonics, harmonics_theta, harmonics_theta_2, r, theta, phi, tree.expansion_order)
    @show length(system[i_body,POTENTIAL][1].tape) # ~ + 600-1500 tape enrtries
    println("starting expensive part!")
    potential = update_potential!(potential,branch.local_expansion,harmonics,tree.expansion_order)
    @show length(system[i_body,POTENTIAL][1].tape)
    potential_jacobian = update_potential_jacobian!(potential_jacobian,branch.local_expansion,harmonics,harmonics_theta,r,tree.expansion_order)
    @show length(system[i_body,POTENTIAL][1].tape)
    potential_hessian = update_potential_hessian!(potential_hessian,branch.local_expansion,harmonics,harmonics_theta,harmonics_theta_2,r,tree.expansion_order)
    #=for n in 0:tree.expansion_order
        # nm = n * n + n + 1 # m = 0
        nms = (n * (n+1)) >> 1 + 1 # m = 0
        for ind in 1:4
            #potential[ind] += real(branch.local_expansion[ind][nms] * harmonics[nms])
            potential[ind] += branch.local_expansion[ind][nms,1] * harmonics[nms,1] - branch.local_expansion[ind][nms,2] * harmonics[nms,2]
            
            # store derivatives of the potential in spherical coordinates here
            potential_jacobian[1,ind] += n/r * (branch.local_expansion[ind][nms,1] * harmonics[nms,1] - branch.local_expansion[ind][nms,2] * harmonics[nms,2])
            potential_jacobian[2,ind] += branch.local_expansion[ind][nms,1] * harmonics_theta[nms,1] - branch.local_expansion[ind][nms,2] * harmonics_theta[nms,2]
            
            # dJ_potential[3,ind] += 0 # dPsi/dphi
            potential_hessian[1,1,ind] += n * (n-1) / r^2 * (branch.local_expansion[ind][nms,1] * harmonics[nms,1] - branch.local_expansion[ind][nms,2] * harmonics[nms,2])
            potential_hessian[2,1,ind] += n/r * (branch.local_expansion[ind][nms,1] * harmonics_theta[nms,1] - branch.local_expansion[ind][nms,2] * harmonics_theta[nms,2])
            # potential_hessian[3,1,ind] += 0 # d2Psi/dphi dr
            potential_hessian[1,2,ind] += n/r * (branch.local_expansion[ind][nms,1] * harmonics_theta[nms,1] - branch.local_expansion[ind][nms,2] * harmonics_theta[nms,2])
            potential_hessian[2,2,ind] += branch.local_expansion[ind][nms,1] * harmonics_theta_2[nms,1] - branch.local_expansion[ind][nms,2] * harmonics_theta_2[nms,2]
            # potential_hessian[3,2,ind] += 0 # d2Psi/dphi dtheta
            # potential_hessian[1,3,ind] += 0 # d2Psi/dr dphi
            # potential_hessian[2,3,ind] += 0 # d2Psi/dtheta dphi
            # potential_hessian[3,3,ind] += 0 # d2Psi/dphi2
        end
        for m in 1:n # m > 0
            # nm = n * n + n + m + 1
            nms = (n * (n + 1)) >> 1 + m + 1
            for ind in 1:4
                potential[ind] += 2 * (branch.local_expansion[ind][nms,1] * harmonics[nms,1] - branch.local_expansion[ind][nms,2] * harmonics[nms,2])

                # store derivatives of the potential in spherical harmonics here
                potential_jacobian[1,ind] += 2 * n/r * (branch.local_expansion[ind][nms,1] * harmonics[nms,1] - branch.local_expansion[ind][nms,2] * harmonics[nms,2]) # dPsi/dr
                potential_jacobian[2,ind] += 2 * (branch.local_expansion[ind][nms,1] * harmonics_theta[nms,1] - branch.local_expansion[ind][nms,2] * harmonics_theta[nms,2]) # dPsi/dtheta
                potential_jacobian[3,ind] += 2 * m * (branch.local_expansion[ind][nms,1] * harmonics[nms,2] + branch.local_expansion[ind][nms,2] * harmonics[nms,1]) # dPsi/dphi
                potential_hessian[1,1,ind] += 2 * n * (n-1) / r^2 * (branch.local_expansion[ind][nms,1] * harmonics[nms,1] - branch.local_expansion[ind][nms,2] * harmonics[nms,2]) # d2Psi/dr2
                potential_hessian[2,1,ind] += 2 * n/r * (branch.local_expansion[ind][nms,1] * harmonics_theta[nms,1] - branch.local_expansion[ind][nms,2] * harmonics_theta[nms,2]) # d2Psi/dtheta dr
                potential_hessian[3,1,ind] += 2 * n * m / r * (branch.local_expansion[ind][nms,1] * harmonics[nms,2] + branch.local_expansion[ind][nms,2] * harmonics[nms,1]) # d2Psi/dphi dr
                potential_hessian[1,2,ind] += 2 * n/r * (branch.local_expansion[ind][nms,1] * harmonics_theta[nms,1] - branch.local_expansion[ind][nms,2] * harmonics_theta[nms,2]) # d2Psi/dr dtheta
                potential_hessian[2,2,ind] += 2 * (branch.local_expansion[ind][nms,1] * harmonics_theta_2[nms,1] - branch.local_expansion[ind][nms,2] * harmonics_theta_2[nms,2]) # d2Psi/dtheta2
                potential_hessian[3,2,ind] += 2 * m * (branch.local_expansion[ind][nms,1] * harmonics_theta[nms,2] + branch.local_expansion[ind][nms,2] * harmonics_theta[nms,1]) # d2Psi/dphi dtheta
                potential_hessian[1,3,ind] += 2 * n * m / r * (branch.local_expansion[ind][nms,1] * harmonics[nms,2] + branch.local_expansion[ind][nms,2] * harmonics[nms,1]) # d2Psi/dr dphi
                potential_hessian[2,3,ind] += 2 * m * (branch.local_expansion[ind][nms,1] * harmonics_theta[nms,2] + branch.local_expansion[ind][nms,2] * harmonics_theta[nms,1]) # d2Psi/dtheta dphi
                potential_hessian[3,3,ind] += 2 * -m^2 * (branch.local_expansion[ind][nms,1] * harmonics[nms,1] - branch.local_expansion[ind][nms,2] * harmonics[nms,2]) # d2Psi/dphi2
            end
        end
    end=#
    @show length(system[i_body,POTENTIAL][1].tape) # + ~8000 tape entries
    potential_hessian = s2c_hess!(potential_jacobian,potential_hessian,workspace,r,theta,phi)
    potential_jacobian = s2c_jac!(potential_jacobian,workspace,r,theta,phi)
    #potential_jacobian, potential_hessian = spherical_2_cartesian!(potential_jacobian, potential_hessian, workspace, r, theta, phi)
    @show length(system[i_body,POTENTIAL][1].tape) # + ~1500 tape entries # 2 tape entries now.
    potential_jacobian = flatten_jacobian!(potential_jacobian)
    potential_hessian = flatten_hessian!(potential_hessian)
    #flatten_derivatives!(potential_jacobian, potential_hessian) # compute velocity and velocity gradient
    system[i_body,POTENTIAL] += potential
    system[i_body,VELOCITY] += potential_jacobian[:,1]
    system[i_body,VELOCITYGRADIENT] += potential_hessian[:,:,1]
    @show length(system[i_body,POTENTIAL][1].tape) # + 80 tape entries, probably corresponding to the += operations. # down to 23 entries for some reason... this is concerning.
end