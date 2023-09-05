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

function spherical_2_cartesian!(potential_jacobian, potential_hessian, workspace, rho, theta, phi)
    # get partial derivatives of the coordinates
    s_theta, c_theta = sincos(theta)
    s_phi, c_phi = sincos(phi)

    drjdxi = @SMatrix [
        s_theta*c_phi c_theta*c_phi/rho -s_phi/rho/s_theta
        s_theta * s_phi c_theta * s_phi / rho c_phi / rho / s_theta
        c_theta -s_theta / rho 0
    ]

    # convert Hessian to cartesian coordinates
    workspace3x3 = view(workspace,:,1:3)
    for ind in 1:4
        workspace3x3 .= potential_hessian[:,:,ind]
        potential_hessian[:,:,ind] .= drjdxi * workspace3x3 * transpose(drjdxi)
        for k_coord in 1:3 # loop over r, theta, and phi to save me some space on the stack
            if k_coord == 1 # r coordinate
                drkdxidxj = @SMatrix [
                    (1-cos(phi)^2 * sin(theta)^2)/rho -sin(theta)^2*cos(phi)*sin(phi)/rho -sin(theta)*cos(phi)*cos(theta)/rho;
                    (-sin(theta)^2*cos(phi)*sin(phi))/rho (1-sin(theta)^2*sin(phi)^2)/rho -sin(theta)*sin(phi)*cos(theta)/rho;
                    -sin(theta)*cos(phi)*cos(theta)/rho -sin(theta)*sin(phi)*cos(theta)/rho sin(theta)^2/rho
                ]
            elseif k_coord == 2 # theta coordinate
                drkdxidxj = @SMatrix [
                    cos(theta)/sin(theta)*(1-cos(phi)^2*(1+2*sin(theta)^2))/rho^2 -cos(theta)/sin(theta)*sin(phi)*cos(phi)*(1+2*sin(theta)^2)/rho^2 cos(phi)*(1-2*cos(theta)^2)/rho^2;
                    -cos(theta)/sin(theta)*sin(phi)*cos(phi)*(1+2*sin(theta)^2)/rho^2 cos(theta)/sin(theta)*(1-sin(phi)^2*(1+2*sin(theta)^2))/rho^2 (2*sin(theta)^2-1)/rho^2*sin(phi);
                    cos(phi)*(1-2*cos(theta)^2)/rho^2 (2*sin(theta)^2-1)/rho^2*sin(phi) 2*sin(theta)*cos(theta)/rho^2
                ]
            else # phi coordinate
                drkdxidxj = @SMatrix [
                    2*cos(phi)*sin(phi)/rho^2/sin(theta)^2 (2*sin(phi)^2-1)/rho^2/sin(theta)^2 0;
                    (2*sin(phi)^2-1)/rho^2/sin(theta)^2 -2*sin(phi)*cos(phi)/rho^2/sin(theta)^2 0;
                    0 0 0
                ]
            end
            potential_hessian[:,:,ind] .+= drkdxidxj * potential_jacobian[k_coord,ind]
        end
    end

    workspace .= potential_jacobian

    mul!(potential_jacobian, drjdxi, workspace)

    return nothing
end

@inline odd_or_even(n::Int) = (n & 1) == 1 ? -1 : 1

@inline ipow2l(n::Int) = n >= 0 ? 1 : odd_or_even(n);

function regular_harmonic!(harmonics, harmonics_theta, harmonics_theta_2, rho, theta, phi, P)
    y,x = sincos(theta)
    invY = y == 0 ? 0 : 1 / y
    fact = 1
    pl = 1
    rhom = 1
    ei = exp(im * phi)
    eim = 1.0
    for m=0:P
        p = pl
        lpl = (m * (m + 1)) >> 1 + m + 1
        # lpl = m * m + 2 * m + 1
        # lml = m * m + 1
        harmonics[lpl] = rhom * p * eim
        # harmonics[lml] = conj(harmonics[lpl])
        p1 = p
        p = x * (2 * m + 1) * p1
        harmonics_theta[lpl] = rhom * (p - (m + 1) * x * p1) * invY * eim
        harmonics_theta_2[lpl] = rhom * (-x * p + (-m + (m+1)^2 * x^2) * p1) * invY^2 * eim
        rhom *= rho
        rhol = rhom
        for l=m+1:P
            lpm = (l * (l + 1)) >> 1 + m + 1
            # lpm = l * l + l + m + 1
            # lmm = l * l + l - m + 1
            rhol /= -(l + m)
            harmonics[lpm] = rhol * p * eim
            # harmonics[lmm] = conj(harmonics[lpm])
            p2 = p1
            p1 = p
            p = (x * (2 * l + 1) * p1 - (l + m) * p2) / (l - m + 1)
            harmonics_theta[lpm] = rhol * ((l - m + 1) * p - (l + 1) * x * p1) * invY * eim
            harmonics_theta_2[lpm] = rhol * ((m-l-1) * x * p + (m^2 - l*(l+1) + (l+1)^2 * x^2) * p1) * invY^2 * eim
            rhol *= rho
        end
        rhom /= -(2 * m + 2) * (2 * m + 1)
        pl = -pl * fact * y
        fact += 2
        eim *= ei
    end
end

function regular_harmonic!(harmonics, rho, theta, phi, P)
    y,x = sincos(theta)
    fact = 1
    pl = 1
    rhom = 1 # rho^l / (l+m)! * (-1)^l
    ei = exp(im * phi)
    eim = 1.0 # e^(i * m * phi)
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
end

function B2M!(branch, system, bodies_index, harmonics, expansion_order)
    @warn "B2M! not implemented for type $(typeof(system)); overload FLOWFMM.B2M!"
end

function B2M!(tree, systems, i_branch, sources_index)
    branch = tree.branches[i_branch]
    harmonics = Vector{eltype(branch.multipole_expansion[1])}(undef, (tree.expansion_order+1)^2) # this is faster than MVector

    # iterate over elements
    for (i_iter, system) in enumerate(systems[sources_index])
        i_type = sources_index[i_iter]
        bodies_index = branch.first_body[i_type]:branch.first_body[i_type] + branch.n_bodies[i_type] - 1
        B2M!(branch, system, bodies_index, harmonics, tree.expansion_order)
    end
end

function M2B!(target_potential, target, i_branch, tree)
    branch = tree.branches[i_branch]
    irregular_harmonics = Vector{eltype(branch.multipole_expansion[1])}(undef, (tree.expansion_order+1)^2)
    dx = target[1:3] - branch.center
    r, theta, phi = cartesian_2_spherical(dx)
    irregular_harmonic!(irregular_harmonics, r, theta, phi, tree.expansion_order)
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

function M2M!(tree, branch, child, harmonics)
    # get distance vector
    dx, dy, dz = branch.center - child.center
    r, theta, phi = cartesian_2_spherical(dx, dy, dz)
    regular_harmonic!(harmonics, r, theta, phi, tree.expansion_order)

    M = zeros(eltype(branch.multipole_expansion[1]), 4)
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
                        M[dim] += child.multipole_expansion[dim][jlkms] * harmonics[lm] * ipow * oddeven
                    end
                end
                for m in k:min(l,j+k-l)
                    jlkms = (((j-l) * (j-l+1)) >> 1) - k + m + 1
                    lm = l * l + l - m + 1
                    oddeven = odd_or_even(k + l + m)
                    for dim in 1:4
                        M[dim] += conj(child.multipole_expansion[dim][jlkms]) * harmonics[lm] * oddeven
                    end
                end
            end
            for dim in 1:4
                branch.multipole_expansion[dim][i_jk] += M[dim]
            end
        end
    end
end

function M2M!(tree, i_branch)
    # expose objects
    branch = tree.branches[i_branch]

    #initialize memory TODO: do this beforehand?
    harmonics = Vector{eltype(branch.multipole_expansion[1])}(undef, (tree.expansion_order+1)^2)

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
                        L[dim] += conj(multipole_expansion[dim][nms]) * Cnm * harmonics[jnkm]
                    end
                end
                for m in 0:n
                    nms = (n * (n+1)) >> 1 + m + 1
                    jnkm = (j + n) * (j + n) + j + n + m - k + 1
                    # jnkm_max = 2P * 2P + 2P + P + P - 0 + 1 = (2P)^2 + 2P + 2P + 1 = 4P^2 + 4P + 1 = (2P + 1)^2
                    Cnm2 = Cnm * odd_or_even((k-m) * (1 >> (k>=m)) + m)
                    for dim in 1:4
                        L[dim] += multipole_expansion[dim][nms] * Cnm2 * harmonics[jnkm]
                    end
                end
            end
            for dim in 1:4
                local_expansion[dim][jks] += L[dim]
            end
        end
    end
end

function M2L!(tree, i_local, j_multipole)
    local_branch = tree.branches[i_local]
    multipole_branch = tree.branches[j_multipole]
    # preallocate
    harmonics = Vector{eltype(local_branch.multipole_expansion[1])}(undef, (2*tree.expansion_order + 1)^2)
    # get separation vector
    dx, dy, dz = local_branch.center - multipole_branch.center
    r, theta, phi = cartesian_2_spherical(dx, dy, dz)
    irregular_harmonic!(harmonics, r, theta, phi, 2*tree.expansion_order)
    L = zeros(eltype(local_branch.local_expansion[1]), 4)
    M2L_loop!(local_branch.local_expansion, L, multipole_branch.multipole_expansion, harmonics, tree.expansion_order)
end

function B2L!(tree, i_branch, source_position, source_strength)
    branch = tree.branches[i_branch]
    irregular_harmonics = Vector{eltype(branch.multipole_expansion[1])}(undef, (tree.expansion_order+1)^2)
    r, theta, phi = cartesian_2_spherical(source_position - branch.center)
    irregular_harmonic!(irregular_harmonics, r, theta, -phi, tree.expansion_order)
    for l in 0:tree.expansion_order
        for m in 0:l
            i_abb = (l * (l+1)) >> 1 + m + 1
            i_exp = l^2 + l + m + 1
            for dim in 1:4
                branch.local_expansion[dim][i_abb] = irregular_harmonics[i_exp] * source_strength[dim]
            end
        end
    end
end

function L2L!(tree, branch, child, harmonics)
    dx, dy, dz = child.center - branch.center
    r, theta, phi = cartesian_2_spherical(dx, dy, dz)
    regular_harmonic!(harmonics, r, theta, phi, tree.expansion_order)
    for j in 0:tree.expansion_order
        for k in 0:j
            jks = (j * (j + 1)) >> 1 + k + 1
            L = zeros(eltype(branch.local_expansion[1]), 4)
            for n in j:tree.expansion_order
                for m in j+k-n:-1
                    jnkm = (n-j) * (n-j) + n - j + m - k + 1
                    nms = (n * (n + 1)) >> 1 - m + 1
                    oddeven = odd_or_even(k)
                    for dim in 1:4
                        L[dim] += conj(branch.local_expansion[dim][nms]) * harmonics[jnkm] * oddeven
                    end
                end
                for m in 0:n
                    if n-j >= abs(m-k)
                        jnkm = (n - j) * (n - j) + n - j + m - k + 1
                        nms = (n * (n + 1)) >> 1 + m + 1
                        oddeven = odd_or_even((m-k) * (1 >> (m >= k)))
                        for dim in 1:4
                            L[dim] += branch.local_expansion[dim][nms] * harmonics[jnkm] * oddeven
                        end
                    end
                end
            end
            for dim in 1:4
                child.local_expansion[dim][jks] += L[dim]
            end
        end
    end
end

function L2L!(tree, j_source)
    # expose branch
    branch = tree.branches[j_source]

    #initialize memory TODO: do this beforehand?
    harmonics = Vector{eltype(branch.multipole_expansion[1])}(undef, (tree.expansion_order+1)^2)

    # iterate over children
    for i_child in branch.first_branch:branch.first_branch + branch.n_branches - 1
        child = tree.branches[i_child]
        L2L!(tree, branch, child, harmonics)
    end
end

"Calculates the potential at all child elements of a branch."
function L2B!(tree, systems, i_branch, targets_index)
    branch = tree.branches[i_branch]
    harmonics = Vector{eltype(branch.multipole_expansion[1])}(undef, ((tree.expansion_order+1) * (tree.expansion_order+2)) >> 1)
    harmonics_theta = zeros(eltype(branch.multipole_expansion[1]), ((tree.expansion_order+1) * (tree.expansion_order+2)) >> 1)
    harmonics_theta_2 = zeros(eltype(branch.multipole_expansion[1]), ((tree.expansion_order+1) * (tree.expansion_order+2)) >> 1)
    workspace = zeros(eltype(branch.multipole_expansion[1]),3,4)
    spherical_potential = zeros(eltype(branch.multipole_expansion[1]),52)
    for (i_target, system) in enumerate(systems[targets_index])
        i_type = targets_index[i_target]
        for i_body in branch.first_body[i_type]:branch.first_body[i_type] + branch.n_bodies[i_type] - 1
            L2B!(system, i_body, harmonics, harmonics_theta, harmonics_theta_2, workspace, spherical_potential, tree, branch)
            spherical_potential .*= 0
        end
    end
end

@inline function L2B!(system, i_body, harmonics, harmonics_theta, harmonics_theta_2, workspace, spherical_potential, tree, branch)
    potential = view(spherical_potential,1:4)
    potential_jacobian = reshape(view(spherical_potential, 5:16),3,4)
    potential_hessian = reshape(view(spherical_potential, 17:52),3,3,4)
    body_position = system[i_body,POSITION]
    dx, dy, dz = body_position - branch.center
    r, theta, phi = cartesian_2_spherical(dx, dy, dz)
    regular_harmonic!(harmonics, harmonics_theta, harmonics_theta_2, r, theta, phi, tree.expansion_order)
    for n in 0:tree.expansion_order
        # nm = n * n + n + 1 # m = 0
        nms = (n * (n+1)) >> 1 + 1 # m = 0
        for ind in 1:4
            potential[ind] += real(branch.local_expansion[ind][nms] * harmonics[nms])
            # store derivatives of the potential in spherical coordinates here
            potential_jacobian[1,ind] += n/r * real(branch.local_expansion[ind][nms] * harmonics[nms]) # dPsi/dr
            potential_jacobian[2,ind] += real(branch.local_expansion[ind][nms] * harmonics_theta[nms]) # dPsi/dtheta
            # dJ_potential[3,ind] += 0 # dPsi/dphi
            potential_hessian[1,1,ind] += n * (n-1) / r^2 * real(branch.local_expansion[ind][nms] * harmonics[nms]) # d2Psi/dr2
            potential_hessian[2,1,ind] += n/r * real(branch.local_expansion[ind][nms] * harmonics_theta[nms]) # d2Psi/dtheta dr
            # potential_hessian[3,1,ind] += 0 # d2Psi/dphi dr
            potential_hessian[1,2,ind] += n/r * real(branch.local_expansion[ind][nms] * harmonics_theta[nms]) # d2Psi/dr dtheta
            potential_hessian[2,2,ind] += real(branch.local_expansion[ind][nms] * harmonics_theta_2[nms]) # d2Psi/dtheta2
            # potential_hessian[3,2,ind] += 0 # d2Psi/dphi dtheta
            # potential_hessian[1,3,ind] += 0 # d2Psi/dr dphi
            # potential_hessian[2,3,ind] += 0 # d2Psi/dtheta dphi
            # potential_hessian[3,3,ind] += 0 # d2Psi/dphi2
        end
        for m in 1:n # m > 0
            # nm = n * n + n + m + 1
            nms = (n * (n + 1)) >> 1 + m + 1
            for ind in 1:4
                potential[ind] += 2 * real(branch.local_expansion[ind][nms] * harmonics[nms])
                # store derivatives of the potential in spherical harmonics here
                potential_jacobian[1,ind] += 2 * n/r * real(branch.local_expansion[ind][nms] * harmonics[nms]) # dPsi/dr
                potential_jacobian[2,ind] += 2 * real(branch.local_expansion[ind][nms] * harmonics_theta[nms]) # dPsi/dtheta
                potential_jacobian[3,ind] += 2 * m * real(im * branch.local_expansion[ind][nms] * harmonics[nms]) # dPsi/dphi
                potential_hessian[1,1,ind] += 2 * n * (n-1) / r^2 * real(branch.local_expansion[ind][nms] * harmonics[nms]) # d2Psi/dr2
                potential_hessian[2,1,ind] += 2 * n/r * real(branch.local_expansion[ind][nms] * harmonics_theta[nms]) # d2Psi/dtheta dr
                potential_hessian[3,1,ind] += 2 * n * m / r * real(im * branch.local_expansion[ind][nms] * harmonics[nms]) # d2Psi/dphi dr
                potential_hessian[1,2,ind] += 2 * n/r * real(branch.local_expansion[ind][nms] * harmonics_theta[nms]) # d2Psi/dr dtheta
                potential_hessian[2,2,ind] += 2 * real(branch.local_expansion[ind][nms] * harmonics_theta_2[nms]) # d2Psi/dtheta2
                potential_hessian[3,2,ind] += 2 * m * real(im * branch.local_expansion[ind][nms] * harmonics_theta[nms]) # d2Psi/dphi dtheta
                potential_hessian[1,3,ind] += 2 * n * m / r * real(im * branch.local_expansion[ind][nms] * harmonics[nms]) # d2Psi/dr dphi
                potential_hessian[2,3,ind] += 2 * m * real(im * branch.local_expansion[ind][nms] * harmonics_theta[nms]) # d2Psi/dtheta dphi
                potential_hessian[3,3,ind] += 2 * -m^2 * real(branch.local_expansion[ind][nms] * harmonics[nms]) # d2Psi/dphi2
            end
        end
    end
    spherical_2_cartesian!(potential_jacobian, potential_hessian, workspace, r, theta, phi)
    system[i_body,POTENTIAL] += potential
    system[i_body,JACOBIAN] += potential_jacobian
    system[i_body,HESSIAN] += potential_hessian
end
