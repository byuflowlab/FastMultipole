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
    return [r, theta, phi]
end

@inline odd_or_even(n::Int) = (n & 1) == 1 ? -1 : 1

@inline ipow2l(n::Int) = n >= 0 ? 1 : odd_or_even(n);

function regular_harmonic!(harmonics, harmonics_theta, rho, theta, phi, P)
    y,x = sincos(theta)
    invY = y == 0 ? 0 : 1 / y
    fact = 1
    pl = 1
    rhom = 1
    ei = exp(im * phi)
    eim = 1.0
    for m=0:P
        p = pl
        lpl = m * m + 2 * m + 1
        lml = m * m + 1
        harmonics[lpl] = rhom * p * eim
        harmonics[lml] = conj(harmonics[lpl])
        p1 = p
        p = x * (2 * m + 1) * p1
        harmonics_theta[lpl] = rhom * (p - (m + 1) * x * p1) * invY * eim
        rhom *= rho
        rhol = rhom
        for l=m+1:P
            lpm = l * l + l + m + 1
            lmm = l * l + l - m + 1
            rhol /= -(l + m)
            harmonics[lpm] = rhol * p * eim
            harmonics[lmm] = conj(harmonics[lpm])
            p2 = p1
            p1 = p
            p = (x * (2 * l + 1) * p1 - (l + m) * p2) / (l - m + 1)
            harmonics_theta[lpm] = rhol * ((l - m + 1) * p - (l + 1) * x * p1) * invY * eim
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

function P2M!(tree, branch, element, harmonics, harmonics_theta, ::Spherical)
    dx = get_x(element) .- branch.center
    cartesian_2_spherical!(dx)
    regular_harmonic!(harmonics, harmonics_theta, dx[1], dx[2], -dx[3], tree.expansion_order[1]) # Ylm^* -> -dx[3]
    # update values
    for l in 0:tree.expansion_order[1]
        for m in 0:l
            i_solid_harmonic = l^2 + l + m + 1
            i_compressed = 1 + (l * (l + 1)) >> 1 + m # only save half as Yl{-m} = conj(Ylm)
            branch.multipole_expansion[i_compressed] += harmonics[i_solid_harmonic] * get_q(element)
        end
    end
end

function P2M!(i_branch, tree, elements, coordinates::Spherical)
    branch = tree.branches[i_branch]

    #initialize memory TODO: do this beforehand?
    harmonics = Vector{Complex{Float64}}(undef, (tree.expansion_order[1]+1)^2)
    harmonics_theta = Vector{Complex{Float64}}(undef, (tree.expansion_order[1]+1)^2)

    # iterate over elements
    for i_element in tree.indices[branch.first_element:branch.first_element + branch.n_elements-1]
        element = elements[i_element]
        P2M!(tree, branch, element, harmonics, harmonics_theta, coordinates)
    end
end

function M2P!(target, i_branch, tree)
    branch = tree.branches[i_branch]
    irregular_harmonics = Vector{Complex{Float64}}(undef, (tree.expansion_order[1]+1)^2)
    dx = get_x(target) - branch.center
    cartesian_2_spherical!(dx)
    irregular_harmonic!(irregular_harmonics, dx..., tree.expansion_order[1])
    for l in 0:tree.expansion_order[1]
        for m in 0:l
            ip = l^2 + l + m + 1
            i_compressed = 1 + (l * (l + 1)) >> 1 + m # only save half as Yl{-m} = conj(Ylm)
            d_potential = real(branch.multipole_expansion[i_compressed] * irregular_harmonics[ip])
            m > 0 && (d_potential *= 2)
            target.potential[1] += d_potential
        end
    end
end

function M2M!(tree, branch, child, harmonics, harmonics_theta, ::Spherical)
    # get distance vector
    dx = branch.center - child.center
    cartesian_2_spherical!(dx)
    regular_harmonic!(harmonics, harmonics_theta, dx..., tree.expansion_order[1])

    for j in 0:tree.expansion_order[1] # iterate over new Multipole coefficients B_j^k
        for k in 0:j
            i_jk = ((j * (j+1)) >> 1) + k + 1 # current index
            M = 0.0 # vectorize later
            for l in 0:j
                for m in max(-l,-j+k+l):min(k-1,l)
                    jlkms = (((j-l) * (j-l+1)) >> 1) + k - m + 1
                    lm = l * l + l - m + 1
                    M += child.multipole_expansion[jlkms] * harmonics[lm] * ipow2l(m) * odd_or_even(l)
                end
                for m in k:min(l,j+k-l)
                    jlkms = (((j-l) * (j-l+1)) >> 1) - k + m + 1
                    lm = l * l + l - m + 1
                    M += conj(child.multipole_expansion[jlkms]) * harmonics[lm] * odd_or_even(k + l + m)
                end
            end
            # if i_jk == 4; println("M2M: l = $j, m = $k\n\tM = $M"); end
            branch.multipole_expansion[i_jk] += M
        end
    end
end

function M2M!(i_branch, tree, coordinates::Spherical)
    # expose objects
    branch = tree.branches[i_branch]

    #initialize memory TODO: do this beforehand?
    harmonics = Vector{Complex{Float64}}(undef, (tree.expansion_order[1]+1)^2)
    harmonics_theta = Vector{Complex{Float64}}(undef, (tree.expansion_order[1]+1)^2)

    # iterate over children
    for i_child in branch.first_branch:branch.first_branch + branch.n_branches - 1
        child = tree.branches[i_child]
        M2M!(tree, branch, child, harmonics, harmonics_theta, coordinates)
    end
end

function M2L!(i_local, j_multipole, tree, elements, derivatives, ::Spherical)
    local_branch = tree.branches[i_local]
    multipole_branch = tree.branches[j_multipole]

    # preallocate
    Ynm2 = Vector{Complex{Float64}}(undef,4 * tree.expansion_order[1]^2)

    # get separation vector
    dx = local_branch.center - multipole_branch.center
    cartesian_2_spherical!(dx)
    evaluate_local!(Ynm2, dx[1], dx[2], dx[3], tree.expansion_order[1])
    for j in 0:tree.expansion_order[1]-1
        Cnm = odd_or_even(j)
        for k in 0:j
            jks = div(j * (j + 1), 2) + k
            L = 0.0
            for n in 0:tree.expansion_order[1]-1
                for m in -n:-1
                    nms = div(n * (n+1), 2) - m
                    jnkm = (j + n) * (j + n) + j + n + m - k
                    L += conj(multipole_branch.multipole_expansion[nms+1]) * Cnm * Ynm2[jnkm+1]
                end
                for m in 0:n
                    nms = div(n * (n+1), 2) + m
                    jnkm = (j + n) * (j + n) + j + n + m - k
                    Cnm2 = Cnm * odd_or_even((k-m) * (k<m) + m) # TODO: fix syntax
                    L += multipole_branch.multipole_expansion[nms+1] * Cnm2 * Ynm2[jnkm+1]
                end
            end
            local_branch.local_expansion[jks+1] += L
        end
    end
end

function L2L!(j_source, tree, ::Spherical)
    # expose branch
    branch = tree.branches[j_source]

    #initialize memory TODO: do this beforehand?
    Ynm = Vector{Complex{Float64}}(undef, tree.expansion_order[1]^2)
    # YnmTheta = Vector{Complex{Float64}}(undef, tree.expansion_order[1]^2)

    # iterate over children
    for i_child in branch.first_branch:branch.first_branch + branch.n_branches - 1
        child = tree.branches[i_child]
        dx = child.center - branch.center
        cartesian_2_spherical!(dx)
        evaluate_multipole!(Ynm, dx[1], dx[2], dx[3], tree.expansion_order[1])
        for j in 0:tree.expansion_order[1]-1
            for k in 0:j
                jks = ((j * (j + 1)) >> 1) + k
                L = 0.0
                for n in j:tree.expansion_order[1]-1
                    for m in j+k-n:-1
                        jnkm = (n-j) * (n-j) + n - j + m - k
                        nms = ((n * (n + 1)) >> 1) - m
                        L += conj(branch.local_expansion[nms+1]) * Ynm[jnkm+1] * odd_or_even(k)
                    end
                    for m in 0:n
                        if n-j >= abs(m-k)
                            jnkm = (n - j) * (n - j) + n - j + m - k
                            nms = ((n * (n + 1)) >> 1) + m
                            L += branch.local_expansion[nms+1] * Ynm[jnkm+1] * odd_or_even((m-k) * (m < k))
                        end
                    end
                end
                child.local_expansion[jks+1] += L
            end
        end
    end
end

"Calculates the potential at all child elements of a branch."
function L2P!(i_branch, tree, elements, ::Spherical)
    branch = tree.branches[i_branch]
    Ynm = Vector{Complex{Float64}}(undef, tree.expansion_order[1]^2)
    # YnmTheta = Vector{Complex{Float64}}(undef, tree.expansion_order[1]^2)
    for i_element in branch.first_element:branch.first_element + branch.n_elements - 1
        element = elements[tree.indices[i_element]]
        dx = get_x(element) - branch.center
        cartesian_2_spherical!(dx)
        evaluate_multipole!(Ynm, dx[1], dx[2], dx[3], tree.expansion_order[1])
        for n in 0:tree.expansion_order[1]-1
            nm = n * n + n
            nms = (n * (n+1)) >> 1
            element.potential .+= real(branch.local_expansion[nms+1] * Ynm[nm+1])
            for m in 1:n
                nm = n * n + n + m
                nms = (n * (n + 1)) >> 1 + m
                element.potential .+= 2 * real(branch.local_expansion[nms+1] * Ynm[nm+1])
            end
        end
    end
end
