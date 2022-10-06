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

function regular_harmonic!(harmonics, rho, theta, phi, p)
    y, x = sincos(theta)
    inv_y = y !== 0.0 ? 1/y : 0.0

    p_m_m = 1.0 # start with l=m=0
    p_m1_m = x

    eim = 1 # e^{i m \phi}
    ei = exp(im * phi)

    mm! = 1 # sqrt((m-m)!/(m+m)!)
    rhom = 1 # rho^l when l=m
    m21 = 1 # 2m + 1

    for m in 0:p
        # l = m iteration
        ip = m^2 + m + m + 1
        im = m^2 + m - m + 1
        harmonics[ip] = mm! * p_m_m * eim * rhom
        harmonics[im] = conj(harmonics[ip])

        # update values for l = m+1

        p_m1_m = x * m21 * p_m_m
        p_l_1_m = p_m_m
        p_l_m = p_m1_m

        # set up variables for inner loop
        rhol = rhom
        lm! = mm! # sqrt((l-m)!/(l+m)!)
        for l in m+1:p
            # update rho^l
            rhol *= rho

            # update sqrt((l-m)!/(l+m)!)
            lm! *= sqrt((l-m)/(l+m))

            ip = l^2 + l + m + 1
            im = l^2 + l - m + 1
            harmonics[ip] = lm! * p_l_m * eim * rhol
            harmonics[im] = conj(harmonics[ip])

            # update legendre polynomials
            p_tmp = ((2*l+1) * x * p_l_m - (l+m) * p_l_1_m) / (l-m+1)
            p_l_1_m = p_l_m
            p_l_m = p_tmp
        end

        # increment legendre polynomials order
        p_m_m *= -m21 * y

        # update e^{i m \phi}
        eim *= ei

        # update rho^m
        rhom *= rho

        # update sqrt((l-m)!/(l+m)!) for l=m
        mm! /= sqrt((m21 + 1) * m21)

        # update 2m + 1
        m21 += 2
    end
end

function irregular_harmonic!(harmonics, rho, theta, phi, p)
    y, x = sincos(theta)
    inv_y = y !== 0.0 ? 1/y : 0.0

    p_m_m = 1.0 # start with l=m=0
    p_m1_m = x

    eim = 1 # e^{i m \phi}
    ei = exp(im * phi)

    mm! = 1 # sqrt((m-m)!/(m+m)!)
    inv_rho = rho !== 0 ? 1/rho : 0.0 # 1/rho
    rhom1 = inv_rho # rho^{-l-1} when l=m
    m21 = 1 # 2m + 1

    for m in 0:p
        # l = m iteration
        ip = m^2 + m + m + 1
        im = m^2 + m - m + 1
        harmonics[ip] = mm! * p_m_m * eim * rhom1
        harmonics[im] = conj(harmonics[ip])

        # update values for l = m+1
        p_m1_m = x * m21 * p_m_m
        p_l_1_m = p_m_m
        p_l_m = p_m1_m

        # set up variables for inner loop
        rhol1 = rhom1
        lm! = mm! # sqrt((l-m)!/(l+m)!)
        for l in m+1:p
            # update rho^{-l-1}
            rhol1 *= inv_rho

            # update sqrt((l-m)!/(l+m)!)
            lm! *= sqrt((l-m)/(l+m))

            ip = l^2 + l + m + 1
            im = l^2 + l - m + 1
            harmonics[ip] = lm! * p_l_m * eim * rhol1
            harmonics[im] = conj(harmonics[ip])

            # update legendre polynomials
            p_tmp = ((2*l+1) * x * p_l_m - (l+m) * p_l_1_m) / (l-m+1)
            p_l_1_m = p_l_m
            p_l_m = p_tmp
        end

        # increment legendre polynomials order
        p_m_m *= -m21 * y

        # update e^{i m \phi}
        eim *= ei

        # update rho^m
        rhom1 *= inv_rho

        # update sqrt((l-m)!/(l+m)!) for l=m
        mm! /= sqrt((m21 + 1) * m21)

        # update 2m + 1
        m21 += 2
    end
end

function P2M!(tree, branch, element, solid_harmonics, ::Spherical)
    dx = get_x(element) .- branch.center
    cartesian_2_spherical!(dx)
    regular_harmonic!(solid_harmonics, dx[1], dx[2], -dx[3], tree.expansion_order[1]) # Ylm^* -> -dx[3]

    # update values
    for l in 0:tree.expansion_order[1]
        for m in 0:l
            i_solid_harmonic = l^2 + l + m + 1
            i_compressed = 1 + (l * (l + 1)) >> 1 + m # only save half as Yl{-m} = conj(Ylm)
            branch.multipole_expansion[i_compressed] += solid_harmonics[i_solid_harmonic] * get_q(element)
        end
    end
end

function P2M!(i_branch, tree, elements, coordinates::Spherical)
    branch = tree.branches[i_branch]

    #initialize memory TODO: do this beforehand?
    solid_harmonics = Vector{Complex{Float64}}(undef, (tree.expansion_order[1]+1)^2)

    # iterate over elements
    for i_element in tree.indices[branch.first_element:branch.first_element + branch.n_elements-1]
        element = elements[i_element]
        P2M!(tree, branch, element, solid_harmonics, coordinates)
    end
end

function M2P!(target, i_branch, tree)
    branch = tree.branches[i_branch]
    irregular_harmonics = Vector{Complex{Float64}}(undef, (tree.expansion_order[1]+1)^2)
    dx = get_x(target) - branch.center
    cartesian_2_spherical!(dx)
    irregular_harmonic!(irregular_harmonics, dx[1], dx[2], dx[3], tree.expansion_order[1])
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

function M2M!(tree, branch, child, Ylm, ::Spherical)
    # get distance vector
    dx = branch.center - child.center
    cartesian_2_spherical!(dx)
    regular_harmonic!(Ylm, dx[1], dx[2], dx[3], tree.expansion_order[1])

    for j in 0:tree.expansion_order[1] # iterate over new Multipole coefficients B_j^k
        for k in 0:j
            i_jk = ((j * (j+1)) >> 1) + k + 1 # current index
            M = 0.0 # vectorize later
            for l in 0:j
                for m in max(-l,-j+k+l):min(k-1,l)
                    jlkms = (((j-l) * (j-l+1)) >> 1) + k - m + 1
                    lm = l * l + l - m + 1
                    M += child.multipole_expansion[jlkms] * Ylm[lm] * ipow2l(m) * odd_or_even(l)
                end
                for m in k:min(l,j+k-l)
                    jlkms = (((j-l) * (j-l+1)) >> 1) - k + m + 1
                    lm = l * l + l - m + 1
                    M += conj(child.multipole_expansion[jlkms]) * Ylm[lm] * odd_or_even(k + l + m)
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
    Ylm = Vector{Complex{Float64}}(undef, (tree.expansion_order[1]+1)^2)
    # YlmTheta = Vector{Complex{Float64}}(undef, tree.expansion_order[1]^2)

    # iterate over children
    for i_child in branch.first_branch:branch.first_branch + branch.n_branches - 1
        child = tree.branches[i_child]
        M2M!(tree, branch, child, Ylm, coordinates)
    end
end

function M2M!(i_branch, tree, ::Spherical)
    # expose objects
    branch = tree.branches[i_branch]

    #initialize memory TODO: do this beforehand?
    Ylm = Vector{Complex{Float64}}(undef, (tree.expansion_order[1]+1)^2)
    # YlmTheta = Vector{Complex{Float64}}(undef, tree.expansion_order[1]^2)

    # iterate over children
    for i_child in branch.first_branch:branch.first_branch + branch.n_branches - 1
        child = tree.branches[i_child]

        # get distance vector
        dx = branch.center - child.center
        cartesian_2_spherical!(dx)
        regular_harmonic!(Ylm, dx[1], dx[2], dx[3], tree.expansion_order[1])

        for j in 0:tree.expansion_order[1]
            pass
        end
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
