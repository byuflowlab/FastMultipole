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

function spherical_2_cartesian!(dx_cartesian, dx_spherical, aux_spherical)
    s_theta, c_theta = sincos(dx_spherical[2])
    s_phi, c_phi = sincos(dx_spherical[3])
    r, theta, phi = dx_spherical
    dx_cartesian[1] = s_theta * c_phi * aux_spherical[1]
        + c_theta * c_phi / dx_spherical[1] * aux_spherical[1]
        - s_phi / r / s_theta * aux_spherical[3]
    dx_cartesian[2] = s_theta * s_phi * aux_spherical[1]
        + c_theta * s_phi / r * aux_spherical[2]
        + c_phi / r / s_theta * aux_spherical[3]
    dx_cartesian[3] = c_theta * aux_spherical[1]
        - s_theta / r * aux_spherical[2]

    return nothing
end

function spherical_2_cartesian(dx_spherical, o_spherical)
    s_theta, c_theta = sincos(dx_spherical[2])
    s_phi, c_phi = sincos(dx_spherical[3])
    dx_cartesian_1 = s_theta * c_phi * o_spherical[1]
        + c_theta * c_phi / dx_spherical[1] * o_spherical[1]
        - s_phi / r / s_theta * o_spherical[3]
    dx_cartesian_2 = s_theta * s_phi * o_spherical[1]
        + c_theta * s_phi / r * o_spherical[2]
        + c_phi / r / s_theta * o_spherical[3]
    dx_cartesian_3 = c_theta * o_spherical[1]
        - s_theta / r * o_spherical[2]

    return [dx_cartesian_1, dx_cartesian_2, dx_cartesian_3]
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

function regular_harmonic!(harmonics, rho, theta, phi, P)
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

function P2M!(tree, elements, i_branch)
    branch = tree.branches[i_branch]

    #initialize memory TODO: do this beforehand?
    harmonics = Vector{Complex{Float64}}(undef, (tree.expansion_order[1]+1)^2)

    # iterate over elements
    for i_element in tree.indices[branch.first_element:branch.first_element + branch.n_elements-1]
        element = elements[i_element]
        tree.P2M!(tree, branch, element, harmonics)
    end
end

function M2P!(target, i_branch, tree)
    branch = tree.branches[i_branch]
    irregular_harmonics = Vector{Complex{Float64}}(undef, (tree.expansion_order[1]+1)^2)
    dx = target.position - branch.center
    cartesian_2_spherical!(dx)
    irregular_harmonic!(irregular_harmonics, dx..., tree.expansion_order[1])
    d_potential = zeros(4)
    for l in 0:tree.expansion_order[1]
        for m in 0:l
            ip = l^2 + l + m + 1
            i_compressed = 1 + (l * (l + 1)) >> 1 + m # only save half as Yl{-m} = conj(Ylm)
            for dim in 1:4
                d_potential[dim] = real(branch.multipole_expansion[dim][i_compressed] * irregular_harmonics[ip])
            end
            m > 0 && (d_potential .*= 2)
            target.potential .+= d_potential
        end
    end
end

function M2M!(tree, branch, child, harmonics, harmonics_theta)
    # get distance vector
    dx = branch.center - child.center
    cartesian_2_spherical!(dx)
    regular_harmonic!(harmonics, harmonics_theta, dx..., tree.expansion_order[1])

    for j in 0:tree.expansion_order[1] # iterate over new Multipole coefficients B_j^k
        for k in 0:j
            i_jk = ((j * (j+1)) >> 1) + k + 1 # current index
            M = zeros(Complex{Float64}, 4) # vectorize later
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
    harmonics = Vector{Complex{Float64}}(undef, (tree.expansion_order[1]+1)^2)
    harmonics_theta = Vector{Complex{Float64}}(undef, (tree.expansion_order[1]+1)^2)

    # iterate over children
    for i_child in branch.first_branch:branch.first_branch + branch.n_branches - 1
        child = tree.branches[i_child]
        M2M!(tree, branch, child, harmonics, harmonics_theta)
    end
end

function M2L!(tree, elements, i_local, j_multipole)
    local_branch = tree.branches[i_local]
    multipole_branch = tree.branches[j_multipole]

    # preallocate
    harmonics = Vector{Complex{Float64}}(undef, (2*tree.expansion_order[1] + 1)^2)

    # get separation vector
    dx = local_branch.center - multipole_branch.center
    cartesian_2_spherical!(dx)
    irregular_harmonic!(harmonics, dx[1], dx[2], dx[3], 2*tree.expansion_order[1])
    for j in 0:tree.expansion_order[1]
        Cnm = odd_or_even(j)
        for k in 0:j
            jks = (j * (j + 1)) >> 1 + k + 1
            L = zeros(Complex{Float64}, 4)
            for n in 0:tree.expansion_order[1]
                for m in -n:-1
                    nms = (n * (n+1)) >> 1 - m + 1
                    jnkm = (j + n)^2 + j + n + m - k + 1
                    # jnkm_max = (P + P)^2 + P + P + -1 - 0 + 1 = (2P)^2 + 2P = 2P(2P+1)
                    for dim in 1:4
                        L[dim] += conj(multipole_branch.multipole_expansion[dim][nms]) * Cnm * harmonics[jnkm]
                    end
                end
                for m in 0:n
                    nms = (n * (n+1)) >> 1 + m + 1
                    jnkm = (j + n) * (j + n) + j + n + m - k + 1
                    # jnkm_max = 2P * 2P + 2P + P + P - 0 + 1 = (2P)^2 + 2P + 2P + 1 = 4P^2 + 4P + 1 = (2P + 1)^2
                    Cnm2 = Cnm * odd_or_even((k-m) * (1 >> (k>=m)) + m)
                    for dim in 1:4
                        L[dim] += multipole_branch.multipole_expansion[dim][nms] * Cnm2 * harmonics[jnkm]
                    end
                end
            end
            for dim in 1:4
                local_branch.local_expansion[dim][jks] += L[dim]
            end
        end
    end
end

function P2L!(tree, i_branch, source)
    branch = tree.branches[i_branch]
    irregular_harmonics = zeros(Complex{Float64},(tree.expansion_order[1]+1)^2)
    dx = cartesian_2_spherical(source.position - branch.center)
    irregular_harmonic!(irregular_harmonics, dx[1], dx[2], -dx[3], tree.expansion_order[1])
    q = source.strength
    for l in 0:tree.expansion_order[1]
        for m in 0:l
            i_abb = (l * (l+1)) >> 1 + m + 1
            i_exp = l^2 + l + m + 1
            for dim in 1:4
                branch.local_expansion[dim][i_abb] = irregular_harmonics[i_exp] * q[dim]
            end
        end
    end
end

function L2L!(tree, branch, child, harmonics, harmonics_theta)
    dx = child.center - branch.center
    cartesian_2_spherical!(dx)
    regular_harmonic!(harmonics, harmonics_theta, dx[1], dx[2], dx[3], tree.expansion_order[1])
    for j in 0:tree.expansion_order[1]
        for k in 0:j
            jks = (j * (j + 1)) >> 1 + k + 1
            L = zeros(Complex{Float64}, 4)
            for n in j:tree.expansion_order[1]
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
    harmonics = Vector{Complex{Float64}}(undef, (tree.expansion_order[1]+1)^2)
    harmonics_theta = Vector{Complex{Float64}}(undef, (tree.expansion_order[1]+1)^2)

    # iterate over children
    for i_child in branch.first_branch:branch.first_branch + branch.n_branches - 1
        child = tree.branches[i_child]
        L2L!(tree, branch, child, harmonics, harmonics_theta)
    end
end

"Calculates the potential at all child elements of a branch."
function L2P!(tree, elements, i_branch)
    branch = tree.branches[i_branch]
    harmonics = Vector{Complex{Float64}}(undef, (tree.expansion_order[1]+1)^2)
    harmonics_theta = Vector{Complex{Float64}}(undef, (tree.expansion_order[1]+1)^2)
    d_potential = zeros(4)
    for i_element in branch.first_element:branch.first_element + branch.n_elements - 1
        element = elements[tree.indices[i_element]]
        L2P!(element, tree, branch, harmonics, harmonics_theta, d_potential) # fix calling syntaxs
    end
end

function L2P!(element, tree, branch, harmonics, harmonics_theta, d_potential)
    dx = element.position - branch.center
    cartesian_2_spherical!(dx)
    regular_harmonic!(harmonics, harmonics_theta, dx[1], dx[2], dx[3], tree.expansion_order[1])
    for n in 0:tree.expansion_order[1]
        nm = n * n + n + 1 # m = 0
        nms = (n * (n+1)) >> 1 + 1 # m = 0
        for dim in 1:4
            d_potential[dim] = real(branch.local_expansion[dim][nms] * harmonics[nm])
        end
        element.potential .+= d_potential
        for m in 1:n
            nm = n * n + n + m + 1
            nms = (n * (n + 1)) >> 1 + m + 1
            for dim in 1:4
                d_potential[dim] = 2 * real(branch.local_expansion[dim][nms] * harmonics[nm])
            end
            element.potential .+= d_potential
        end
    end
end

# function L2P!(element::VectorSource, tree, branch, harmonics, harmonics_theta)
#     dx = fmm2.element.position - branch.center
#     spherical_1 = zeros(3)
#     spherical_2 = zeros(3)
#     spherical_3 = zeros(3)
#     cartesian = zeros(3)
#     cartesian_2_spherical!(dx)
#     regular_harmonic!(harmonics, harmonics_theta, dx..., tree.expansion_order[1]) # use the conjugate
#     d_potential = zeros(4)
#     for n in 0:tree.expansion_order[1]
#         nm = n * n + n + 1 # m = 0
#         nms = (n * (n+1)) >> 1 + 1 # m = 0
#         for dim in 1:4
#             d_potential[dim] = real(branch.local_expansion[dim][nms] * harmonics[nm])
#         end
#         update!(element, d_potential)
#         spherical_1[1] += real(branch.local_expansion[1][nms] * harmonics[nm]) / dx[1] * n
#         spherical_1[2] += real(branch.local_expansion[1][nms] * harmonics_theta[nm])
#         spherical_2[1] += real(branch.local_expansion[2][nms] * harmonics[nm]) / dx[1] * n
#         spherical_2[2] += real(branch.local_expansion[2][nms] * harmonics_theta[nm])
#         spherical_3[1] += real(branch.local_expansion[3][nms] * harmonics[nm]) / dx[1] * n
#         spherical_3[2] += real(branch.local_expansion[3][nms] * harmonics_theta[nm])
#         for m in 1:n
#             nm = n * n + n + m + 1
#             nms = (n * (n + 1)) >> 1 + m + 1
#             for dim in 1:4
#                 d_potential[dim] = 2 * real(branch.local_expansion[dim][nms] * harmonics[nm])
#             end
#             update!(element, d_potential)
#             spherical_1[1] += 2 * real(branch.local_expansion[1][nms] * harmonics[nm]) / dx[1] * n
#             spherical_1[2] += 2 * real(branch.local_expansion[1][nms] * harmonics_theta[nm])
#             spherical_1[3] += 2 * real(branch.local_expansion[1][nms] * harmonics[nm] * im) * m
#             spherical_2[1] += 2 * real(branch.local_expansion[2][nms] * harmonics[nm]) / dx[1] * n
#             spherical_2[2] += 2 * real(branch.local_expansion[2][nms] * harmonics_theta[nm])
#             spherical_2[3] += 2 * real(branch.local_expansion[2][nms] * harmonics[nm] * im) * m
#             spherical_3[1] += 2 * real(branch.local_expansion[3][nms] * harmonics[nm]) / dx[1] * n
#             spherical_3[2] += 2 * real(branch.local_expansion[3][nms] * harmonics_theta[nm])
#             spherical_3[3] += 2 * real(branch.local_expansion[3][nms] * harmonics[nm] * im) * m
#         end
#     end
#     spherical_2_cartesian!(cartesian, dx, spherical_1)
#     for ind in 1:3
#         element.Jexa[3*0 + ind] += cartesian[ind]
#     end
#     spherical_2_cartesian!(cartesian, dx, spherical_2)
#     for ind in 1:3
#         element.Jexa[3*1 + ind] += cartesian[ind]
#     end
#     spherical_2_cartesian!(cartesian, dx, spherical_3)
#     for ind in 1:3
#         element.Jexa[3*2 + ind] += cartesian[ind]
#     end
#   end
  