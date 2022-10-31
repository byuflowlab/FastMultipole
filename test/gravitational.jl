import FLOWFMM as fmm

#####
##### gravitational kernel and mass elements
#####
struct Mass
    position
    strength
    velocity
    potential
    J_potential # jacobian of the vector potential
    H_potential
end

function Mass(dims)
    Mass(zeros(dims), zeros(4), zeros(dims), zeros(4), zeros(3,4), zeros(3,4,3))
end

function P2M!(tree, branch, element::Mass, harmonics)
    dx = element.position - branch.center
    fmm.cartesian_2_spherical!(dx)
    fmm.regular_harmonic!(harmonics, dx[1], dx[2], -dx[3], tree.expansion_order[1]) # Ylm^* -> -dx[3]
    # update values
    for l in 0:tree.expansion_order[1]
        for m in 0:l
            i_solid_harmonic = l^2 + l + m + 1
            i_compressed = 1 + (l * (l + 1)) >> 1 + m # only save half as Yl{-m} = conj(Ylm)
            q = element.strength
            for dim in 1:4
                branch.multipole_expansion[dim][i_compressed] += harmonics[i_solid_harmonic] * q[dim]
            end
        end
    end
end

function P2P!(target::Mass, source::Mass)
    dx = target.position - source.position
    r = sqrt(dx' * dx)
    if r > 0
        dV = source.strength[1] / r
        target.potential[1] += dV
    end
end

function gravitational_potential_3D(Rho, m_source=1.0)
    rho_sq = Rho' * Rho
    rho = sqrt(rho_sq)
    return m_source / rho
end

function gravitational_dx(Rho, m_source=1.0)
    rho_sq = Rho' * Rho
    dVdx = -m_source * Rho[1] / rho_sq^1.5
end

function gravitational_dy(Rho, m_source=1.0)
    rho_sq = Rho' * Rho
    dVdx = -m_source * Rho[2] / rho_sq^1.5
end

function gravitational_dz(Rho, m_source=1.0)
    rho_sq = Rho' * Rho
    dVdx = -m_source * Rho[3] / rho_sq^1.5
end

function gravitational_dx2(Rho, m_source=1.0)
    rho_sq = Rho' * Rho
    d2Vdx2 = -m_source * (-2 * Rho[1]^2 + Rho[2]^2 + Rho[3]^2) / rho_sq^2.5
end

function gravitational_dy2(Rho, m_source=1.0)
    rho_sq = Rho' * Rho
    d2Vdy2 = -m_source * (Rho[1]^2 - 2*Rho[2]^2 + Rho[3]^2) / rho_sq^2.5
end

function gravitational_dz2(Rho, m_source=1.0)
    rho_sq = Rho' * Rho
    d2Vdz2 = -m_source * (Rho[1]^2 + Rho[2]^2 - 2*Rho[3]^2) / rho_sq^2.5
end

function gravitational_dxdy(Rho, m_source=1.0)
    rho_sq = Rho' * Rho
    d2Vdxdy = 3 * m_source * Rho[1] * Rho[2] / rho_sq^2.5
end

function gravitational_dydz(Rho, m_source=1.0)
    rho_sq = Rho' * Rho
    d2Vdydz = 3 * m_source * Rho[2] * Rho[3] / rho_sq^2.5
end

function gravitational_dxdz(Rho, m_source=1.0)
    rho_sq = Rho' * Rho
    d2Vdxdz = 3 * m_source * Rho[3] * Rho[1] / rho_sq^2.5
end

function gravitational_dx3(Rho, m_source=1.0)
    rho_sq = Rho' * Rho
    d3Vdx3 = -3 * m_source * Rho[1] * (2 * Rho[1]^2 - 3 * (Rho[2]^2 + Rho[3]^2)) / rho_sq^3.5
end

function gravitational_dy3(Rho, m_source=1.0)
    rho_sq = Rho' * Rho
    d3Vdy3 = -3 * m_source * Rho[2] * (2 * Rho[2]^2 - 3 * (Rho[3]^2 + Rho[1]^2)) / rho_sq^3.5
end

function gravitational_dz3(Rho, m_source=1.0)
    rho_sq = Rho' * Rho
    d3Vdz3 = -3 * m_source * Rho[3] * (2 * Rho[3]^2 - 3 * (Rho[1]^2 + Rho[2]^2)) / rho_sq^3.5
end

function gravitational_dx2dy(Rho, m_source=1.0)
    rho_sq = Rho' * Rho
    d3Vdx2dy = 3 * m_source * Rho[2] * (-4 * Rho[1]^2 + Rho[2]^2 + Rho[3]^2) / rho_sq^3.5
end

function gravitational_dx2dz(Rho, m_source=1.0)
    rho_sq = Rho' * Rho
    d3Vdx2dy = 3 * m_source * Rho[3] * (-4 * Rho[1]^2 + Rho[2]^2 + Rho[3]^2) / rho_sq^3.5
end

function gravitational_dy2dz(Rho, m_source=1.0)
    rho_sq = Rho' * Rho
    d3Vdy2dz = 3 * m_source * Rho[3] * (-4 * Rho[2]^2 + Rho[3]^2 + Rho[1]^2) / rho_sq^3.5
end

function gravitational_dxdy2(Rho, m_source=1.0)
    rho_sq = Rho' * Rho
    d3Vdy2dz = 3 * m_source * Rho[1] * (-4 * Rho[2]^2 + Rho[3]^2 + Rho[1]^2) / rho_sq^3.5
end

function gravitational_dxdz2(Rho, m_source=1.0)
    rho_sq = Rho' * Rho
    d3Vdz2dx = 3 * m_source * Rho[1] * (-4 * Rho[3]^2 + Rho[1]^2 + Rho[2]^2) / rho_sq^3.5
end

function gravitational_dydz2(Rho, m_source=1.0)
    rho_sq = Rho' * Rho
    d3Vdz2dx = 3 * m_source * Rho[2] * (-4 * Rho[3]^2 + Rho[1]^2 + Rho[2]^2) / rho_sq^3.5
end

function gravitational_dxdydz(Rho, m_source=1.0)
    rho_sq = Rho' * Rho
    d3Vdxdydz = -15 * m_source * Rho[1] * Rho[2] * Rho[3] / rho_sq^3.5
end

function gravitational_dx4(Rho, m_source=1.0)
    rho_sq = Rho' * Rho
    dx4 = 3 * m_source * (8 * Rho[1]^4 - 24 * Rho[1]^2 * (Rho[2]^2 + Rho[3]^2) + 3 * (Rho[2]^2 + Rho[3]^2)^2) / rho_sq^4.5
end

function gravitational_dy4(Rho, m_source=1.0)
    rho_sq = Rho' * Rho
    dy4 = 3 * m_source * (8 * Rho[2]^4 - 24 * Rho[2]^2 * (Rho[3]^2 + Rho[1]^2) + 3 * (Rho[3]^2 + Rho[1]^2)^2) / rho_sq^4.5
end

function gravitational_dz4(Rho, m_source=1.0)
    rho_sq = Rho' * Rho
    dz4 = 3 * m_source * (8 * Rho[3]^4 - 24 * Rho[3]^2 * (Rho[1]^2 + Rho[2]^2) + 3 * (Rho[1]^2 + Rho[2]^2)^2) / rho_sq^4.5
end

function gravitational_dx3dy(Rho, m_source=1.0)
    rho_sq = Rho' * Rho
    dx3dy = -15 * m_source * Rho[2] * (3 * Rho[1] * (Rho[2]^2 + Rho[3]^2) - 4 * Rho[1]^3) / rho_sq^4.5
end

function gravitational_dx3dz(Rho, m_source=1.0)
    rho_sq = Rho' * Rho
    dx3dz = -15 * m_source * Rho[3] * (3 * Rho[1] * (Rho[3]^2 + Rho[2]^2) - 4 * Rho[1]^3) / rho_sq^4.5
end

function gravitational_dy3dz(Rho, m_source=1.0)
    rho_sq = Rho' * Rho
    dy3dz = -15 * m_source * Rho[3] * (3 * Rho[2] * (Rho[3]^2 + Rho[1]^2) - 4 * Rho[2]^3) / rho_sq^4.5
end

function gravitational_dxdy3(Rho, m_source=1.0)
    rho_sq = Rho' * Rho
    dxdy3 = -15 * m_source * Rho[1] * (3 * Rho[2] * (Rho[1]^2 + Rho[3]^2) - 4 * Rho[2]^3) / rho_sq^4.5
end

function gravitational_dxdz3(Rho, m_source=1.0)
    rho_sq = Rho' * Rho
    dxdz3 = -15 * m_source * Rho[1] * (3 * Rho[3] * (Rho[1]^2 + Rho[2]^2) - 4 * Rho[3]^3) / rho_sq^4.5
end

function gravitational_dydz3(Rho, m_source=1.0)
    rho_sq = Rho' * Rho
    dydz3 = -15 * m_source * Rho[2] * (3 * Rho[3] * (Rho[2]^2 + Rho[1]^2) - 4 * Rho[3]^3) / rho_sq^4.5
end

function gravitational_dx2dy2(Rho, m_source=1.0)
    rho_sq = Rho' * Rho
    dx2dy2 = -3 * m_source * (4 * Rho[1]^4 + 3 * Rho[1]^2 * (Rho[3]^2 - 9 * Rho[2]^2) +
             4 * Rho[2]^4 + 3 * Rho[2]^2 * Rho[3]^2 - Rho[3]^4) / rho_sq^4.5
end

function gravitational_dy2dz2(Rho, m_source=1.0)
    rho_sq = Rho' * Rho
    dy2dz2 = -3 * m_source * (4 * Rho[3]^4 + 3 * Rho[3]^2 * (Rho[1]^2 - 9 * Rho[2]^2) +
            4 * Rho[2]^4 + 3 * Rho[2]^2 * Rho[1]^2 - Rho[1]^4) / rho_sq^4.5
end

function gravitational_dx2dz2(Rho, m_source=1.0)
    rho_sq = Rho' * Rho
    dx2dz2 = -3 * m_source * (4 * Rho[1]^4 + 3 * Rho[1]^2 * (Rho[2]^2 - 9 * Rho[3]^2) +
            4 * Rho[3]^4 + 3 * Rho[3]^2 * Rho[2]^2 - Rho[2]^4) / rho_sq^4.5
end

function gravitational_dx2dydz(Rho, m_source=1.0)
    rho_sq = Rho' * Rho
    dx2dydz = -15 * m_source * Rho[2] * Rho[3] * (-6 * Rho[1]^2 + Rho[2]^2 + Rho[3]^2) / rho_sq^4.5
end

function gravitational_dxdy2dz(Rho, m_source=1.0)
    rho_sq = Rho' * Rho
    dxdy2dz = -15 * m_source * Rho[1] * Rho[3] * (-6 * Rho[2]^2 + Rho[1]^2 + Rho[3]^2) / rho_sq^4.5
end

function gravitational_dxdydz2(Rho, m_source=1.0)
    rho_sq = Rho' * Rho
    dxdydz2 = -15 * m_source * Rho[2] * Rho[1] * (-6 * Rho[3]^2 + Rho[2]^2 + Rho[1]^2) / rho_sq^4.5
end

order = 4
dims = 3
derivatives = Vector{Function}(undef, fmm.n_terms(order, dims)) # default through 4th derivatives
# raw potential
derivatives[1] = gravitational_potential_3D
# derivatives[1,1,1] = gravitational_potential_3D
# first derivatives
derivatives[2] = gravitational_dx
# derivatives[2,1,1] = gravitational_dx
derivatives[3] = gravitational_dy
# derivatives[1,2,1] = gravitational_dy
derivatives[4] = gravitational_dz
# derivatives[1,1,2] = gravitational_dz
# second derivatives
derivatives[5] = gravitational_dx2
# derivatives[3,1,1] = gravitational_dx2
derivatives[6] = gravitational_dxdy
# derivatives[2,2,1] = gravitational_dxdy
derivatives[7] = gravitational_dxdz
# derivatives[2,1,2] = gravitational_dxdz
derivatives[8] = gravitational_dy2
# derivatives[1,3,1] = gravitational_dy2
derivatives[9] = gravitational_dydz
# derivatives[1,2,2] = gravitational_dydz
derivatives[10] = gravitational_dz2
# derivatives[1,1,3] = gravitational_dz2
# third derivatives
derivatives[11] = gravitational_dx3
# derivatives[4,1,1] = gravitational_dx3
derivatives[12] = gravitational_dx2dy
# derivatives[3,2,1] = gravitational_dx2dy
derivatives[13] = gravitational_dx2dz
# derivatives[3,1,2] = gravitational_dx2dz
derivatives[14] = gravitational_dxdy2
# derivatives[2,3,1] = gravitational_dxdy2
derivatives[15] = gravitational_dxdydz
# derivatives[2,2,2] = gravitational_dxdydz
derivatives[16] = gravitational_dxdz2
# derivatives[2,1,3] = gravitational_dxdz2
derivatives[17] = gravitational_dy3
# derivatives[1,4,1] = gravitational_dy3
derivatives[18] = gravitational_dy2dz
# derivatives[1,3,2] = gravitational_dy2dz
derivatives[19] = gravitational_dydz2
# derivatives[1,2,3] = gravitational_dydz2
derivatives[20] = gravitational_dz3
# derivatives[1,1,4] = gravitational_dz3
# fourth derivatives
derivatives[21] = gravitational_dx4
# derivatives[5,1,1] = gravitational_dx4
derivatives[22] = gravitational_dx3dy
# derivatives[4,2,1] = gravitational_dx3dy
derivatives[23] = gravitational_dx3dz
# derivatives[4,1,2] = gravitational_dx3dz
derivatives[24] = gravitational_dx2dy2
# derivatives[3,3,1] = gravitational_dx2dy2
derivatives[25] = gravitational_dx2dydz
# derivatives[3,2,2] = gravitational_dx2dydz
derivatives[26] = gravitational_dx2dz2
# derivatives[3,1,3] = gravitational_dx2dz2
derivatives[27] = gravitational_dxdy3
# derivatives[2,4,1] = gravitational_dxdy3
derivatives[28] = gravitational_dxdy2dz
# derivatives[2,3,2] = gravitational_dxdy2dz
derivatives[29] = gravitational_dxdydz2
# derivatives[2,2,3] = gravitational_dxdydz2
derivatives[30] = gravitational_dxdz3
# derivatives[2,1,4] = gravitational_dxdz3
derivatives[31] = gravitational_dy4
# derivatives[1,5,1] = gravitational_dy4
derivatives[32] = gravitational_dy3dz
# derivatives[1,4,2] = gravitational_dy3dz
derivatives[33] = gravitational_dy2dz2
# derivatives[1,3,3] = gravitational_dy2dz2
derivatives[34] = gravitational_dydz3
# derivatives[1,2,4] = gravitational_dydz3
derivatives[35] = gravitational_dz4
# derivatives[1,1,5] = gravitational_dz4