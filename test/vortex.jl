import FLOWFMM
fmm = FLOWFMM

#####
##### classic vortex particle method
#####

struct Vorton
    position
    strength
    velocity
    potential
    J_potential # jacobian of the vector potential
    H_potential # hessian of the vector potential
end

function Vorton(position, strength;
    velocity = zeros(3),
    potential = zeros(4),
    J_potential = zeros(3,4),
    H_potential = zeros(3,3,3)
)
    Vorton(position, strength, velocity, potential, J_potential, H_potential)
end

function P2P!(source::Vorton, target::Vorton)
    dx = target.position - source.position
    r = sqrt(dx' * dx)
    if r > 0
        # calculate induced potential
        d_potential = source.strength / r
        target.potential .+= d_potential

        # calculate induced jacobian
        # off diagonal elements are formed from the vector potential
        # diagonal elements are omitted (not needed?)
        gamma_over_R3 = -source.strength / r^3
        for i_r in 1:3
            for j_potential in 1:4
                target.J_potential[i_r,j_potential] += gamma_over_R3[j_potential] * dx[i_r]
            end
        end
    end
end

function P2M!(tree, branch, element::Vorton, harmonics)
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
