import FLOWFMM
fmm = FLOWFMM

#####
##### classic vortex particle method
#####

struct VortexParticles
    bodies
    potential
    velocity
    direct!
    B2M!
end

function direct_vortex!(target_potential, target_position, sources)
    n_targets = size(target_potential)[2]
    n_sources = size(sources)[2]
    for j_source in 1:n_sources
        x_source = sources[i_POSITION,j_source]
        for i_target in 1:n_targets
            target_jacobian = reshape(view(target_potential,i_POTENTIAL_JACOBIAN,i_target),3,4)
            x_target = target_position[:,i_target]
            dx = x_target - x_source
            r = sqrt(dx' * dx)
            if r > 0
                # calculate induced potential
                d_potential = sources[i_STRENGTH,j_source] / r
                target_potential[i_POTENTIAL,i_target] .+= d_potential

                # calculate induced jacobian
                # off diagonal elements are formed from the vector potential
                # diagonal elements are omitted (not needed?)
                gamma_over_R3 = -sources[i_STRENGTH,j_source] / r^3
                for j_potential in 1:4
                    for i_r in 1:3
                        target_jacobian[i_r,j_potential] += gamma_over_R3[j_potential] * dx[i_r]
                    end
                end
            end
        end
    end
end

function B2M_vortex!(tree, branch, bodies, n_bodies, harmonics)
    for i_body in 1:n_bodies
        dx = bodies[i_POSITION,i_body] - branch.center
        q = bodies[i_STRENGTH,i_body]
        fmm.cartesian_2_spherical!(dx)
        fmm.regular_harmonic!(harmonics, dx[1], dx[2], -dx[3], tree.expansion_order[1]) # Ylm^* -> -dx[3]
        # update values
        for l in 0:tree.expansion_order[1]
            for m in 0:l
                i_solid_harmonic = l^2 + l + m + 1
                i_compressed = 1 + (l * (l + 1)) >> 1 + m # only save half as Yl{-m} = conj(Ylm)
                for dim in 1:4
                    branch.multipole_expansion[dim][i_compressed] += harmonics[i_solid_harmonic] * q[dim]
                end
            end
        end
    end
end

function VortexParticles(position, strength;
    N = size(position)[2],
    potential = zeros(52,N),
    velocity = zeros(3,N),
)
    bodies = vcat(position, zeros(1,N), strength)
    return VortexParticles(bodies, potential, velocity, direct_vortex!, B2M_vortex!)
end

function VortexParticles(bodies;
    N = size(bodies)[2],
    potential = zeros(52,N),
    velocity = zeros(3,N),
)
    return VortexParticles(bodies, potential, velocity, direct_vortex!, B2M_vortex!)
end
