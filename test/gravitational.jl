import FLOWFMM as fmm
i_POSITION = fmm.i_POSITION
i_STRENGTH = 4:7
i_POTENTIAL = fmm.i_POTENTIAL
i_POTENTIAL_JACOBIAN = fmm.i_POTENTIAL_JACOBIAN
i_POTENTIAL_HESSIAN = fmm.i_POTENTIAL_HESSIAN

#####
##### gravitational kernel and mass elements
#####
struct Gravitational
    bodies
    potential
    velocity
    direct!
    B2M!
end

function B2M_gravitational!(tree, branch, bodies, n_bodies, harmonics)
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

function direct_gravitational!(target_potential, target_position, sources)
    n_targets = size(target_potential)[2]
    n_sources = size(sources)[2]
    for j_source in 1:n_sources
        x_source = sources[i_POSITION,j_source]
        for i_target in 1:n_targets
            x_target = target_position[:,i_target]
            dx = x_target - x_source
            r = sqrt(dx'*dx)
            if r > 0
                dV = sources[i_STRENGTH[1],j_source] / r
                target_potential[i_POTENTIAL[1],i_target] += dV
            end
        end
    end
end

function direct_gravitational!(gravitational::Gravitational)
    target_potential = gravitational.potential
    target_position = view(gravitational.bodies,1:3,:)
    sources = gravitational.bodies
    direct_gravitational!(target_potential, target_position, sources)
end

function Gravitational(bodies)
    N = size(bodies)[2]
    potential = zeros(52,N)
    velocity = zeros(3,N)
    return Gravitational(bodies, potential, velocity, direct_gravitational!, B2M_gravitational!)
end
