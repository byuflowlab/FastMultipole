import FLOWFMM as fmm
using WriteVTK
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
    index
    potential
    velocity
    direct!
    B2M!
end

function B2M_gravitational!(branch, bodies, harmonics, expansion_order)
    for i_body in 1:size(bodies)[2]
        dx = bodies[i_POSITION,i_body] - branch.center
        q = bodies[i_STRENGTH,i_body]
        fmm.cartesian_2_spherical!(dx)
        fmm.regular_harmonic!(harmonics, dx[1], dx[2], -dx[3], expansion_order) # Ylm^* -> -dx[3]
        # update values
        for l in 0:expansion_order
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
    # n_targets = length(target_potential.indices[2])
    n_targets = size(target_potential)[2]
    # n_sources = length(sources.indices[2])
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
    target_potential = view(gravitational.potential,:,:)
    target_position = view(gravitational.bodies,1:3,:)
    sources = view(gravitational.bodies,:,:)
    direct_gravitational!(target_potential, target_position, sources)
end

function Gravitational(bodies)
    N = size(bodies)[2]
    index = zeros(Int32,N)
    potential = zeros(52,N)
    velocity = zeros(3,N)
    return Gravitational(bodies, index, potential, velocity, direct_gravitational!, B2M_gravitational!)
end

function save_vtk(filename, element::Gravitational, nt=0; compress=false, extra_fields=nothing)
    _, n = size(element.bodies)
    WriteVTK.vtk_grid(filename*"_point_masses."*string(nt)*".vts", reshape(view(element.bodies,1:3,:),3,n,1,1); compress) do vtk
        vtk["strength"] = reshape(view(element.bodies,4,:), 1, n, 1, 1)
        vtk["velocity"] = reshape(element.velocity, 3, n, 1, 1)
        vtk["scalar potential"] = reshape(view(element.potential,1,:), n, 1, 1)
        vtk["vector potential"] = reshape(view(element.potential,2:4,:), 3, n, 1, 1)
        if !isnothing(extra_fields)
            for i in 1:length(extra_fields)
                vtk[extra_fields[i][1]] = extra_fields[i][2]
            end
        end
    end
end