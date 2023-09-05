# module VPM

import FLOWFMM as fmm
import WriteVTK
using SpecialFunctions:erf

#####
##### classic vortex particle method
#####
const ONE_OVER_4PI = 1/4/pi
const sqrt2 = sqrt(2)
const i_POSITION_vortex = 1:3
const i_STRENGTH_vortex = 4:6
const i_POTENTIAL_SCALAR = 1:1
const i_POTENTIAL_VECTOR = 2:4
i_POTENTIAL_JACOBIAN = 5:16
i_POTENTIAL_HESSIAN = 17:52
const i_VELOCITY_vortex = 1:3
const i_STRETCHING_vortex = 4:6

struct Vorton{TF}
    position::SVector{3,TF}
    strength::SVector{3,TF}
    sigma::TF
end

struct VortexParticles{TF}
    bodies::Vector{Vorton{TF}}
    potential::Matrix{TF}
    velocity_stretching::Matrix{TF}
end

abstract type IntegrationScheme end

struct Euler{TF} <: IntegrationScheme
    dt::TF
end

"""
Classical formulation so far.
"""
function fmm.direct!(target_system, target_index, source_system::VortexParticles, source_index)
    jacobian = zeros(eltype(target_system),3,4)
    hessian = zeros(eltype(target_system),3,3,4)
    for j_source in source_index
        x_source = source_system[j_source,fmm.POSITION]
        for i_target in target_index
            jacobian .*= 0.0
            hessian .*= 0.0
            x_target = target_system[i_target,fmm.POSITION]
            dx = x_target - x_source
            r = sqrt(dx' * dx)
            if r > 0
                # calculate induced potential
                gamma_over_R = source_system.bodies[j_source].strength / r
                target_potential[i_target,fmm.VECTOR_POTENTIAL] .+= gamma_over_R

                # calculate induced jacobian
                # off diagonal elements are formed from the vector potential
                # diagonal elements are omitted (not needed?)
                # gamma_over_R3 = source_bodies[i_STRENGTH_vortex,j_source] / r^3
                gamma_over_R /= r^2
                for j_potential in 2:4
                    for i_r in 1:3
                        jacobian[i_r,j_potential] -= gamma_over_R[j_potential-1] * dx[i_r]
                    end
                end
                gamma_over_R /= r^2
                hessian[1,1,i_POTENTIAL_VECTOR] += gamma_over_R .* (2 * dx[1]^2 - dx[2]^2 - dx[3]^2) # dx2
                hessian[1,2,i_POTENTIAL_VECTOR] += gamma_over_R * 3 * dx[1] * dx[2] # dxdy
                hessian[1,3,i_POTENTIAL_VECTOR] += gamma_over_R * 3 * dx[1] * dx[3] # dxdz
                hessian[2,1,i_POTENTIAL_VECTOR] += gamma_over_R * 3 * dx[1] * dx[2] # dxdy
                hessian[2,2,i_POTENTIAL_VECTOR] += gamma_over_R .* (-dx[1]^2 + 2 * dx[2]^2 - dx[3]^2) # dy2
                hessian[2,3,i_POTENTIAL_VECTOR] += gamma_over_R * 3 * dx[2] * dx[3] # dydz
                hessian[3,1,i_POTENTIAL_VECTOR] += gamma_over_R * 3 * dx[1] * dx[3] # dxdz
                hessian[3,2,i_POTENTIAL_VECTOR] += gamma_over_R * 3 * dx[2] * dx[3] # dydz
                hessian[3,3,i_POTENTIAL_VECTOR] += gamma_over_R .* (-dx[1]^2 - dx[2]^2 + 2 * dx[3]^2) # dz2
            end
            target_system[i_target,fmm.JACOBIAN] .+= jacobian
            target_system[i_target,fmm.HESSIAN] .+= hessian
        end
    end
end

"""
Assumes singular particles.
"""
function fmm.B2M!(branch, system::VortexParticles, index, harmonics, expansion_order)
    for i_body in index
        dx = system[i_body,fmm.POSITION] - branch.center
        q = system.bodies[i_body].strength
        fmm.cartesian_2_spherical!(dx)
        fmm.regular_harmonic!(harmonics, dx[1], dx[2], -dx[3], expansion_order) # Ylm^* -> -dx[3]
        # update values
        for l in 0:expansion_order
            for m in 0:l
                i_solid_harmonic = l^2 + l + m + 1
                i_compressed = 1 + (l * (l + 1)) >> 1 + m # only save half as Yl{-m} = conj(Ylm)
                for dim in 2:4 # does not affect scalar potential
                    branch.multipole_expansion[dim][i_compressed] += harmonics[i_solid_harmonic] * q[dim-1]
                end
            end
        end
    end
end

function VortexParticles(position, strength;
    N = size(position)[2],
    potential = zeros(52,N),
    velocity_stretching = zeros(3+3,N),
)
    @assert size(position)[1] == 3
    @assert size(strength)[1] == 3
    bodies = [Vorton(SVector{3}(position[:,i]), SVector{3}(strength[:,i]), 0.0) for i in 1:size(position)[2]]
    return VortexParticles(bodies, potential, velocity_stretching)
end

function VortexParticles(bodies;
    potential = zeros(52,N),
    velocity_stretching = zeros(3+3,N),
)
    bodies = [Vorton(SVector{3}(bodies[1:3,i]), SVector{3}(bodies[4:6,i])) for i in 1:size(bodies)[2]]
    return VortexParticles(bodies, potential, velocity_stretching)
end

@inline function update_velocity_stretching!(system, i_body)
    # velocity is the curl of the vector potential
    # minus the gradient of the scalar potential
    jacobian = system[i_body,fmm.JACOBIAN]
    system.velocity_stretching[1,i_body] = (-jacobian[1,1] + jacobian[2,4] - jacobian[3,3]) * ONE_OVER_4PI
    system.velocity_stretching[2,i_body] = (-jacobian[2,1] + jacobian[3,2] - jacobian[1,4]) * ONE_OVER_4PI
    system.velocity_stretching[3,i_body] = (-jacobian[3,1] + jacobian[1,3] - jacobian[2,2]) * ONE_OVER_4PI
    # jacobian of the velocity
    hessian = system[i_body,fmm.HESSIAN]
    duidxj = fmm.SMatrix{3,3}([
        -hessian[1,1,1]+hessian[2,1,4]-hessian[3,1,3] -hessian[1,2,1]+hessian[2,2,4]-hessian[3,2,3] -hessian[1,3,1]+hessian[2,3,4]-hessian[3,3,3];
        -hessian[2,1,1]+hessian[3,1,2]-hessian[1,1,4] -hessian[2,2,1]+hessian[3,2,2]-hessian[1,2,4] -hessian[2,3,1]+hessian[3,3,2]-hessian[1,3,4];
        -hessian[3,1,1]+hessian[1,1,3]-hessian[2,1,2] -hessian[3,2,1]+hessian[1,2,3]-hessian[2,2,2] -hessian[3,3,1]+hessian[1,3,3]-hessian[2,3,2];
    ])
    # vorticity = @SVector [
    #     jacobian[2,3] - jacobian[3,2],
    #     jacobian[3,1] - jacobian[1,3],
    #     jacobian[1,2] - jacobian[2,1]
    # ]
    # stretching term (omega dot nabla)

    # total derivative of the strength due to vortex stretching
    fmm.mul!(view(system.velocity_stretching,i_STRETCHING_vortex,i_body), duidxj, system.bodies[i_body].strength * ONE_OVER_4PI)

    return nothing
end

function update_velocity_stretching!(vortex_particles::VortexParticles)
    for i_body in 1:size(bodies)[2]
        update_velocity_stretching!(vortex_particles, i_body)
    end
end

"Assumes velocities and potential and derivatives are already calculated."
function (euler::Euler)(vortex_particles::VortexParticles, fmm_options, direct)
    # reset potential
    vortex_particles.potential .*= 0
    # reset velocity and stretching
    vortex_particles.velocity_stretching .*= 0

    # calculate influences
    if direct
        fmm.direct!(vortex_particles)
    else
        fmm.fmm!((vortex_particles,), fmm_options)
    end

    # convect bodies
    bodies = vortex_particles.bodies
    potential = vortex_particles.potential
    velocity_stretching = vortex_particles.velocity_stretching
    update_velocity_stretching!(vortex_particles)
    for i_body in 1:length(vortex_particles)
        vortex_particles[i_body, fmm.POSITION] .+= vortex_particles.velocity_stretching[1:3,i_body] * euler.dt
        vortex_particles.bodies[i_body] .+= vortex_particles.velocity_stretching[4:6] * euler.dt
    end
end

function convect!(vortex_particles::VortexParticles, nsteps;
        # integration options
        integrate!::IntegrationScheme=Euler(1.0),
        # fmm options
        fmm_p=4, fmm_ncrit=50, fmm_theta=4.0, fmm_targets=SVector{1}(Int8(1)),
        direct::Bool=false,
        # save options
        save::Bool=true, filename::String="default", compress::Bool=false,
    )
    fmm_options = fmm.Options(fmm_p, fmm_ncrit, fmm_theta, fmm_targets)
    save && save_vtk(filename, vortex_particles; compress)
    for istep in 1:nsteps
        integrate!(vortex_particles, fmm_options, direct)
        save && save_vtk(filename, vortex_particles, istep; compress)
    end
    return nothing
end

function save_vtk(filename, vortex_particles::VortexParticles, nt=0; compress=false)
    n_bodies = size(vortex_particles.bodies)[2]
    WriteVTK.vtk_grid(filename*"."*string(nt)*".vts", reshape(vortex_particles.bodies[i_POSITION_vortex,:], 3, n_bodies,1,1); compress) do vtk
        vtk["vector strength"] = reshape(vortex_particles.bodies[i_STRENGTH_vortex,:],3,n_bodies,1,1)
        vtk["scalar potential"] = reshape(vortex_particles.potential[i_POTENTIAL_SCALAR,:],n_bodies,1,1)
        vtk["vector potential"] = reshape(vortex_particles.potential[i_POTENTIAL_VECTOR,:],3,n_bodies,1,1)
        vtk["velocity"] = reshape(vortex_particles.velocity_stretching[i_VELOCITY_vortex,:],3,n_bodies,1,1)
        vtk["stretching"] = reshape(vortex_particles.velocity_stretching[i_STRETCHING_vortex,:],3,n_bodies,1,1)
    end
end
