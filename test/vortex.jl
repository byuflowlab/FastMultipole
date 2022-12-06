# module VPM

import FLOWFMM as fmm
import WriteVTK
using SpecialFunctions:erf

#####
##### classic vortex particle method
#####
const ONE_OVER_4PI = 1/4/pi
const sqrt2 = sqrt(2)
const i_POSITION_vortex = fmm.i_POSITION
const i_STRENGTH_vortex = 4:6
const i_POTENTIAL_SCALAR = fmm.i_POTENTIAL[1] # 1:4
const i_POTENTIAL_VECTOR = fmm.i_POTENTIAL[2:4]
i_POTENTIAL_JACOBIAN = fmm.i_POTENTIAL_JACOBIAN # 5:16
i_POTENTIAL_HESSIAN = fmm.i_POTENTIAL_HESSIAN # 17:52
const i_VELOCITY_vortex = 1:3
const i_STRETCHING_vortex = 4:6

struct VortexParticles
    bodies
    index
    potential
    velocity_stretching # velocity and stretching term
    direct!
    B2M!
end

abstract type IntegrationScheme end

struct Euler{TF} <: IntegrationScheme
    dt::TF
end

"""
Classical formulation so far.
"""
function direct_vortex!(target_potential, target_positions, source_bodies)
    n_targets = size(target_potential)[2]
    n_sources = size(source_bodies)[2]
    for j_source in 1:n_sources
        x_source = source_bodies[i_POSITION_vortex,j_source]
        for i_target in 1:n_targets
            target_jacobian = reshape(view(target_potential,i_POTENTIAL_JACOBIAN,i_target),3,4)
            target_hessian = reshape(view(target_potential,i_POTENTIAL_HESSIAN,i_target),3,3,4)
            x_target = target_positions[:,i_target]
            dx = x_target - x_source
            r = sqrt(dx' * dx)
            if r > 0
                # calculate induced potential
                gamma_over_R = source_bodies[i_STRENGTH_vortex,j_source] / r
                target_potential[i_POTENTIAL_VECTOR,i_target] .+= gamma_over_R

                # calculate induced jacobian
                # off diagonal elements are formed from the vector potential
                # diagonal elements are omitted (not needed?)
                # gamma_over_R3 = source_bodies[i_STRENGTH_vortex,j_source] / r^3
                gamma_over_R /= r^2
                for j_potential in 2:4
                    for i_r in 1:3
                        target_jacobian[i_r,j_potential] -= gamma_over_R[j_potential-1] * dx[i_r]
                    end
                end
                gamma_over_R /= r^2
                target_hessian[1,1,i_POTENTIAL_VECTOR] += gamma_over_R .* (2 * dx[1]^2 - dx[2]^2 - dx[3]^2) # dx2
                target_hessian[1,2,i_POTENTIAL_VECTOR] += gamma_over_R * 3 * dx[1] * dx[2] # dxdy
                target_hessian[1,3,i_POTENTIAL_VECTOR] += gamma_over_R * 3 * dx[1] * dx[3] # dxdz
                target_hessian[2,1,i_POTENTIAL_VECTOR] += gamma_over_R * 3 * dx[1] * dx[2] # dxdy
                target_hessian[2,2,i_POTENTIAL_VECTOR] += gamma_over_R .* (-dx[1]^2 + 2 * dx[2]^2 - dx[3]^2) # dy2
                target_hessian[2,3,i_POTENTIAL_VECTOR] += gamma_over_R * 3 * dx[2] * dx[3] # dydz
                target_hessian[3,1,i_POTENTIAL_VECTOR] += gamma_over_R * 3 * dx[1] * dx[3] # dxdz
                target_hessian[3,2,i_POTENTIAL_VECTOR] += gamma_over_R * 3 * dx[2] * dx[3] # dydz
                target_hessian[3,3,i_POTENTIAL_VECTOR] += gamma_over_R .* (-dx[1]^2 - dx[2]^2 + 2 * dx[3]^2) # dz2
            end
        end
    end
end

"""
Gaussian smoothing.
"""
function direct_vortex_gaussian!(target_potential, target_positions, source_bodies)
    n_targets = size(target_potential)[2]
    n_sources = size(source_bodies)[2]
    for j_source in 1:n_sources
        x_source = source_bodies[i_POSITION_vortex,j_source]
        for i_target in 1:n_targets
            target_jacobian = reshape(view(target_potential,i_POTENTIAL_JACOBIAN,i_target),3,4)
            target_hessian = reshape(view(target_potential,i_POTENTIAL_HESSIAN,i_target),3,3,4)
            x_target = target_positions[:,i_target]
            dx = x_target - x_source
            r2 = dx' * dx
            if r2 > 0
                r = sqrt(r2)
                # calculate induced potential
                erfrs2 = erf(r/sqrt2)
                sigma = source_bodies[i_SIGMA,j_source]
                gamma_over_R = source_bodies[i_STRENGTH_vortex,j_source] / r * erfrs2 / sigma
                target_potential[i_POTENTIAL_VECTOR,i_target] .+= gamma_over_R

                # calculate induced jacobian
                # off diagonal elements are formed from the vector potential
                # diagonal elements are omitted (not needed?)
                # gamma_over_R3 = source_bodies[i_STRENGTH_vortex,j_source] / r^3
                dGdr = -erfrs2/r2 + sqrt2/sqrt(pi)/r*exp(-r2/2)
                # gamma_over_R /= r^2
                # gamma_over_R *=
                for j_potential in 2:4
                    for i_r in 1:3
                        target_jacobian[i_r,j_potential] -= gamma_over_R[j_potential-1] * dx[i_r]
                    end
                end
                gamma_over_R /= r^2
                target_hessian[1,1,i_POTENTIAL_VECTOR] += gamma_over_R .* (2 * dx[1]^2 - dx[2]^2 - dx[3]^2) # dx2
                target_hessian[1,2,i_POTENTIAL_VECTOR] += gamma_over_R * 3 * dx[1] * dx[2] # dxdy
                target_hessian[1,3,i_POTENTIAL_VECTOR] += gamma_over_R * 3 * dx[1] * dx[3] # dxdz
                target_hessian[2,1,i_POTENTIAL_VECTOR] += gamma_over_R * 3 * dx[1] * dx[2] # dxdy
                target_hessian[2,2,i_POTENTIAL_VECTOR] += gamma_over_R .* (-dx[1]^2 + 2 * dx[2]^2 - dx[3]^2) # dy2
                target_hessian[2,3,i_POTENTIAL_VECTOR] += gamma_over_R * 3 * dx[2] * dx[3] # dydz
                target_hessian[3,1,i_POTENTIAL_VECTOR] += gamma_over_R * 3 * dx[1] * dx[3] # dxdz
                target_hessian[3,2,i_POTENTIAL_VECTOR] += gamma_over_R * 3 * dx[2] * dx[3] # dydz
                target_hessian[3,3,i_POTENTIAL_VECTOR] += gamma_over_R .* (-dx[1]^2 - dx[2]^2 + 2 * dx[3]^2) # dz2
            end
        end
    end
end

"""
Assumes singular particles.
"""
function B2M_vortex!(tree, branch, bodies, n_bodies, harmonics)
    for i_body in 1:n_bodies
        dx = bodies[i_POSITION_vortex,i_body] - branch.center
        q = bodies[i_STRENGTH_vortex,i_body]
        fmm.cartesian_2_spherical!(dx)
        fmm.regular_harmonic!(harmonics, dx[1], dx[2], -dx[3], tree.expansion_order) # Ylm^* -> -dx[3]
        # update values
        for l in 0:tree.expansion_order
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
    bodies = vcat(position, strength)
    return VortexParticles(bodies; N, potential, velocity_stretching)
end

function VortexParticles(bodies;
    N = size(bodies)[2],
    index = zeros(Int32,N),
    potential = zeros(52,N),
    velocity_stretching = zeros(3+3,N),
)
    @assert size(bodies)[1] == 6 "bodies size first index is incorrect; got $(size(bodies)[1]); expected 6"
    return VortexParticles(bodies, index, potential, velocity_stretching, direct_vortex!, B2M_vortex!)
end

@inline function update_velocity_stretching!(velocity, body, potential)
    # velocity is the curl of the vector potential
    # minus the gradient of the scalar potential
    jacobian = reshape(potential[i_POTENTIAL_JACOBIAN],3,4)
    velocity[1] = (-jacobian[1,1] + jacobian[2,4] - jacobian[3,3]) * ONE_OVER_4PI
    velocity[2] = (-jacobian[2,1] + jacobian[3,2] - jacobian[1,4]) * ONE_OVER_4PI
    velocity[3] = (-jacobian[3,1] + jacobian[1,3] - jacobian[2,2]) * ONE_OVER_4PI
    # jacobian of the velocity
    hessian = reshape(potential[i_POTENTIAL_HESSIAN],3,3,4)
    duidxj = fmm.@SMatrix [
        -hessian[1,1,1]+hessian[2,1,4]-hessian[3,1,3] -hessian[1,2,1]+hessian[2,2,4]-hessian[3,2,3] -hessian[1,3,1]+hessian[2,3,4]-hessian[3,3,3];
        -hessian[2,1,1]+hessian[3,1,2]-hessian[1,1,4] -hessian[2,2,1]+hessian[3,2,2]-hessian[1,2,4] -hessian[2,3,1]+hessian[3,3,2]-hessian[1,3,4];
        -hessian[3,1,1]+hessian[1,1,3]-hessian[2,1,2] -hessian[3,2,1]+hessian[1,2,3]-hessian[2,2,2] -hessian[3,3,1]+hessian[1,3,3]-hessian[2,3,2];
    ]
    # vorticity = @SVector [
    #     jacobian[2,3] - jacobian[3,2],
    #     jacobian[3,1] - jacobian[1,3],
    #     jacobian[1,2] - jacobian[2,1]
    # ]
    # stretching term (omega dot nabla)

    # total derivative of the strength due to vortex stretching
    fmm.mul!(view(velocity,i_STRETCHING_vortex), duidxj, body[i_STRENGTH_vortex])
    velocity[i_STRETCHING_vortex] .*= ONE_OVER_4PI

    return nothing
end

function update_velocity_stretching!(vortex_particles::VortexParticles)
    bodies = vortex_particles.bodies
    potential = vortex_particles.potential
    velocity_stretching = vortex_particles.velocity_stretching
    for i_body in 1:size(bodies)[2]
        this_velocity_stretching = view(velocity_stretching,:,i_body)
        this_body = view(bodies,:,i_body)
        this_potential = potential[:,i_body]
        update_velocity_stretching!(this_velocity_stretching, this_body, this_potential)
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
    for i_body in 1:size(bodies)[2]
        this_velocity_stretching = view(velocity_stretching,:,i_body)
        this_body = view(bodies,:,i_body)
        this_potential = potential[:,i_body]
        update_velocity_stretching!(this_velocity_stretching, this_body, this_potential)
        this_body[i_POSITION_vortex] .+= this_velocity_stretching[i_VELOCITY_vortex] * euler.dt
        this_body[i_STRENGTH_vortex] .+= this_velocity_stretching[i_STRETCHING_vortex] * euler.dt
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
