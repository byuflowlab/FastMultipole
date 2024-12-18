# module VPM

import FastMultipole as fmm
import FastMultipole.WriteVTK
using FastMultipole
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
const i_VELOCITY_GRADIENT_vortex = 5:13
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

function generate_vortex(seed, n_bodies; strength_scale=1/n_bodies)
    Random.seed!(seed)
    position = rand(3, n_bodies)
    strength = rand(3, n_bodies) .* strength_scale
    return VortexParticles(position, strength)
end

abstract type IntegrationScheme end

struct Euler{TF} <: IntegrationScheme
    dt::TF
end

#####
##### overload compatibility functions
#####

Base.getindex(vp::VortexParticles, i, ::Position) = vp.bodies[i].position
Base.getindex(vp::VortexParticles, i, ::Radius) = vp.bodies[i].sigma
#Base.getindex(vp::VortexParticles, i, ::fmm.VectorPotential) = view(vp.potential,2:4,i)
Base.getindex(vp::VortexParticles, i, ::ScalarPotential) = vp.potential[1,i]
Base.getindex(vp::VortexParticles, i, ::Velocity) = view(vp.velocity_stretching,i_VELOCITY_vortex,i)
Base.getindex(vp::VortexParticles, i, ::VelocityGradient) = reshape(view(vp.potential,i_VELOCITY_GRADIENT_vortex,i),3,3)
Base.getindex(vp::VortexParticles, i, ::Strength) = vp.bodies[i].strength
Base.getindex(vp::VortexParticles, i, ::FastMultipole.Body) = vp.bodies[i], vp.potential[:,i], vp.velocity_stretching[:,i]

function Base.setindex!(vp::VortexParticles, val, i, ::Position)
    vp.bodies[i] = Vorton(val, vp.bodies[i].strength, vp.bodies[i].sigma)
    return nothing
end
function Base.setindex!(vp::VortexParticles, val, i, ::Strength)
    vp.bodies[i] = Vorton(vp.bodies[i].position, val, vp.bodies[i].sigma)
    return nothing
end
function Base.setindex!(vp::VortexParticles, val, i, ::FastMultipole.Body)
    body, potential, velocity = val
    vp.bodies[i] = body
    vp.potential[:,i] .= potential
    vp.velocity_stretching[:,i] .= velocity
    return nothing
end
function Base.setindex!(vp::VortexParticles, val, i, ::ScalarPotential)
    # vp.potential[i_POTENTIAL[1],i] = val
    return nothing
end
#function Base.setindex!(vp::VortexParticles, val, i, ::fmm.VectorPotential)
#    vp.potential[i_POTENTIAL_VECTOR,i] .= val
#end
function Base.setindex!(vp::VortexParticles, val, i, ::Velocity)
    vp.velocity_stretching[i_VELOCITY_vortex,i] .= val
end
function Base.setindex!(vp::VortexParticles, val, i, ::VelocityGradient)
    vp.potential[i_VELOCITY_GRADIENT_vortex,i] .= reshape(val,9)
end
FastMultipole.get_n_bodies(vp::VortexParticles) = length(vp.bodies)
Base.eltype(::VortexParticles{TF}) where TF = TF

function flatten_derivatives!(jacobian, hessian, derivatives_switch::DerivativesSwitch{PS,VS,GS}) where {PS,VS,GS}
    if VS
        # velocity
        jacobian[1,1] = -jacobian[1,1] + jacobian[2,4] - jacobian[3,3]
        jacobian[2,1] = -jacobian[2,1] + jacobian[3,2] - jacobian[1,4]
        jacobian[3,1] = -jacobian[3,1] + jacobian[1,3] - jacobian[2,2]
    end

    if GS
        # velocity gradient
        hessian[1,1,1] = -hessian[1,1,1] + hessian[2,1,4] - hessian[3,1,3]
        hessian[2,1,1] = -hessian[2,1,1] + hessian[3,1,2] - hessian[1,1,4]
        hessian[3,1,1] = -hessian[3,1,1] + hessian[1,1,3] - hessian[2,1,2]
        hessian[1,2,1] = -hessian[1,2,1] + hessian[2,2,4] - hessian[3,2,3]
        hessian[2,2,1] = -hessian[2,2,1] + hessian[3,2,2] - hessian[1,2,4]
        hessian[3,2,1] = -hessian[3,2,1] + hessian[1,2,3] - hessian[2,2,2]
        hessian[1,3,1] = -hessian[1,3,1] + hessian[2,3,4] - hessian[3,3,3]
        hessian[2,3,1] = -hessian[2,3,1] + hessian[3,3,2] - hessian[1,3,4]
        hessian[3,3,1] = -hessian[3,3,1] + hessian[1,3,3] - hessian[2,3,2]
    end
end

"""
Classical formulation so far.
"""
function fmm.direct!(target_system, target_index, derivatives_switch, source_system::VortexParticles, source_index)
    jacobian = zeros(eltype(target_system),3,4)
    hessian = zeros(eltype(target_system),3,3,4)
    for j_source in source_index
        x_source = source_system[j_source,Position()]
        for i_target in target_index
            jacobian .*= 0.0
            hessian .*= 0.0
            x_target = target_system[i_target,Position()]
            dx = x_target - x_source
            r = sqrt(dx' * dx)
            if r > 0
                # calculate induced potential
                gamma_over_R = source_system.bodies[j_source].strength / r * ONE_OVER_4PI
                #target_system[i_target,fmm.VECTOR_POTENTIAL] .+= gamma_over_R

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
            flatten_derivatives!(jacobian, hessian, derivatives_switch)

            target_system[i_target,fmm.VELOCITY] .+= view(jacobian,:,1)
            target_system[i_target,fmm.VELOCITY_GRADIENT] .+= view(hessian,:,:,1)
        end

    end
end

fmm.buffer_element(system::VortexParticles) = (deepcopy(system.bodies[1]),zeros(eltype(system),52),zeros(eltype(system),6))

fmm.body_to_multipole!(system::VortexParticles, args...) = body_to_multipole!(Point{Vortex}, system, args...)

function VortexParticles(position, strength;
    N = size(position)[2],
    potential = zeros(i_VELOCITY_GRADIENT_vortex[end],N),
    velocity_stretching = zeros(3+3,N),
)
    @assert size(position)[1] == 3
    @assert size(strength)[1] == 3
    bodies = [Vorton(SVector{3}(position[:,i]), SVector{3}(strength[:,i]), 0.0) for i in 1:size(position)[2]]
    return VortexParticles(bodies, potential, velocity_stretching)
end

function VortexParticles(bodies;
    N = size(bodies)[2],
    potential = zeros(i_VELOCITY_GRADIENT_vortex[end],N),
    velocity_stretching = zeros(3+3,N)
)
    bodies = [Vorton(SVector{3}(bodies[1:3,i]), SVector{3}(bodies[5:7,i]), bodies[4,i]) for i in 1:size(bodies)[2]]
    return VortexParticles(bodies, potential, velocity_stretching)
end

@inline function update_velocity_stretching!(system, i_body)
    # vorticity = @SVector [
    #     jacobian[2,3] - jacobian[3,2],
    #     jacobian[3,1] - jacobian[1,3],
    #     jacobian[1,2] - jacobian[2,1]
    # ]
    # stretching term (omega dot nabla)

    # total derivative of the strength due to vortex stretching
    duidxj = reshape(view(system.potential,i_VELOCITY_GRADIENT_vortex,i_body),3,3)
    # @show duidxj system.bodies[i_body].strength
    fmm.mul!(view(system.velocity_stretching,i_STRETCHING_vortex,i_body), duidxj, system.bodies[i_body].strength)

    return nothing
end

function update_velocity_stretching!(vortex_particles::VortexParticles)
    for i_body in 1:fmm.get_n_bodies(vortex_particles)
        update_velocity_stretching!(vortex_particles, i_body)
    end
end

@inline function update_velocity_stretching_new!(system, i_body)
    # vorticity = @SVector [
    #     jacobian[2,3] - jacobian[3,2],
    #     jacobian[3,1] - jacobian[1,3],
    #     jacobian[1,2] - jacobian[2,1]
    # ]
    # stretching term (omega dot nabla)

    velocity_gradient = reshape(view(system.potential,i_VELOCITY_GRADIENT_vortex,i_body),3,3)

    system.velocity_stretching[4] = dot(velocity_gradient[1,1:3], system.bodies[i_body].strength)
    system.velocity_stretching[5] = dot(velocity_gradient[2,1:3], system.bodies[i_body].strength)
    system.velocity_stretching[6] = dot(velocity_gradient[3,1:3], system.bodies[i_body].strength)

    return nothing
end

function update_velocity_stretching_new!(vortex_particles::VortexParticles)
    for i_body in 1:fmm.get_n_bodies(vortex_particles)
        update_velocity_stretching_new!(vortex_particles, i_body)
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
        fmm.fmm!((vortex_particles,); fmm_options...)
    end

    # convect bodies
    bodies = vortex_particles.bodies
    potential = vortex_particles.potential
    velocity_stretching = vortex_particles.velocity_stretching
    update_velocity_stretching!(vortex_particles)
    for i_body in 1:fmm.get_n_bodies(vortex_particles)
        vortex_particles[i_body, Position()] = vortex_particles[i_body, Position()] + vortex_particles.velocity_stretching[i_VELOCITY_vortex,i_body] * euler.dt
        vortex_particles[i_body, Strength()] = vortex_particles[i_body, Strength()] + vortex_particles.velocity_stretching[i_STRETCHING_vortex] * euler.dt
    end
end

function convect!(vortex_particles::VortexParticles, nsteps;
        # integration options
        integrate::IntegrationScheme=Euler(1.0),
        # fmm options
        fmm_p=4, fmm_ncrit=50, fmm_multipole_threshold=0.5,
        direct::Bool=false,
        # save options
        save::Bool=true, filename::String="default", compress::Bool=false,
    )
    fmm_options = (; expansion_order=fmm_p, leaf_size=fmm_ncrit, multipole_threshold=fmm_multipole_threshold, lamb_helmholtz=true)
    save && save_vtk(filename, vortex_particles; compress)
    for istep in 1:nsteps
        integrate(vortex_particles, fmm_options, direct)
        save && save_vtk(filename, vortex_particles, istep; compress)
    end
    return nothing
end

function save_vtk(filename, vortex_particles::VortexParticles, nt=0; compress=false)

    n_bodies = length(vortex_particles.bodies)

    positions = reshape([vortex_particles[j, Position()][i] for i in 1:3, j in 1:n_bodies], 3, n_bodies, 1, 1)
    vectorstrength = reshape([vortex_particles[j, Strength()][i] for i in 1:3, j in 1:n_bodies], 3, n_bodies, 1, 1)
    scalarpotential = reshape([vortex_particles[j, ScalarPotential()] for j in 1:n_bodies], n_bodies, 1, 1)
    vectorpotential = reshape(vortex_particles.potential[i_POTENTIAL_VECTOR,:],3,n_bodies,1,1)
    velocity = reshape(vortex_particles.velocity_stretching[i_VELOCITY_vortex,:],3,n_bodies,1,1)
    stretching = reshape(vortex_particles.velocity_stretching[i_STRETCHING_vortex,:],3,n_bodies,1,1)

    WriteVTK.vtk_grid(filename*"."*string(nt)*".vts", positions; compress) do vtk
        vtk["vector strength"] = vectorstrength
        vtk["scalar potential"] = scalarpotential
        vtk["vector potential"] = vectorpotential
        vtk["velocity"] = velocity
        vtk["stretching"] = stretching
    end

end
