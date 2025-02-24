# module VPM

import FastMultipole as fmm
import FastMultipole.WriteVTK
using FastMultipole
using SpecialFunctions:erf

#------- classic vortex particle method -------#

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
    strength = 2 .* rand(3, n_bodies) .- 1.0
    strength .-= 0.5
    strength .*= 2
    strength .*= strength_scale
    return VortexParticles(position, strength)
end

function reset!(system::VortexParticles{TF}) where TF
    system.potential .= zero(TF)
    system.velocity_stretching .= zero(TF)
end

abstract type IntegrationScheme end

struct Euler{TF} <: IntegrationScheme
    dt::TF
end

#------- overload compatibility functions -------#

function FastMultipole.source_system_to_buffer!(buffer, i_buffer, system::VortexParticles, i_body)
    buffer[1:3, i_buffer] .= system.bodies[i_body].position
    buffer[4, i_buffer] = system.bodies[i_body].sigma
    buffer[5:7, i_buffer] .= system.bodies[i_body].strength
end

function FastMultipole.data_per_body(system::VortexParticles)
    return 7
end

function FastMultipole.get_position(system::VortexParticles, i_body)
    return system.bodies[i_body].position
end

function FastMultipole.strength_dims(system::VortexParticles)
    return 3
end

function FastMultipole.get_n_bodies(system::VortexParticles)
    return length(system.bodies)
end

fmm.body_to_multipole!(system::VortexParticles, args...) = body_to_multipole!(Point{Vortex}, system, args...)

"""
Classical formulation so far.
"""
function fmm.direct!(target_system, target_index, derivatives_switch::FastMultipole.DerivativesSwitch{S,V,VG}, source_system::VortexParticles, source_buffer, source_index) where {S,V,VG}
    for j_source in source_index
        x_source = FastMultipole.get_position(source_buffer, j_source)
        Γx, Γy, Γz = FastMultipole.get_strength(source_buffer,  source_system, j_source)
        for i_target in target_index
            x_target = FastMultipole.get_position(target_system, i_target)
            dx, dy, dz = x_target - x_source
            r2 = dx * dx + dy * dy + dz * dz
            if FastMultipole.DEBUG[]
            println("\n\tdirect!: x_source = $x_source, x_target = $x_target")
        end
            if r2 > 0
                # distance away
                r = sqrt(r2)
                rinv = 1 / r
                rinv2 = rinv * rinv

                # useful denominator
                denom = FastMultipole.ONE_OVER_4π * rinv * rinv2

                # induced velocity
                vx, vy, vz = zero(r2), zero(r2), zero(r2)
                if V
                    vx = (dz * Γy - dy * Γz) * denom
                    vy = (dx * Γz - dz * Γx) * denom
                    vz = (dy * Γx - dx * Γy) * denom
                    FastMultipole.set_velocity!(target_system, i_target, SVector{3}(vx,vy,vz))
                end

                # velocity gradient
                if VG
                    denom *= rinv2
                    vxx = -3 * dx * (Γy * dz - Γz * dy) * denom
                    vxy = (-3 * dx * (Γz * dx - Γx * dz) + Γz * r2) * denom
                    vxz = (-3 * dx * (Γx * dy - Γy * dx) - Γy * r2) * denom
                    vyx = (-3 * dy * (Γy * dz - Γz * dy) - Γz * r2) * denom
                    vyy = -3 * dy * (Γz * dx - Γx * dz) * denom
                    vyz = (-3 * dy * (Γx * dy - Γy * dx) + Γx * r2) * denom
                    vzx = (-3 * dz * (Γy * dz - Γz * dy) + Γy * r2) * denom
                    vzy = (-3 * dz * (Γz * dx - Γx * dz) - Γx * r2) * denom
                    vzz = -3 * dz * (Γx * dy - Γy * dx) * denom
                    FastMultipole.set_velocity_gradient!(target_system, i_target, SMatrix{3,3}(vxx, vxy, vxz, vyx, vyy, vyz, vzx, vzy, vzz))
                end

            end
        end
    end
end

function FastMultipole.buffer_to_target_system!(target_system::VortexParticles, i_target, ::FastMultipole.DerivativesSwitch{PS,VS,GS}, target_buffer, i_buffer) where {PS,VS,GS}
    # retrieve fields
    TF = eltype(target_system)
    velocity = VS ? FastMultipole.get_velocity(target_buffer, i_buffer) : zero(SVector{3,TF})
    velocity_gradient = GS ? FastMultipole.get_velocity_gradient(target_buffer, i_buffer) : zero(SMatrix{3,3,TF,9})

    # update system
    if VS
        target_system.velocity_stretching[i_VELOCITY_vortex, i_target] .+= velocity
    end
    if GS
        target_system.potential[i_VELOCITY_GRADIENT_vortex, i_target] .+= reshape(velocity_gradient, 9)
    end
end

FastMultipole.get_n_bodies(vp::VortexParticles) = length(vp.bodies)

Base.eltype(::VortexParticles{TF}) where TF = TF

#------- additional functions -------#

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
