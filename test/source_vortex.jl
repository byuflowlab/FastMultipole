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
const i_STRENGTH_scalr = 7
const i_POTENTIAL_SCALAR = 1:1
const i_POTENTIAL_VECTOR = 2:4
const i_VELOCITY_GRADIENT_vortex = 5:13
const i_VELOCITY_vortex = 1:3
const i_STRETCHING_vortex = 4:6

struct SourceVorton{TF}
    position::SVector{3,TF}
    strength::SVector{4,TF}
    sigma::TF
end

struct SourceVortons{TF}
    bodies::Vector{SourceVorton{TF}}
    potential::Matrix{TF}
    velocity_stretching::Matrix{TF}
end

function generate_source_vortex(seed, n_bodies; strength_scale=1/n_bodies)
    Random.seed!(seed)
    position = rand(3, n_bodies)
    strength = 2 .* rand(4, n_bodies) .- 1.0
    strength .-= 0.5
    strength .*= 2
    strength .*= strength_scale
    return SourceVortons(position, strength)
end

abstract type IntegrationScheme end

struct Euler{TF} <: IntegrationScheme
    dt::TF
end

#####
##### overload compatibility functions
#####

Base.getindex(vp::SourceVortons, i, ::Position) = vp.bodies[i].position
Base.getindex(vp::SourceVortons, i, ::Radius) = vp.bodies[i].sigma
#Base.getindex(vp::SourceVortons, i, ::fmm.VectorPotential) = view(vp.potential,2:4,i)
Base.getindex(vp::SourceVortons, i, ::ScalarPotential) = vp.potential[1,i]
Base.getindex(vp::SourceVortons, i, ::Velocity) = view(vp.velocity_stretching,i_VELOCITY_vortex,i)
Base.getindex(vp::SourceVortons, i, ::VelocityGradient) = reshape(view(vp.potential,i_VELOCITY_GRADIENT_vortex,i),3,3)
Base.getindex(vp::SourceVortons, i, ::Strength) = vp.bodies[i].strength
Base.getindex(vp::SourceVortons, i, ::FastMultipole.Body) = vp.bodies[i], vp.potential[:,i], vp.velocity_stretching[:,i]

function Base.setindex!(vp::SourceVortons, val, i, ::Position)
    vp.bodies[i] = Vorton(val, vp.bodies[i].strength, vp.bodies[i].sigma)
    return nothing
end
function Base.setindex!(vp::SourceVortons, val, i, ::Strength)
    vp.bodies[i] = Vorton(vp.bodies[i].position, val, vp.bodies[i].sigma)
    return nothing
end
function Base.setindex!(vp::SourceVortons, val, i, ::FastMultipole.Body)
    body, potential, velocity = val
    vp.bodies[i] = body
    vp.potential[:,i] .= potential
    vp.velocity_stretching[:,i] .= velocity
    return nothing
end
function Base.setindex!(vp::SourceVortons, val, i, ::ScalarPotential)
    # vp.potential[i_POTENTIAL[1],i] = val
    return nothing
end
#function Base.setindex!(vp::SourceVortons, val, i, ::fmm.VectorPotential)
#    vp.potential[i_POTENTIAL_VECTOR,i] .= val
#end
function Base.setindex!(vp::SourceVortons, val, i, ::Velocity)
    vp.velocity_stretching[i_VELOCITY_vortex,i] .= val
end
function Base.setindex!(vp::SourceVortons, val, i, ::VelocityGradient)
    vp.potential[i_VELOCITY_GRADIENT_vortex,i] .= reshape(val,9)
end
FastMultipole.get_n_bodies(vp::SourceVortons) = length(vp.bodies)
Base.eltype(::SourceVortons{TF}) where TF = TF

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
function fmm.direct!(target_system, target_index, derivatives_switch::FastMultipole.DerivativesSwitch{S,V,VG}, source_system::SourceVortons, source_index) where {S,V,VG}
    for j_source in source_index
        x_source = source_system[j_source,Position()]
        σ, Γx, Γy, Γz = source_system[j_source,Strength()]
        for i_target in target_index
            x_target = target_system[i_target,Position()]
            dx, dy, dz = x_target - x_source
            r2 = dx * dx + dy * dy + dz * dz
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
                    # vector influence
                    vx = (dz * Γy - dy * Γz) * denom
                    vy = (dx * Γz - dz * Γx) * denom
                    vz = (dy * Γx - dx * Γy) * denom

                    # scalar influence
                    vx += dx * σ * denom
                    vy += dy * σ * denom
                    vz += dz * σ * denom
                end

                # # velocity gradient
                # vxx, vxy, vxz = zero(r), zero(r), zero(r)
                # vyx, vyy, vyz = zero(r), zero(r), zero(r)
                # vzx, vzy, vzz = zero(r), zero(r), zero(r)
                # if VG
                #     denom *= rinv2
                #     vxx = -3 * dx * (Γy * dz - Γz * dy) * denom
                #     vxy = (-3 * dx * (Γz * dx - Γx * dz) + Γz * r2) * denom
                #     vxz = (-3 * dx * (Γx * dy - Γy * dx) - Γy * r2) * denom
                #     vyx = (-3 * dy * (Γy * dz - Γz * dy) - Γz * r2) * denom
                #     vyy = -3 * dy * (Γz * dx - Γx * dz) * denom
                #     vyz = (-3 * dy * (Γx * dy - Γy * dx) + Γx * r2) * denom
                #     vzx = (-3 * dz * (Γy * dz - Γz * dy) + Γy * r2) * denom
                #     vzy = (-3 * dz * (Γz * dx - Γx * dz) - Γx * r2) * denom
                #     vzz = -3 * dz * (Γx * dy - Γy * dx) * denom
                #     @warn "velocity gradient only partially implemented for SourceVortons"
                # end

                v = target_system[i_target, Velocity()]
                target_system[i_target, Velocity()] = v + SVector{3}(vx,vy,vz)

                # vgxx, vgxy, vgxz, vgyx, vgyy, vgyz, vgzx, vgzy, vgzz = target_system[i_target, VelocityGradient()]
                # target_system[i_target, VelocityGradient()] = SMatrix{3,3}(vxx+vgxx, vxy+vgxy, vxz+vgxz, vyx+vgyx, vyy+vgyy, vyz+vgyz, vzx+vgzx, vzy+vgzy, vzz+vgzz)
            end
        end
    end
end

fmm.buffer_element(system::SourceVortons) = (deepcopy(system.bodies[1]),zeros(eltype(system),52),zeros(eltype(system),6))

function fmm.body_to_multipole!(system::SourceVortons, args...)
    body_to_multipole!(Point{SourceVortex}, system, args...)
end

function SourceVortons(position, strength;
    N = size(position)[2],
    potential = zeros(i_VELOCITY_GRADIENT_vortex[end],N),
    velocity_stretching = zeros(3+3,N),
)
    @assert size(position)[1] == 3
    @assert size(strength)[1] == 4
    bodies = [SourceVorton(SVector{3}(position[:,i]), SVector{4}(strength[:,i]), 0.0) for i in 1:size(position)[2]]
    return SourceVortons(bodies, potential, velocity_stretching)
end

function SourceVortons(bodies;
    N = size(bodies)[2],
    potential = zeros(i_VELOCITY_GRADIENT_vortex[end],N),
    velocity_stretching = zeros(3+3,N)
)
    bodies = [SourceVorton(SVector{3}(bodies[1:3,i]), SVector{3}(bodies[5:7,i]), bodies[8,i],  bodies[4,i]) for i in 1:size(bodies)[2]]
    return SourceVortons(bodies, potential, velocity_stretching)
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

function update_velocity_stretching!(vortex_particles::SourceVortons)
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

function update_velocity_stretching_new!(vortex_particles::SourceVortons)
    for i_body in 1:fmm.get_n_bodies(vortex_particles)
        update_velocity_stretching_new!(vortex_particles, i_body)
    end
end
#=
"Assumes velocities and potential and derivatives are already calculated."
function (euler::Euler)(vortex_particles::SourceVortons, fmm_options, direct)
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

function convect!(vortex_particles::SourceVortons, nsteps;
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

function save_vtk(filename, vortex_particles::SourceVortons, nt=0; compress=false)

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
=#
