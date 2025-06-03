function ProbeSystem(n_bodies, TF=Float64)
    position = zeros(SVector{3,TF}, n_bodies)
    scalar_potential = zeros(TF, n_bodies)
    velocity = zeros(SVector{3,TF}, n_bodies)
    velocity_gradient = zeros(SMatrix{3,3,TF,9}, n_bodies)
    return ProbeSystem{TF}(position, scalar_potential, velocity, velocity_gradient)
end

function reset!(system::ProbeSystem{TF}) where TF
    system.scalar_potential .= zero(TF)
    for i in eachindex(system.velocity)
        system.velocity[i] = zero(SVector{3,TF})
    end
    for i in eachindex(system.velocity_gradient)
        system.velocity_gradient[i] = zero(SMatrix{3,3,TF,9})
    end
end

#------- FastMultipole compatibility functions -------#

Base.eltype(::ProbeSystem{TF}) where TF = TF

function FastMultipole.source_system_to_buffer!(buffer, i_buffer, system::ProbeSystem, i_body)
    buffer[1:3, i_buffer] .= system.bodies[i_body].position
    buffer[4, i_buffer] = zero(TF)
    buffer[5, i_buffer] = zero(TF)
end

function FastMultipole.data_per_body(system::ProbeSystem)
    return 5
end

function FastMultipole.get_position(system::ProbeSystem, i)
    return system.position[i]
end

function FastMultipole.strength_dims(system::ProbeSystem)
    return 1
end

FastMultipole.get_n_bodies(system::ProbeSystem) = length(system.position)

function FastMultipole.body_to_multipole!(system::ProbeSystem, args...)
    return nothing
end

function FastMultipole.direct!(target_system, target_index, derivatives_switch, source_system::ProbeSystem, source_buffer, source_index)
    return nothing
end

function FastMultipole.buffer_to_target_system!(target_system::ProbeSystem, i_target, ::FastMultipole.DerivativesSwitch{PS,VS,GS}, target_buffer, i_buffer) where {PS,VS,GS}
    TF = eltype(target_buffer)
    scalar_potential = PS ? FastMultipole.get_scalar_potential(target_buffer, i_buffer) : zero(TF)
    velocity = VS ? FastMultipole.get_velocity(target_buffer, i_buffer) : zero(SVector{3,TF})
    velocity_gradient = GS ? FastMultipole.get_velocity_gradient(target_buffer, i_buffer) : zero(SMatrix{3,3,TF,9})

    target_system.scalar_potential[i_target] = scalar_potential
    target_system.velocity[i_target] = velocity
    target_system.velocity_gradient[i_target] = velocity_gradient
end
