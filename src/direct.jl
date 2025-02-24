"""
    direct!(systems; derivatives_switches)

Applies all interactions of `systems` acting on itself without multipole acceleration.

# Arguments

- `systems`: either

    - a system object for which compatibility functions have been overloaded, or
    - a tuple of system objects for which compatibility functions have been overloaded

# Optional Arguments

- `scalar_potential::Bool`: either a `::Bool` or a `::AbstractVector{Bool}` of length `length(target_systems)` indicating whether each system should receive a scalar potential from `source_systems`
- `velocity::Bool`: either a `::Bool` or a `::AbstractVector{Bool}` of length `length(target_systems)` indicating whether each system should receive a velocity from `source_systems`
- `velocity_gradient::Bool`: either a `::Bool` or a `::AbstractVector{Bool}` of length `length(target_systems)` indicating whether each system should receive a velocity gradient from `source_systems`

"""
function direct!(systems::Tuple; args...)
    direct!(systems, systems, args...)
end

function direct!(system; args...)
    direct!((system,); args...)
end

function direct!(target_system, source_system; args...)
    target_system = to_tuple(target_system)
    source_system = to_tuple(source_system)
    direct!((target_system,), (source_system,); args...)
end

function direct!(target_systems::Tuple, source_systems::Tuple; target_buffers=nothing, source_buffers=nothing, scalar_potential=fill(true, length(target_systems)), velocity=fill(true, length(target_systems)), velocity_gradient=fill(true, length(target_systems)))
    # set up target buffers
    if isnothing(target_buffers)
        target_buffers = allocate_buffers(target_systems, true)
        target_to_buffer!(target_buffers, target_systems)
    end

    # set up source buffers
    if isnothing(source_buffers)
        source_buffers = allocate_buffers(source_systems, false)
        system_to_buffer!(source_buffers, source_systems)
    end

    # ensure derivative switch information is a vector
    scalar_potential = to_vector(scalar_potential, length(target_systems))
    velocity = to_vector(velocity, length(target_systems))
    velocity_gradient = to_vector(velocity_gradient, length(target_systems))
    derivatives_switches = DerivativesSwitch(scalar_potential, velocity, velocity_gradient)

    for (source_system, source_buffer) in zip(source_systems, source_buffers)
        for (target_system, target_buffer, derivatives_switch) in zip(target_systems, target_buffers, derivatives_switches)
            direct!(target_buffer, 1:get_n_bodies(target_system), derivatives_switch, source_system, source_buffer, 1:get_n_bodies(source_system))
        end
    end

    # update target systems
    if DEBUG[]
        println("\nASLAN")
        @show target_systems[1].velocity_stretching[1:3,:]
    end
    buffer_to_target!(target_systems, target_buffers, derivatives_switches)
    if DEBUG[]
        println("\nthere:")
        @show target_systems[1].velocity_stretching[1:3,:]
    end

end

