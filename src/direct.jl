const DEBUG_TOGGLE = Array{Bool,0}(undef)
DEBUG_TOGGLE[] = false

"""
    direct!(systems; derivatives_switches)

Applies all interactions of `systems` acting on itself without multipole acceleration.

# Arguments

- `systems`: either

    - a system object for which compatibility functions have been overloaded, or
    - a tuple of system objects for which compatibility functions have been overloaded

# Optional Arguments

- `scalar_potential::Bool`: either a `::Bool` or a `::AbstractVector{Bool}` of length `length(target_systems)` indicating whether each system should receive a scalar potential from `source_systems`
- `vector_potential::Bool`: either a `::Bool` or a `::AbstractVector{Bool}` of length `length(target_systems)` indicating whether each system should receive a vector potential from `source_systems`
- `velocity::Bool`: either a `::Bool` or a `::AbstractVector{Bool}` of length `length(target_systems)` indicating whether each system should receive a velocity from `source_systems`
- `velocity_gradient::Bool`: either a `::Bool` or a `::AbstractVector{Bool}` of length `length(target_systems)` indicating whether each system should receive a velocity gradient from `source_systems`

"""
function direct!(systems::Tuple; scalar_potential=fill(true, length(systems)), vector_potential=fill(true, length(systems)), velocity=fill(true, length(systems)), velocity_gradient=fill(true, length(systems)))
    derivatives_switches = DerivativesSwitch(scalar_potential, vector_potential, velocity, velocity_gradient)
    for source_system in systems
        for (target_system, derivatives_switch) in zip(systems, derivatives_switches)
            direct!(target_system, (1:get_n_bodies(target_system)), derivatives_switch, source_system, 1:get_n_bodies(source_system))
        end
    end
end

function direct!(system; scalar_potential=true, vector_potential=true, velocity=true, velocity_gradient=true)
    direct!(system, system; scalar_potential, vector_potential, velocity, velocity_gradient)
end

"""
    direct!(target_system, source_system; derivatives_switches)

Applies all interactions of `source_system` acting on `target_system` without multipole acceleration.

# Arguments

- `target_system`: either

    - a system object for which compatibility functions have been overloaded, or
    - a tuple of system objects for which compatibility functions have been overloaded

- `source_system`: either

    - a system object for which compatibility functions have been overloaded, or
    - a tuple of system objects for which compatibility functions have been overloaded

# Optional Arguments

- `scalar_potential::Bool`: either a `::Bool` or a `::AbstractVector{Bool}` of length `length(target_systems)` indicating whether each system should receive a scalar potential from `source_systems`
- `vector_potential::Bool`: either a `::Bool` or a `::AbstractVector{Bool}` of length `length(target_systems)` indicating whether each system should receive a vector potential from `source_systems`
- `velocity::Bool`: either a `::Bool` or a `::AbstractVector{Bool}` of length `length(target_systems)` indicating whether each system should receive a velocity from `source_systems`
- `velocity_gradient::Bool`: either a `::Bool` or a `::AbstractVector{Bool}` of length `length(target_systems)` indicating whether each system should receive a velocity gradient from `source_systems`

"""
@inline function direct!(target_system, source_system; scalar_potential=true, vector_potential=true, velocity=true, velocity_gradient=true)
    derivatives_switch = DerivativesSwitch(scalar_potential, vector_potential, velocity, velocity_gradient)
    _direct!(target_system, (1:get_n_bodies(target_system)), derivatives_switch, source_system, 1:get_n_bodies(source_system))
end

function direct!(target_systems::Tuple, source_systems::Tuple; scalar_potential=fill(true, length(target_systems)), vector_potential=fill(true, length(target_systems)), velocity=fill(true, length(target_systems)), velocity_gradient=fill(true, length(target_systems)))
    for source_system in source_systems
        for (target_system, sp, vp, v, vg) in zip(target_systems, scalar_potential, vector_potential, velocity, velocity_gradient)
            direct!(target_system, source_system; scalar_potential=sp, vector_potential=vp, velocity=v, velocity_gradient=vg)
        end
    end
end

function direct!(target_systems::Tuple, source_system; scalar_potential=fill(true,length(target_systems)), vector_potential=fill(true,length(target_systems)), velocity=fill(true,length(target_systems)), velocity_gradient=fill(true,length(target_systems)))
    for (target_system, sp, vp, v, vg) in zip(target_systems, scalar_potential, vector_potential, velocity, velocity_gradient)
        direct!(target_system, source_system; scalar_potential=sp, vector_potential=vp, velocity=v, velocity_gradient=vg)
    end
end

function direct!(target_system, source_systems::Tuple; scalar_potential=true, vector_potential=true, velocity=true, velocity_gradient=true)
    for source_system in source_systems
        direct!(target_system, source_system; scalar_potential, vector_potential, velocity, velocity_gradient)
    end
end

#####
##### private methods for dispatch
#####
function _direct!(target_system, target_indices, derivatives_switch, source_system, source_index)
    direct!(target_system, target_indices, derivatives_switch, source_system, source_index)
end

function _direct!(target_system, target_indices, derivatives_switch, source_system::SortWrapper, source_index)
    for i in source_index
        direct!(target_system, target_indices, derivatives_switch, source_system.system, source_system.index[i])
    end
end

function _direct!(target_system::SortWrapper, target_indices, derivatives_switch, source_system, source_index)
    for target_index in target_indices
        for i in target_index
            direct!(target_system.system, (i:i), derivatives_switch, source_system, source_index)
        end
    end
end

function _direct!(target_system::SortWrapper, target_indices, derivatives_switch, source_system::SortWrapper, source_index)
    for i in source_index
        direct!(target_system, target_indices, derivatives_switch, source_system.system, source_system.index[i])
    end
end
