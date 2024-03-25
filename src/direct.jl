"""
    direct!(systems; derivatives_switches)

Applies all interactions of `systems` acting on itself without multipole acceleration.

# Arguments

- `systems`: either

    - a system object for which compatibility functions have been overloaded, or
    - a tuple of system objects for which compatibility functions have been overloaded

# Optional Arguments

- `derivatives_switches::DerivativesSwitch`: determines whether to include the scalar potential, vector potential, velocity, and/or velocity gradient when solving the n-body problem; either

    - (if `systems` is not a `::Tuple` of systems) an instance of [`DerivativesSwitch`](@ref)
    - (if `systems` is a `::Tuple` of systems) a `::Tuple` of [`DerivativesSwitch`](@ref) of length `length(systems)`

"""
function direct!(systems::Tuple; derivatives_switches=DerivativesSwitch(true,true,true,true,systems))
    for source_system in systems
        for (target_system, derivatives_switch) in zip(systems, derivatives_switches)
            _direct!(target_system, 1:get_n_bodies(target_system), derivatives_switch, source_system, 1:get_n_bodies(source_system))
        end
    end
end

function direct!(system; derivatives_switch=DerivativesSwitch{true,true,true,true}())
    direct!(system, system; derivatives_switch)
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

- `derivatives_switches::DerivativesSwitch`: determines whether to include the scalar potential, vector potential, velocity, and/or velocity gradient for the `target_system`; either

    - (if `target_system` is not a `::Tuple` of systems) an instance of [`DerivativesSwitch`](@ref)
    - (if `target_system` is a `::Tuple` of systems) a `::Tuple` of [`DerivativesSwitch`](@ref) of length `length(target_systems)`

"""
@inline function direct!(target_system, source_system; derivatives_switch=DerivativesSwitch{true,true,true,true}())
    _direct!(target_system, 1:get_n_bodies(target_system), derivatives_switch, source_system, 1:get_n_bodies(source_system))
end

function direct!(target_systems::Tuple, source_systems::Tuple; derivatives_switches=Tuple(DerivativesSwitch{true,true,true,true}() for _ in target_systems))
    for source_system in source_systems
        for (target_system, derivatives_switch) in zip(target_systems, derivatives_switches)
            direct!(target_system, source_system; derivatives_switch)
        end
    end
end

function direct!(target_systems::Tuple, source_system; derivatives_switches=Tuple(DerivativesSwitch{true,true,true,true}() for _ in target_systems))
    for (target_system, derivatives_switch) in zip(target_systems, derivatives_switches)
        direct!(target_system, source_system; derivatives_switch)
    end
end

function direct!(target_systems, source_system::Tuple; derivatives_switch=DerivativesSwitch{true,true,true,true}())
    for source_system in source_systems
        direct!(target_system, source_system; derivatives_switch)
    end
end

#####
##### private methods for dispatch
#####
@inline function _direct!(target_system, target_bodies_index, derivatives_switch, source_system, source_bodies_index)
    direct!(target_system, target_bodies_index, derivatives_switch, source_system, source_bodies_index)
end

@inline function _direct!(target_system::SortWrapper, target_bodies_index, derivatives_switch, source_system, source_bodies_index)
    direct!(target_system.system, target_system.index[target_bodies_index], derivatives_switch, source_system, source_bodies_index)
end

@inline function _direct!(target_system, target_bodies_index, derivatives_switch, source_system::SortWrapper, source_bodies_index)
    direct!(target_system, target_bodies_index, source_system.system, derivatives_switch, source_system.index[source_bodies_index])
end

@inline function _direct!(target_system::SortWrapper, target_bodies_index, derivatives_switch, source_system::SortWrapper, source_bodies_index)
    direct!(target_system.system, target_system.index[target_bodies_index], derivatives_switch, source_system.system, source_system.index[source_bodies_index])
end

