# """
#     direct!(targets, source)

# Uses a naive direct algorithm to evaluate the influence of all sources on all targets.

# # Inputs

# - `targets::Vector{element}`- a vector of target elements
# - `sources::Vector{element}`- a vector of source elements

# """
# function direct!(targets::Vector, sources::Vector, P2P!::Function)
#     for target in targets
#         for source in sources
#             P2P!(target, source)
#         end
#     end
# end

"""
    direct!(elements)

Direct calculation of induced potential (no FMM acceleration).
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

