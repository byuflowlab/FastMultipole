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
function direct!(systems::Tuple)
    for source_system in systems
        for target_system in systems
            _direct!(target_system, 1:get_n_bodies(target_system), source_system, 1:get_n_bodies(source_system))
        end
    end
end

function direct!(system)
    direct!(system, system)
end

@inline function direct!(target_system, source_system)
    _direct!(target_system, 1:get_n_bodies(target_system), source_system, 1:get_n_bodies(source_system))
end

function direct!(target_systems::Tuple, source_systems::Tuple)
    for source_system in source_systems
        for target_system in target_systems
            direct!(target_system, source_system)
        end
    end
end

function direct!(target_systems::Tuple, source_system)
    for target_system in target_systems
        direct!(target_system, source_system)
    end
end

function direct!(target_systems, source_system::Tuple)
    for source_system in source_systems
        direct!(target_system, source_system)
    end
end

#####
##### private methods for dispatch
#####
@inline function _direct!(target_system, target_bodies_index, source_system, source_bodies_index)
    direct!(target_system, target_bodies_index, source_system, source_bodies_index)
end

@inline function _direct!(target_system::SortWrapper, target_bodies_index, source_system, source_bodies_index)
    direct!(target_system.system, target_system.index[target_bodies_index], source_system, source_bodies_index)
end

@inline function _direct!(target_system, target_bodies_index, source_system::SortWrapper, source_bodies_index)
    direct!(target_system, target_bodies_index, source_system.system, source_system.index[source_bodies_index])
end

@inline function _direct!(target_system::SortWrapper, target_bodies_index, source_system::SortWrapper, source_bodies_index)
    direct!(target_system.system, target_system.index[target_bodies_index], source_system.system, source_system.index[source_bodies_index])
end

