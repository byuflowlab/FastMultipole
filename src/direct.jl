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
            direct!(target_system, 1:length(target_system), source_system, 1:length(source_system))
        end
    end
end

function direct!(systems::Tuple, targets_index, sources_index)
    for source_system in systems[sources_index]
        for target_system in systems[targets_index]
            direct!(target_system, 1:length(target_system), source_system, 1:length(source_system))
        end
    end
end

function direct!(system)
    direct!(system, 1:length(system), system, 1:length(system))
end

function direct!(target_system, target_indices, source_system, source_indices)
    @warn "direct! function not implemented for source type $(typeof(source_system)) and will do nothing; overload FLOWFMM.direct!"
    return nothing
end
