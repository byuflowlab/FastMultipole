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
function direct!(elements_tuple::Tuple)
    for source_elements in elements_tuple
        for target_elements in elements_tuple
            source_elements.direct!(target_elements.potential, target_elements.bodies[i_POSITION,:], source_elements.bodies)
        end
    end
end

function direct!(elements_tuple::Tuple, options::Options, targets_index, sources_index)
    for source_elements in elements_tuple[sources_index]
        for target_elements in elements_tuple[targets_index]
            source_elements.direct!(target_elements.potential, target_elements.bodies[i_POSITION,:], source_elements.bodies)
        end
    end
end
