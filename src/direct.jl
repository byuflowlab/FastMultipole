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
function direct!(elements)
    elements.direct!(view(elements.potential,:,:), elements.bodies[i_POSITION,:], elements.bodies)
end
