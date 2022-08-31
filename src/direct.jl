"""
    direct!(targets, source)

Uses a naive direct algorithm to evaluate the influence of all sources on all targets.

# Inputs

- `targets::Vector{element}`- a vector of target elements
- `sources::Vector{element}`- a vector of source elements

"""
function direct!(targets::Vector, sources::Vector)
    for (i,target) in enumerate(targets)
        for (j,source) in enumerate(sources)
            kernel!(target, source)
        end
    end
end

"""
    direct!(elements)

Assumes sources and targets are identical sets.

# Optional Arguments

- `reflex`- if `true`, includes interactions of each element on itself

"""
function direct!(elements::Vector; reflex=false)
    for (i,target) in enumerate(elements)
        for (j,source) in enumerate(elements)
            if i!=j || reflex; kernel!(target, source); end
        end
    end
end
