#####
##### abstract Element definition
#####
abstract type Element{TF} end

abstract type VectorPotential{TF} <: Element{TF} end

abstract type ScalarPotential{TF} <: Element{TF} end

##
## Example Point <: Element
##
struct ScalarPoint{TF} <: ScalarPotential{TF}
    V::Vector{TF} # scalar potential
    X::Vector{TF} # locations
end

function get_X(point::ScalarPoint)
    return point.X
end

function get_V(point::ScalarPoint)
    return point.V[1]
end

struct VectorPoint{TF} <: VectorPotential{TF}
    V::Vector{TF} # Vector potential
    X::Vector{TF} # locations
end

function get_X(point::VectorPoint)
    return point.X
end

function get_V(point::VectorPoint)
    return point.V
end

##
## general functions; perhaps I should remove these
##
function get_X(elements::AbstractArray{e}, i) where e<:Element
    return get_X(elements[i])
end

function get_q(elements::AbstractArray{e}, i) where e<:Element
    return get_q(elements[i])
end

function get_V(elements::AbstractArray{e}, i) where e<:Element
    return get_V(elements[i])
end

function get_dims(elements::AbstractArray{e}) where e<:Element
    return length(get_X(elements, 1))
end
