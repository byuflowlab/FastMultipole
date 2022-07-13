#####
##### abstract Element definition
#####
"For each element introduced, the following functions must be defined:

* `get_X(::Element)` - returns the element location

"
abstract type Element{TF} end

##
## Example Point <: Element
##
struct Point{TF} <: Element{TF}
    X::Vector{TF} # locations
end

function get_X(point::Point)
    return point.X
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

function eltype(elements::AbstractArray{e}) where e<:Element
    return eltype(get_X(elements,1))
end

# function iterate(elements::Vector{e}) where e <: Element
#     return get_X(elements,1), 1
# end

# function iterate(elements::Vector{e}, state) where e <: Element
#     if state > length(elements)
#         return nothing
#     else
#         return get_X(elements,state+1), state+1
#     end
# end

##
## Next steps: implement sorting using sort! to build tree; update Element struct usage in the rest of the code; test FMM with multiple particles; add M2M and L2L functions
##
