#####
##### kernel related functions
#####
"""
    kernel!(target, source)

Describes the interaction of a source element on a target element.

# Inputs

- `target`- a user-defined element for which the functions `FLOWFMM.get_x` and `FLOWFMM.get_q` have been defined.

"""
function kernel! end

"""
    get_x(element)

Returns the location of an element.

# Inputs

- `element`- a user-defined object for which the functions `FLOWFMM.kernel!` and `FLOWFMM.get_q` have been defined.

Outputs:

- `x::Vector{Float64}`- a vector of length 3 containing the Cartesian coordinates of the element's location

"""
function get_x end

"""
get_q(element)

Returns the strength of an element.

# Inputs

- `element`- a user-defined object for which the functions `FLOWFMM.kernel!` and `FLOWFMM.get_x` have been defined.

Outputs:

- `Float64`- the elements strength

"""
function get_q end

function derivatives end
