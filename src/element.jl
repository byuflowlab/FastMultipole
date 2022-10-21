#####
##### abstract Element definition
#####
"""
    type ::Element

Instances of ::Element contain the following properties:

- `strength`: a vector of length 4 whose first index describes the source potential, and whose final 3 indices represent the vector source
- `potential`: a vector of length 4 whose first index is the scalar potential, and whose final 3 indices represent a vector potential
- `velocity`: a vector of length 3 describing the induced velocity, obtained as the negative gradient of the scalar potential
- `jacobian`: a 3x3 matrix containing the spatial Jacobian of the vector potential
- `hessian`: a 3x3x3 array containing the spatial Hessian of the vector potential
"""
abstract type Element end

# abstract type VectorSource <: Element end

# abstract type ScalarSource <: Element end

# abstract type MultiSource <: Element end
