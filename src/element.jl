abstract type Element end

struct Source <: Element end

struct Dipole <: Element end

abstract type Dimension end

struct Point <: Dimension end

struct Line <: Dimension end

struct Panel <: Dimension end

struct QuadPanel <: Dimension end

struct Volume <: Dimension end

abstract type Distribution end

struct Uniform <: Distribution end

struct Linear <: Distribution end


#####
##### code architecture
#####
# FMM

#     * inputs: 
        
#         - provide a tuple of vectors of elements (bodies)

#             - position (3-tuple for point, 5 3-tuples for quad panels, etc.)
#             - strength (scalar for <:Uniform, vector for <:Linear )

#         - element dimension
#         - element type
#         - function for computing multipole coefficients? (at least make it optional)

# function fmm!(bodies::Vector, potential, options::Options, targets, sources; reset_tree, local_P2P)

# end
    
#     * outputs:

#         - potential, gradient, and gradient derivatives
#         - sorts elements into octree (make copies?)
#         - allow AD derivatives
#         - calculates potential, gradient, and gradient derivatives

# FLOWSolve

#     * user:

#         - provide a tuple of vectors of boundary elements
#         - provide a tuple of vectors of wake elements
#         - somehow describe boundary conditions
    
#     * code:

#         - solve for the strengths of each boundary element
