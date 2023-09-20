struct Body end
const BODY = Body()

struct Position end
const POSITION = Position()

struct Radius end
const RADIUS = Radius()

struct Potential end
const POTENTIAL = Potential()

struct ScalarPotential end
const SCALAR_POTENTIAL = ScalarPotential()

struct VectorPotential end
const VECTOR_POTENTIAL = VectorPotential()

struct Velocity end
const VELOCITY = Velocity()

struct VelocityGradient end
const VELOCITYGRADIENT = VelocityGradient()

struct Branch{TF,N}
    n_branches::Int8        # number of child branches
    n_bodies::SVector{N,Int32}         # number of descendent bodies
    first_branch::Int32     # index of the first branch
    first_body::SVector{N,Int32}       # index of the first element
    center::SVector{3,TF}   # center of the branch
    radius::TF              # side lengths of the cube encapsulating the branch
    multipole_expansion::NTuple{4,Vector{Complex{TF}}} # multipole expansion coefficients
    local_expansion::NTuple{4,Vector{Complex{TF}}}     # local expansion coefficients
    lock::ReentrantLock
    child_lock::ReentrantLock
end

"""
bodies[index_list] is the same sort operation as performed by the tree
sorted_bodies[inverse_index_list] undoes the sort operation performed by the tree
"""
struct Tree{TF,N}
    branches::Vector{Branch{TF,N}}        # a vector of `Branch` objects composing the tree
    expansion_order::Int16
    n_per_branch::Int32    # max number of bodies in a leaf
    index_list::NTuple{N,Vector{Int}}
    inverse_index_list::NTuple{N,Vector{Int}}
    leaf_index::Vector{Int}
    cumulative_count::Vector{Int} # starting with 0, a cumulative accounting of how many bodies are in leaf branches
end

struct Val{T} end