struct Branch
    n_branches          # number of child branches
    n_bodies          # number of descendent bodies
    first_branch        # index of the first branch
    first_body       # index of the first element
    center              # center of the branch
    radius              # side lengths of the rectangle encapsulating the branch
    multipole_expansion # multipole expansion coefficients
    local_expansion     # local expansion coefficients
end

struct Tree
    branches        # a vector of `Branch` objects composing the tree
    expansion_order
    n_per_branch    # max number of bodies in a leaf
end

"""
    Tree(elements; expansion_order=2, n_per_branch=1)

Constructs an octree of the provided element objects.

# Inputs

- `elements`- a struct containing the following members:

    * `bodies::Array{Float64,2}`- 3+4+mxN array containing element positions, strengths, and m other values that must be sorted into the octree
    * `potential::Array{Float64,2}`- 4+12+36xN array of the potential, Jacobian, and Hessian that are reset every iteration (and don't require sorting)
    * `velocity::Array{Float64,2}`- 3xN array of the velocity vectors at each element, reset every iteration, and calculated in post-processing
    * `direct!::Function`- function calculates the direct influence of the body at the specified location
    * `B2M!::Function`- function converts the body's influence into a multipole expansion
"""
function Tree(elements; expansion_order=2, n_per_branch=1)
    # initialize objects
    bodies = elements.bodies
    buffer = similar(bodies)
    branches = Vector{Branch}(undef,1)

    # recursively build branches
    i_start = 1
    i_end = size(bodies)[2]
    i_branch = 1
    center, radius = center_radius(elements; scale_radius = 1.00001)
    level = 0
    branch!(branches, bodies, buffer, i_start, i_end, i_branch, center, radius, level, expansion_order, n_per_branch)

    # assemble tree
    tree = Tree(branches, [expansion_order], n_per_branch)

    return tree
end

function branch!(branches, bodies, buffer, i_start, i_end, i_branch, center, radius, level, expansion_order, n_per_branch)
    n_branches = 0
    n_bodies = i_end - i_start + 1
    multipole_expansion = initialize_expansion(expansion_order)
    local_expansion = initialize_expansion(expansion_order)
    if n_bodies <= n_per_branch # branch is a leaf; no new branches needed
        i_child = -1
        branch = Branch(n_branches, n_bodies, i_child, i_start, center, radius, multipole_expansion, local_expansion)
        branches[i_branch] = branch
        return nothing
    else # not a leaf; branch children
        # count elements in each octant
        octant_attendance = zeros(Int64,8)
        for i_body in i_start:i_end
            x = bodies[1:3,i_body]
            i_octant = get_octant(x, center)
            octant_attendance[i_octant] += 1
        end

        # determine number of children
        for i_octant in 1:8
            n_branches += octant_attendance[i_octant] > 0
        end

        # create branch
        i_child = length(branches) + 1
        branch = Branch(n_branches, n_bodies, i_child, i_start, center, radius, multipole_expansion, local_expansion)
        branches[i_branch] = branch

        # sort bodies
        ## write offsets
        offsets = zeros(Int64,8)
        offsets[1] = i_start
        for i=2:8
            offsets[i] = offsets[i-1] + octant_attendance[i-1]
        end

        ## create counter
        counter = deepcopy(offsets)

        ## sort element indices into the buffer
        for i_body in i_start:i_end
            x = bodies[1:3,i_body]
            i_octant = get_octant(x, center)
            buffer[:,counter[i_octant]] .= bodies[:,i_body]
            bodies[1,i_body] += 45
            if buffer[1,i_body] == bodies[1,i_body]
                println("Warning: reference error")
            end
            bodies[1,i_body] -= 45
            counter[i_octant] += 1
        end

        bodies[:,i_start:i_end] .= buffer[:,i_start:i_end]

        # recursively build new branches
        prev_branches = length(branches)
        resize!(branches, prev_branches + n_branches)

        i_tape = 1
        for i_octant in Int8(0):Int8(7)
            if octant_attendance[i_octant+Int8(1)] > 0
                child_i_start = offsets[i_octant+Int8(1)]
                child_i_end = child_i_start + octant_attendance[i_octant+Int8(1)] - 1
                child_radius = radius / 2#(1 << (level + 1))
                child_center = deepcopy(center)
                for d in Int8(0):Int8(2)
                    child_center[d + Int8(1)] += child_radius * (((i_octant & Int8(1) << d) >> d) * Int8(2) - Int8(1))
                end
                branch!(branches, bodies, buffer, child_i_start, child_i_end, i_tape + prev_branches, child_center, child_radius, level + 1, expansion_order, n_per_branch)
                i_tape += 1
            end
        end
    end
end

function center_radius(elements; scale_radius = 1.00001)
    x_min = deepcopy(elements.bodies[i_POSITION,1])
    x_max = deepcopy(elements.bodies[i_POSITION,1])
    for i in 1:size(elements.bodies)[2]
        x = elements.bodies[i_POSITION,i]
        for dim in 1:3
            if x[dim] < x_min[dim]
                x_min[dim] = x[dim]
            elseif x[dim] > x_max[dim]
                x_max[dim] = x[dim]
            end
        end
    end
    center = (x_max + x_min) ./ 2
    # TODO: add element smoothing radius here? Note this method creates cubic cells
    radius = max(x_max-center...) * scale_radius # get half of the longest side length of the rectangle
    return center, radius
end

@inline function get_octant(x, center)
    i_octant = (UInt8((x[1] > center[1])) + UInt8(x[2] > center[2]) << 0b1 + UInt8(x[3] > center[3]) << 0b10) + 0b1
end

function n_terms(expansion_order, dimensions)
    n = 0
    for order in 0:expansion_order
        # leverage stars and bars theorem
        n += binomial(order + dimensions - 1, dimensions - 1)
    end
    return n
end

function initialize_expansion(expansion_order)
    return [zeros(Complex{Float64}, ((expansion_order+1) * (expansion_order+2)) >> 1) for _ in 1:4]
end

function change_expansion_order!(tree::Tree, new_order)
    old_n = length(tree.branches[1].multipole_expansion[1])
    new_n = n_terms(new_order, 3)
    old_order = tree.expansion_order[1]
    for branch in tree.branches
        if  old_order < new_order
            for i in 1:new_n - old_n
                for dim in 1:4
                    push!(branch.multipole_expansion[dim],0.0)
                    push!(branch.local_expansion[dim],0.0)
                end
            end
        elseif old_order > new_order
            for i in 1:new_n - old_n
                for dim in 1:4
                    pop!(branch.multipole_expansion[dim])
                    pop!(branch.local_expansion[dim])
                end
            end
        end
    end
    tree.expansion_order[1] = new_order
end

function reset_expansions!(tree)
    for branch in tree.branches
        branch.multipole_expansion .*= 0
        branch.local_expansion .*= 0
    end
end
