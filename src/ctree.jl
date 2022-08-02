abstract type Element end

struct Branch
    n_branches          # number of child branches
    n_elements          # number of descendent elements
    first_branch        # index of the first branch
    first_element       # index of the first element
    center              # center of the branch
    radius              # side lengths of the rectangle encapsulating the branch
    multipole_expansion # multipole expansion coefficients
    local_expansion     # local expansion coefficients
end

struct Tree
    indices         # a vector of integer indices referring to the associated element list
    branches        # a vector of `Branch` objects composing the tree
    expansion_order # order of the multipole/local expansions
    n_per_branch    # max number of elements in a leaf
end

function Tree(elements; expansion_order=8, n_per_branch=1)
    # initialize objects
    n_elements = length(elements)
    indices = collect(1:n_elements)
    buffer = similar(indices)
    branches = Vector{Branch}(undef,1)

    # recursively build branches
    i_start = 1
    i_end = n_elements
    i_branch = 1
    center, radius = center_radius(elements; scale_radius = 1.00001)
    level = 0
    branch!(branches, indices, buffer, elements, i_start, i_end, i_branch, center, radius, level, expansion_order, n_per_branch)

    # assemble tree
    tree = Tree(indices, branches, expansion_order, n_per_branch)

    return tree
end

function branch!(branches, indices, buffer, elements, i_start, i_end, i_branch, center, radius, level, expansion_order, n_per_branch)
    n_branches = 0
    n_elements = i_end - i_start + 1
    multipole_expansion = zeros(expansion_order)
    local_expansion = zeros(expansion_order)
    if n_elements <= n_per_branch # branch is a leaf; no new branches needed
        i_child = -1
        branch = Branch(n_branches, n_elements, i_child, i_start, center, radius, multipole_expansion, local_expansion)
        branches[i_branch] = branch
        return nothing
    else # not a leaf; branch children
        # count elements in each octant
        octant_attendance = zeros(Int64,8)
        for i_element in indices[i_start:i_end]
            x = get_x(elements[i_element])
            i_octant = get_octant(x, center)
            # @show x center (i_octant-0b1) indices
            octant_attendance[i_octant] += 1
        end

        # determine number of children
        for i_octant in 1:8
            n_branches += octant_attendance[i_octant] > 0
        end

        # create branch
        i_child = length(branches) + 1
        branch = Branch(n_branches, n_elements, i_child, i_start, center, radius, multipole_expansion, local_expansion)
        branches[i_branch] = branch

        # sort element indices
        ## write offsets
        offsets = zeros(Int64,8)
        offsets[1] = i_start
        for i=2:8
            offsets[i] = offsets[i-1] + octant_attendance[i-1]
        end

        ## create counter
        counter = deepcopy(offsets)

        ## sort element indices into the buffer
        for i_element in indices[i_start:i_end]
            x = get_x(elements[i_element])
            i_octant = get_octant(x, center)
            buffer[counter[i_octant]] = i_element
            counter[i_octant] += 1
        end
        # println("Sherlock!")
        # @show indices
        indices[i_start:i_end] .= buffer[i_start:i_end]

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
                branch!(branches, indices, buffer, elements, child_i_start, child_i_end, i_tape + prev_branches, child_center, child_radius, level + 1, expansion_order, n_per_branch)
                i_tape += 1
            end
        end
    end
end

function center_radius(elements::AbstractArray{e}; scale_radius = 1.00001) where e<:Element
    x_min = deepcopy(get_x(elements[1]))
    x_max = deepcopy(get_x(elements[1]))
    for element in elements
        x = get_x(element)
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

#####
##### test script
#####

# build list of elements to sort
struct Mass{TF} <: Element
    x::Array{TF,1}
    mass::Array{TF,1}
    potential::Array{TF,1}
    force::Array{TF,1}
end

function get_x(mass::Mass)
    mass.x
end

@inline function get_octant(x, center)
    i_octant = (UInt8((x[1] > center[1])) + UInt8(x[2] > center[2]) << 0b1 + UInt8(x[3] > center[3]) << 0b10) + 0b1
end

# n_masses = 10
# masses = Vector{Mass}(undef,n_masses)
# for i in 1:length(masses)
#     x = 2 .* rand(3) .- 1.0
#     mass = rand(1)
#     potential = zeros(1)
#     force = zeros(3)
#     masses[i] = Mass(x,mass,potential,force)
# end

xs = [
    1.2 1.1 0.8;
    0.8 0.9 0.2;
    0.1 0.2 0.9;
    0.1 0.3 0.2;
    0.2 0.25 0.4
]
ms = rand(size(xs)[1])
masses = Vector{Mass}(undef,length(ms))
for i in 1:length(ms)
    x = xs[i,:]
    mass = rand(1)
    potential = zeros(1)
    force = zeros(3)
    masses[i] = Mass(x,mass,potential,force)
end

# test center_radius function
center, radius = center_radius(masses; scale_radius = 1.00001)
test_center = [0.65, 0.65, 0.55]
test_radius = 0.5500055

for i in 1:3
    @show isapprox(center[i], test_center[i]; atol=1e-4)
end
@show isapprox(radius, test_radius; atol=1e-4)

# test branch! function
tree = Tree(masses)

test_branches = [
    5 0.65 0.65 0.55 0.5500055;
    2 0.37499725 0.37499725 0.27499725 0.27500275;
    1 0.92500275 0.92500275 0.27499725 0.27500275;
    1 0.37499725 0.37499725 0.82500275 0.27500275;
    1 0.92500275 0.92500275 0.82500275 0.27500275;
    1 0.237495875 0.237495875 0.137495875 0.137501375;
    1 0.237495875 0.237495875 0.412498625 0.137501375;
]

@show length(tree.branches) == size(test_branches)[1]

for i_branch in 1:length(tree.branches)
    @show isapprox(tree.branches[i_branch].n_elements, test_branches[i_branch,1]; atol=1e-8)
    for i in 1:3
        @show isapprox(tree.branches[i_branch].center[i], test_branches[i_branch,1+i]; atol=1e-7)
    end
    @show isapprox(tree.branches[i_branch].radius, test_branches[i_branch,5]; atol=1e-7)
end
