function P2P!(system, target_branch::SingleBranch, source_branch::SingleBranch)
    direct!(system, target_branch.bodies_index, system, source_branch.bodies_index)
end

function P2P!(systems, target_branch::MultiBranch, source_branch::MultiBranch)
    for (i,source_system) in enumerate(systems)
        source_indices = source_branch.bodies_index[i]
        for (j,target_system) in enumerate(systems)
            target_indices = target_branch.bodies_index[j]
            direct!(target_system, target_indices, source_system, source_indices)
        end
    end
end

function P2P!(target_system, target_branch::SingleBranch, source_system, source_branch::SingleBranch)
    direct!(target_system, target_branch.bodies_index, source_system, source_branch.bodies_index)
end

function P2P!(target_system, target_branch::SingleBranch, source_systems, source_branch::MultiBranch)
    target_indices = target_branch.bodies_index
    for (i,source_system) in enumerate(source_systems)
        source_indices = source_branch.bodies_index[i]
        direct!(target_system, target_indices, source_system, source_indices)
    end
end

function P2P!(target_systems, target_branch::MultiBranch, source_system, source_branch::SingleBranch)
    source_indices = source_branch.bodies_index
    for (j,target_system) in enumerate(target_systems)
        target_indices = target_branch.bodies_index[j]
        direct!(target_system, target_indices, source_system, source_indices)
    end
end

function P2P!(target_systems, target_branch::MultiBranch, source_systems, source_branch::MultiBranch)
    for (i,source_system) in enumerate(source_systems)
        source_indices = source_branch.bodies_index[i]
        for (j,target_system) in enumerate(target_systems)
            target_indices = target_branch.bodies_index[j]
            direct!(target_system, target_indices, source_system, source_indices)
        end
    end
end

#####
##### upward pass
#####
function upward_pass_single_thread!(branches, systems, expansion_order)
    # initialize memory
    harmonics = zeros(eltype(branches[1].multipole_expansion), (expansion_order+1)*(expansion_order+1))
    M = zeros(eltype(branches[1].multipole_expansion), 4)

    # loop over branches
    for branch in view(branches,length(branches):-1:2) # no need to create a multipole expansion at the very top level
        if branch.n_branches == 0 # branch is a leaf
            B2M!(branch, systems, harmonics, expansion_order)
        else # not a leaf
            # iterate over children
            for child_branch in view(branches, branch.branch_index)
                M2M!(branch, child_branch, harmonics, M, expansion_order)
            end
        end
    end
end

function body_2_multipole_multi_thread!(branches, systems, expansion_order, leaf_index)
    n_threads = Threads.nthreads()
    n_per_chunk, rem = divrem(length(leaf_index),n_threads)
    rem > 0 && (n_per_chunk += 1)

    Threads.@threads for i_start in 1:n_per_chunk:length(leaf_index)
        harmonics = zeros(eltype(branches[i_start].multipole_expansion), (expansion_order+1)*(expansion_order+1))
        M = zeros(eltype(branches[i_start].multipole_expansion), 4)
        for leaf in view(branches,i_start:i_start+n_per_chunk-1)
            B2M!(leaf, systems, harmonics, expansion_order)
        end
    end
end

function upward_pass_multithread!(branches, systems, expansion_order, levels_indices, leaf_index, thread_pool)
    # create multipole expansions
    body_2_multipole_multi_thread!(branches, systems, expansion_order, leaf_index, thread_pool)

    # m2m translation
    translate_multipoles_multi_thread!(branches, systems, expansion_order, levels_indices, thread_pool)
end

#####
##### horizontal pass
#####
function nearfield_single_thread!(system, branches, direct_list)
    for (i_target, j_source) in direct_list
        P2P!(system, branches[i_target], branches[j_source])
    end
end

function nearfield_single_thread!(target_system, target_branches, source_system, source_branches, direct_list)
    for (i_target, j_source) in direct_list
        P2P!(target_system, target_branches[i_target], source_system, source_branches[j_source])
    end
end

# function nearfield_multithread!(system, branches, direct_list)
#     # divide chunks
#     n_P2P = length(direct_list)
#     n_chunks = Threads.nthreads() - 1
#     n_per_chunk, rem = divrem(n_P2P, n_chunks)

#     # spread remainder
#     for i_start in 1:n_per_chunk+1:(n_per_chunk+1)*rem
#         @show i_start:i_start+n_per_chunk
#         for (i_target, j_source) in view(direct_list,i_start:i_start+n_per_chunk)
#             P2P!(system, branches[i_target], branches[j_source])
#         end
#     end
#     # regular chunks
#     for i_start in (n_per_chunk+1)*rem+1:n_per_chunk:n_P2P
#         @show i_start:i_start+n_per_chunk-1
#         for (i_target, j_source) in view(direct_list,i_start:i_start+n_per_chunk-1)
#             P2P!(system, branches[i_target], branches[j_source])
#         end
#     end
# end

# function nearfield_multithread!(target_system, target_branches, source_system, source_branches, direct_list)
#     # divide chunks
#     n_P2P = length(direct_list)
#     n_chunks = Threads.nthreads() - 1
#     n_per_chunk, rem = divrem(n_P2P, n_chunks)

#     # spread remainder
#     for i_start in 1:n_per_chunk+1:(n_per_chunk+1)*rem
#         @show i_start:i_start+n_per_chunk
#         for (i_target, j_source) in view(direct_list,i_start:i_start+n_per_chunk)
#             P2P!(target_system, target_branches[i_target], source_system, source_branches[j_source])
#         end
#     end
#     # regular chunks
#     for i_start in (n_per_chunk+1)*rem+1:n_per_chunk:n_P2P
#         @show i_start:i_start+n_per_chunk-1
#         for (i_target, j_source) in view(direct_list,i_start:i_start+n_per_chunk-1)
#             P2P!(target_system, target_branches[i_target], source_system, source_branches[j_source])
#         end
#     end
# end

function horizontal_pass_single_thread!(branches, m2l_list, expansion_order)
    harmonics = zeros(eltype(branches[1].multipole_expansion), (expansion_order<<1 + 1)*(expansion_order<<1 + 1))
    for (i_target, j_source) in m2l_list
        M2L!(branches[i_target], branches[j_source], harmonics, expansion_order)
    end
end

function horizontal_pass_single_thread!(target_branches, source_branches, m2l_list, expansion_order)
    harmonics = zeros(eltype(target_branches[1].multipole_expansion), (expansion_order<<1 + 1)*(expansion_order<<1 + 1))
    for (i_target, j_source) in m2l_list
        M2L!(target_branches[i_target], source_branches[j_source], harmonics, expansion_order)
    end
end

function horizontal_pass_multi_thread!(branches, m2l_list, expansion_order)
    # divide chunks
    n_M2L = length(m2l_list)
    n_chunks = Threads.nthreads() - 1
    n_per_chunk, rem = divrem(n_P2P, n_chunks)

    harmonics = zeros(eltype(branches[1].multipole_expansion), (expansion_order<<1 + 1)*(expansion_order<<1 + 1))
    for (i_target, j_source) in m2l_list
        M2L!(branches[i_target], branches[j_source], harmonics, expansion_order)
    end
end

function horizontal_pass_multi_thread!(target_branches, source_branches, m2l_list, expansion_order)
    harmonics = zeros(eltype(target_branches[1].multipole_expansion), (expansion_order<<1 + 1)*(expansion_order<<1 + 1))
    for (i_target, j_source) in m2l_list
        M2L!(target_branches[i_target], source_branches[j_source], harmonics, expansion_order)
    end
end

# function horizontal_pass_single_thread!(branches, m2l_list, expansion_order)
#     harmonics = zeros(eltype(branches[1].multipole_expansion), (expansion_order<<1 + 1)*(expansion_order<<1 + 1))
#     for (i_target, j_source) in m2l_list
#         M2L!(branches[i_target], branches[j_source], harmonics, expansion_order)
#     end
# end

# function horizontal_pass_single_thread!(target_branches, source_branches, m2l_list, expansion_order)
#     harmonics = zeros(eltype(target_branches[1].multipole_expansion), (expansion_order<<1 + 1)*(expansion_order<<1 + 1))
#     for (i_target, j_source) in m2l_list
#         M2L!(target_branches[i_target], source_branches[j_source], harmonics, expansion_order)
#     end
# end

#####
##### downward pass
#####
function downward_pass_single_thread!(branches, systems, expansion_order)
    regular_harmonics = zeros(eltype(branches[1].multipole_expansion), (expansion_order+1)*(expansion_order+1))
    vector_potential = zeros(eltype(branches[1]),3)
    potential_jacobian = zeros(eltype(branches[1]),3,4)
    potential_hessian = zeros(eltype(branches[1]),3,3,4)
    derivative_harmonics = zeros(eltype(branches[1].multipole_expansion), ((expansion_order+1) * (expansion_order+2)) >> 1)
    derivative_harmonics_theta = zeros(eltype(branches[1].multipole_expansion), ((expansion_order+1) * (expansion_order+2)) >> 1)
    derivative_harmonics_theta_2 = zeros(eltype(branches[1].multipole_expansion), ((expansion_order+1) * (expansion_order+2)) >> 1)
    workspace = zeros(eltype(branches[1]),3,4)
    L = zeros(eltype(branches[1].multipole_expansion),4)
    for branch in branches
        if branch.n_branches == 0 # leaf level
            L2B!(systems, branch, expansion_order, vector_potential, potential_jacobian, potential_hessian, derivative_harmonics, derivative_harmonics_theta, derivative_harmonics_theta_2, workspace)
        else
            for child_branch in view(branches,branch.branch_index)
                L2L!(branch, child_branch, regular_harmonics, L, expansion_order)
            end
        end
    end
end

# function downward_pass_multi_thread!(branches, systems, level_indices, expansion_order)
#     for level_index in level_indices
#         Threads.@threads for branch in view(branches,level_index)
#             harmonics = zeros(eltype(branch.multipole_expansion), (expansion_order+1)*(expansion_order+1))
#             if branch.n_branches == 0 # leaf level
#                 L2B!(systems, branch, expansion_order)
#             else
#                 for child_branch in view(branches,branch.branch_index)
#                     L2L!(branch, child_branch, harmonics, expansion_order)
#                 end
#             end
#         end
#     end
# end

#####
##### create interaction lists
#####
function build_interaction_lists(branches, theta, farfield, nearfield)
    m2l_list = Vector{SVector{2,Int32}}(undef,0)
    direct_list = Vector{SVector{2,Int32}}(undef,0)
    build_interaction_lists!(m2l_list, direct_list, 1, 1, branches, branches, theta, farfield, nearfield)
    
    return m2l_list, direct_list
end

# function build_interaction_lists!(m2l_list, direct_list, i_target, j_source, branches, theta, farfield, nearfield)
#     source_branch = branches[j_source]
#     target_branch = branches[i_target]

#     spacing = source_branch.center - target_branch.center
#     center_spacing_squared = spacing[1]*spacing[1] + spacing[2]*spacing[2] + spacing[3]*spacing[3]
#     summed_radii_squared = target_branch.radius + source_branch.radius
#     summed_radii_squared *= summed_radii_squared
#     if center_spacing_squared * theta * theta >= summed_radii_squared && farfield# meet M2L criteria
#         push!(m2l_list, SVector{2}(i_target, j_source))
#     elseif source_branch.n_branches == target_branch.n_branches == 0 && nearfield # both leaves
#         push!(direct_list, SVector{2}(i_target, j_source))
#     elseif source_branch.n_branches == 0 || (target_branch.radius >= source_branch.radius && target_branch.n_branches != 0)
#         for i_child in target_branch.branch_index
#             build_interaction_lists!(m2l_list, direct_list, i_child, j_source, branches, theta, farfield, nearfield)
#         end
#     else
#         for j_child in source_branch.branch_index
#             build_interaction_lists!(m2l_list, direct_list, i_target, j_child, branches, theta, farfield, nearfield)
#         end
#     end
# end

function build_interaction_lists(target_branches, source_branches, theta, farfield, nearfield)
    m2l_list = Vector{SVector{2,Int32}}(undef,0)
    direct_list = Vector{SVector{2,Int32}}(undef,0)
    build_interaction_lists!(m2l_list, direct_list, 1, 1, target_branches, source_branches, theta, farfield, nearfield)
    
    return m2l_list, direct_list
end

function build_interaction_lists!(m2l_list, direct_list, i_target, j_source, target_branches, source_branches, theta, farfield, nearfield)
    source_branch = source_branches[j_source]
    target_branch = target_branches[i_target]

    spacing = source_branch.center - target_branch.center
    center_spacing_squared = spacing[1]*spacing[1] + spacing[2]*spacing[2] + spacing[3]*spacing[3]
    summed_radii_squared = target_branch.radius + source_branch.radius
    summed_radii_squared *= summed_radii_squared
    if center_spacing_squared * theta * theta >= summed_radii_squared && farfield# meet M2L criteria
        push!(m2l_list, SVector{2}(i_target, j_source))
    elseif source_branch.n_branches == target_branch.n_branches == 0 && nearfield # both leaves
        push!(direct_list, SVector{2}(i_target, j_source))
    elseif source_branch.n_branches == 0 || (target_branch.radius >= source_branch.radius && target_branch.n_branches != 0)
        for i_child in target_branch.branch_index
            build_interaction_lists!(m2l_list, direct_list, i_child, j_source, target_branches, source_branches, theta, farfield, nearfield)
        end
    else
        for j_child in source_branch.branch_index
            build_interaction_lists!(m2l_list, direct_list, i_target, j_child, target_branches, source_branches, theta, farfield, nearfield)
        end
    end
end

#####
##### running FMM
#####
function fmm!(tree::Tree, systems; theta=0.4, reset_tree=true, nearfield=true, farfield=true, unsort_bodies=true)
    # reset multipole/local expansions
    reset_tree && (reset_expansions!(tree))

    # create interaction lists
    # println("interaction list")
    # @time m2l_list, direct_list = build_interaction_lists(tree.branches, theta, farfield, nearfield)
    m2l_list, direct_list = build_interaction_lists(tree.branches, theta, farfield, nearfield)
    
    # run FMM
    # println("nearfield")
    # @time nearfield && (nearfield!(systems, tree.branches, direct_list))
    nearfield && (nearfield_single_thread!(systems, tree.branches, direct_list))
    if farfield
        # println("upward pass")
        upward_pass_single_thread!(tree.branches, systems, tree.expansion_order)
        # @time upward_pass_single_thread!(tree.branches, systems, tree.expansion_order)
        # println("horizontal pass")
        # @time horizontal_pass!(tree.branches, m2l_list, tree.expansion_order)
        horizontal_pass_single_thread!(tree.branches, m2l_list, tree.expansion_order)
        # println("downward pass")
        downward_pass_single_thread!(tree.branches, systems, tree.expansion_order)
        # @time downward_pass_single_thread!(tree.branches, systems, tree.expansion_order)
    end
    
    # unsort bodies
    # println("unsort bodies")
    # @time unsort_bodies && (unsort!(systems, tree))
    unsort_bodies && (unsort!(systems, tree))
end

function fmm!(target_tree::Tree, target_systems, source_tree::Tree, source_systems; theta=0.4, reset_source_tree=true, reset_target_tree=true, nearfield=true, farfield=true, unsort_source_bodies=true, unsort_target_bodies=true)
    # reset multipole/local expansions
    reset_target_tree && (reset_expansions!(source_tree))
    reset_source_tree && (reset_expansions!(source_tree))

    # create interaction lists
    m2l_list, direct_list = build_interaction_lists(target_tree.branches, source_tree.branches, theta, farfield, nearfield)

    # run FMM
    nearfield && (nearfield_single_thread!(target_systems, target_tree.branches, source_systems, source_tree.branches, direct_list))
    if farfield
        upward_pass_single_thread!(source_tree.branches, source_systems, source_tree.expansion_order)
        horizontal_pass_single_thread!(target_tree.branches, source_tree.branches, m2l_list, source_tree.expansion_order)
        downward_pass_single_thread!(target_tree.branches, target_systems, target_tree.expansion_order)
    end

    # unsort bodies
    unsort_source_bodies && (unsort!(source_systems, source_tree))
    unsort_target_bodies && (unsort!(target_systems, target_tree))
end

"""
    fmm!(elements::Tuple{Element1, Element2,...}, options::Options)

Calculates the influence of all scalar and/or vector potential elements using the fast multipole method in spherical harmonics.

# Inputs

- `elements::Tuple{Element1, Element2, ...}`- a tuple of structs, each containing the following members:

    * `bodies::Array{Float64,2}`- 3+4+mxN array containing element positions, strengths, and m other values that must be sorted into the octree
    * `index::Vector{Int32}`- length N; `index[i]` contains the new index of the element that used to occupy the ith position before being sorted into the octree
    * `potential::Array{Float64,2}`- 4+12+36xN array of the potential, Jacobian, and Hessian that are reset every iteration (and don't require sorting)
    * `direct!::Function`- function calculates the direct influence of the body at the specified location
    * `B2M!::Function`- function converts the body's influence into a multipole expansion

Note: this function merely adds to existing potential of its elements to avoid overwriting potential calculated due to other sources.

The user must reset the potential manually.
"""
function fmm!(systems; expansion_order=5, n_per_branch=50, theta=0.4, ndivisions=7, nearfield=true, farfield=true, unsort_bodies=true, shrink_recenter=true, save_tree=false, save_name="tree")
    tree = Tree(systems; expansion_order=expansion_order, n_per_branch=n_per_branch, ndivisions=ndivisions, shrink_recenter=shrink_recenter)
    fmm!(tree, systems; theta=theta, reset_tree=false, nearfield=nearfield, farfield=farfield, unsort_bodies=unsort_bodies)
    save_tree && (visualize(save_name, systems, tree))
    return tree
end

function fmm!(target_systems, source_systems; expansion_order=5, n_per_branch_source=50, n_per_branch_target=50, theta=0.4, ndivisions_source=7, ndivisions_target=7, nearfield=true, farfield=true, unsort_source_bodies=true, unsort_target_bodies=true, source_shrink_recenter=true, target_shrink_recenter=true, save_tree=false, save_name="tree")
    source_tree = Tree(source_systems; expansion_order=expansion_order, n_per_branch=n_per_branch_source, shrink_recenter=source_shrink_recenter, ndivisions=ndivisions_source)
    target_tree = Tree(target_systems; expansion_order=expansion_order, n_per_branch=n_per_branch_target, shrink_recenter=target_shrink_recenter, ndivisions=ndivisions_target)
    fmm!(target_tree, target_systems, source_tree, source_systems; theta=theta, reset_source_tree=false, reset_target_tree=false, nearfield=nearfield, farfield=farfield, unsort_source_bodies=unsort_source_bodies, unsort_target_bodies=unsort_target_bodies)
    save_tree && (visualize(save_name, systems, tree))
    return source_tree, target_tree
end
