#####
##### direct interactions
#####
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
    # divide chunks
    n_threads = Threads.nthreads()
    n_per_chunk, rem = divrem(length(leaf_index),n_threads)
    if n_per_chunk == 0
        n_per_chunk = 1 # 1 task per thread, though not all threads get a task
        rem = 0 # no reason to spread a remainder
    end

    # spread remainder across rem chunks
    Threads.@threads for i_start in (range(1,step=n_per_chunk+1,length=rem)...,range(1+(n_per_chunk+1)*rem,step=n_per_chunk,stop=length(leaf_index))...)
        harmonics = zeros(eltype(branches[1].multipole_expansion), (expansion_order+1)*(expansion_order+1))
        chunk_size = i_start > (n_per_chunk+1)*rem ? n_per_chunk : n_per_chunk + 1
        for leaf in view(branches,view(leaf_index,i_start:i_start+chunk_size-1))
            B2M!(leaf, systems, harmonics, expansion_order)
        end
    end
end

function body_2_multipole_single_thread!(branches, systems, expansion_order, leaf_index)
    # divide chunks
    n_threads = Threads.nthreads()
    n_per_chunk, rem = divrem(length(leaf_index),n_threads)
    if n_per_chunk == 0
        n_per_chunk = 1 # 1 task per thread, though not all threads get a task
        rem = 0 # no reason to spread a remainder
    end

    # spread remainder across rem chunks
    for i_start in (range(1,step=n_per_chunk+1,length=rem)...,range(1+(n_per_chunk+1)*rem,step=n_per_chunk,stop=length(leaf_index))...)
        harmonics = zeros(eltype(branches[1].multipole_expansion), (expansion_order+1)*(expansion_order+1))
        chunk_size = i_start > (n_per_chunk+1)*rem ? n_per_chunk : n_per_chunk + 1
        for leaf in view(branches,view(leaf_index,i_start:i_start+chunk_size-1))
            B2M!(leaf, systems, harmonics, expansion_order)
        end
    end
end

@inline function translate_multipoles(branch, branches, harmonics, M, expansion_order)
    if branch.n_branches > 0 # branch is not a leaf

        # iterate over children
        for child_branch in view(branches, branch.branch_index)
            M2M!(branch, child_branch, harmonics, M, expansion_order)
        end
    end
end

function translate_multipoles_multi_thread!(branches, expansion_order, levels_index)
    # initialize memory
    n_threads = Threads.nthreads()
    
    # iterate over levels
    for level_index in view(levels_index,length(levels_index):-1:1)
        
        # if its too fine, no sense multithreading
        if Threads.nthreads() < 5 # length(level_index) > 100 # this just isn't efficient for high numbers of threads

            # divide chunks
            n_per_chunk, rem = divrem(length(level_index),n_threads)
            if n_per_chunk == 0
                n_per_chunk = 1 # 1 task per thread, though not all threads get a task
                rem = 0 # no reason to spread a remainder
            end

            # loop over branches
            Threads.@threads for i_start in (range(1,step=n_per_chunk+1,length=rem)...,range(1+(n_per_chunk+1)*rem,step=n_per_chunk,stop=length(level_index))...)
                chunk_size = i_start > (n_per_chunk+1)*rem ? n_per_chunk : n_per_chunk + 1

                # initialize memory
                harmonics = zeros(eltype(branches[1].multipole_expansion), (expansion_order+1)*(expansion_order+1))
                M = zeros(eltype(branches[1].multipole_expansion), 4)

                # loop over branches
                for branch in view(branches,view(level_index,i_start:i_start+chunk_size-1)) # no need to create a multipole expansion at the very top level
                    translate_multipoles(branch, branches, harmonics, M, expansion_order)
                end
            end
        else
            harmonics = zeros(eltype(branches[1].multipole_expansion), (expansion_order+1)*(expansion_order+1))
            M = zeros(eltype(branches[1].multipole_expansion), 4)
            for branch in view(branches,level_index)
                translate_multipoles(branch, branches, harmonics, M, expansion_order)
            end
        end
    end
    return nothing
end

function translate_multipoles_single_thread!(branches, expansion_order, levels_index)
    # initialize memory
    n_threads = Threads.nthreads()
    
    # iterate over levels
    for level_index in view(levels_index,length(levels_index):-1:1)
        harmonics = zeros(eltype(branches[1].multipole_expansion), (expansion_order+1)*(expansion_order+1))
        M = zeros(eltype(branches[1].multipole_expansion), 4)
        for branch in view(branches,level_index)
            translate_multipoles(branch, branches, harmonics, M, expansion_order)
        end
    end
    return nothing
end

function upward_pass_multi_thread!(branches, systems, expansion_order, levels_index, leaf_index)
    # create multipole expansions
    println("b2m")
    @time body_2_multipole_multi_thread!(branches, systems, expansion_order, leaf_index)

    # m2m translation
    println("tm")
    @time translate_multipoles_multi_thread!(branches, expansion_order, levels_index)
end

#####
##### horizontal pass
#####
function nearfield_single_thread!(target_system, target_branches, source_system, source_branches, direct_list)
    for (i_target, j_source) in direct_list
        P2P!(target_system, target_branches[i_target], source_system, source_branches[j_source])
    end
end

function nearfield_multi_thread!(target_system, target_branches, source_system, source_branches, direct_list)
    # divide chunks
    n_threads = Threads.nthreads()
    n_per_chunk, rem = divrem(length(direct_list),n_threads)
    if n_per_chunk == 0
        n_per_chunk = 1 # 1 task per thread, though not all threads get a task
        rem = 0 # no reason to spread a remainder
    end

    # spread remainder over rem chunks
    Threads.@threads for i_start in (range(1,step=n_per_chunk+1,length=rem)..., range(1+(n_per_chunk+1)*rem,step=n_per_chunk,stop=length(direct_list))...)
        chunk_size = i_start > (n_per_chunk+1)*rem ? n_per_chunk : n_per_chunk + 1
        for (i_target, j_source) in view(direct_list,i_start:i_start+chunk_size-1)
            P2P!(target_system, target_branches[i_target], source_system, source_branches[j_source])
        end
    end
end

function horizontal_pass_single_thread!(target_branches, source_branches, m2l_list, expansion_order)
    harmonics = zeros(eltype(target_branches[1].multipole_expansion), (expansion_order<<1 + 1)*(expansion_order<<1 + 1))
    L = zeros(eltype(target_branches[1].local_expansion), 4)
    for (i_target, j_source) in m2l_list
        M2L!(target_branches[i_target], source_branches[j_source], harmonics, L, expansion_order)
    end
end

function horizontal_pass_multi_thread!(target_branches, source_branches, m2l_list, expansion_order)
    # divide chunks
    n_threads = Threads.nthreads()
    n_per_chunk, rem = divrem(length(m2l_list),n_threads)
    if n_per_chunk == 0
        n_per_chunk = 1
        rem = 0
    end

    # spread remainder over rem chunks
    Threads.@threads for i_start in (range(1,step=n_per_chunk+1,length=rem)...,range(1+(n_per_chunk+1)*rem,step=n_per_chunk,stop=length(m2l_list))...)
        chunk_size = i_start > (n_per_chunk+1)*rem ? n_per_chunk : n_per_chunk + 1
        harmonics = zeros(eltype(target_branches[1].multipole_expansion), (expansion_order<<1 + 1)*(expansion_order<<1 + 1))
        L = zeros(eltype(target_branches[1].local_expansion), 4)
        for (i_target, j_source) in view(m2l_list,i_start:i_start+chunk_size-1)
            target_branch = target_branches[i_target]
            lock(target_branch.lock)
            try
                M2L!(target_branch, source_branches[j_source], harmonics, L, expansion_order)
            finally
                unlock(target_branch.lock)
            end
        end
    end
end

#####
##### downward pass
#####
function preallocate_l2b(float_type, expansion_type, expansion_order)
    vector_potential = zeros(float_type,3)
    potential_jacobian = zeros(float_type,3,4)
    potential_hessian = zeros(float_type,3,3,4)
    derivative_harmonics = zeros(expansion_type, ((expansion_order+1) * (expansion_order+2)) >> 1)
    derivative_harmonics_theta = zeros(expansion_type, ((expansion_order+1) * (expansion_order+2)) >> 1)
    derivative_harmonics_theta_2 = zeros(expansion_type, ((expansion_order+1) * (expansion_order+2)) >> 1)
    workspace = zeros(float_type,3,4)
    return vector_potential, potential_jacobian, potential_hessian, derivative_harmonics, derivative_harmonics_theta, derivative_harmonics_theta_2, workspace
end

function downward_pass_single_thread!(branches, systems, expansion_order)
    regular_harmonics = zeros(eltype(branches[1].multipole_expansion), (expansion_order+1)*(expansion_order+1))
    L = zeros(eltype(branches[1].multipole_expansion),4)
    vector_potential, potential_jacobian, potential_hessian, derivative_harmonics, derivative_harmonics_theta, derivative_harmonics_theta_2, workspace = preallocate_l2b(eltype(branches[1]), eltype(branches[1].multipole_expansion), expansion_order)
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

@inline function translate_locals(branch, branches, harmonics, L, expansion_order)
    if branch.n_branches > 0 # branch is not a leaf

        # iterate over children
        for child_branch in view(branches, branch.branch_index)
            L2L!(branch, child_branch, harmonics, L, expansion_order)
        end
    end
end

function translate_locals_multi_thread!(branches, expansion_order, levels_index)
    # initialize memory
    n_threads = Threads.nthreads()
    
    # iterate over levels
    for level_index in levels_index
        
        # if its too fine, no sense multithreading
        if Threads.nthreads() < 5 # length(level_index) > 100

            # divide chunks
            n_per_chunk, rem = divrem(length(level_index),n_threads)
            if n_per_chunk == 0
                n_per_chunk = 1 # 1 task per thread, though not all threads get a task
                rem = 0 # no reason to spread a remainder
            end

            # loop over branches
            Threads.@threads for i_start in (range(1,step=n_per_chunk+1,length=rem)...,range(1+(n_per_chunk+1)*rem,step=n_per_chunk,stop=length(level_index))...)
                chunk_size = i_start > (n_per_chunk+1)*rem ? n_per_chunk : n_per_chunk + 1 # some chunks are 1 branch longer to evenly spread the remainder in `n_per_chunk, rem = divrem(...)`

                # initialize memory
                harmonics = zeros(eltype(branches[1].multipole_expansion), (expansion_order+1)*(expansion_order+1))
                L = zeros(eltype(branches[1].multipole_expansion), 4)

                # loop over branches
                for branch in view(branches,view(level_index,i_start:i_start+chunk_size-1)) # no need to create a multipole expansion at the very top level
                    translate_locals(branch, branches, harmonics, L, expansion_order)
                end
            end
        else
            harmonics = zeros(eltype(branches[1].multipole_expansion), (expansion_order+1)*(expansion_order+1))
            L = zeros(eltype(branches[1].multipole_expansion), 4)
            for branch in view(branches,level_index)
                translate_locals(branch, branches, harmonics, L, expansion_order)
            end
        end
    end
    return nothing
end

function local_2_body_multi_thread!(branches, systems, expansion_order, leaf_index)
    # divide chunks
    n_threads = Threads.nthreads()
    n_per_chunk, rem = divrem(length(leaf_index),n_threads)
    if n_per_chunk == 0
        n_per_chunk = 1 # 1 task per thread, though not all threads get a task
        rem = 0 # no reason to spread a remainder
    end

    # spread remainder across rem chunks
    Threads.@threads for i_start in (range(1,step=n_per_chunk+1,length=rem)...,range(1+(n_per_chunk+1)*rem,step=n_per_chunk,stop=length(leaf_index))...)
        vector_potential, potential_jacobian, potential_hessian, derivative_harmonics, derivative_harmonics_theta, derivative_harmonics_theta_2, workspace = preallocate_l2b(eltype(branches[1]), eltype(branches[1].multipole_expansion), expansion_order)
        chunk_size = i_start > (n_per_chunk+1)*rem ? n_per_chunk : n_per_chunk + 1
        for leaf in view(branches,view(leaf_index,i_start:i_start+chunk_size-1))
            L2B!(systems, leaf, expansion_order, vector_potential, potential_jacobian, potential_hessian, derivative_harmonics, derivative_harmonics_theta, derivative_harmonics_theta_2, workspace)
        end
    end
end

function local_2_body_single_thread!(branches, systems, expansion_order, leaf_index)
    # divide chunks
    n_threads = Threads.nthreads()
    n_per_chunk, rem = divrem(length(leaf_index),n_threads)
    if n_per_chunk == 0
        n_per_chunk = 1 # 1 task per thread, though not all threads get a task
        rem = 0 # no reason to spread a remainder
    end

    # spread remainder across rem chunks
    for i_start in (range(1,step=n_per_chunk+1,length=rem)...,range(1+(n_per_chunk+1)*rem,step=n_per_chunk,stop=length(leaf_index))...)
        vector_potential, potential_jacobian, potential_hessian, derivative_harmonics, derivative_harmonics_theta, derivative_harmonics_theta_2, workspace = preallocate_l2b(eltype(branches[1]), eltype(branches[1].multipole_expansion), expansion_order)
        chunk_size = i_start > (n_per_chunk+1)*rem ? n_per_chunk : n_per_chunk + 1
        for leaf in view(branches,view(leaf_index,i_start:i_start+chunk_size-1))
            L2B!(systems, leaf, expansion_order, vector_potential, potential_jacobian, potential_hessian, derivative_harmonics, derivative_harmonics_theta, derivative_harmonics_theta_2, workspace)
        end
    end
end

function downward_pass_multi_thread!(branches, systems, expansion_order, levels_index, leaf_index)
    # m2m translation
    println("tl")
    @time translate_locals_multi_thread!(branches, expansion_order, levels_index)
    
    # create multipole expansions
    println("l2b")
    @time local_2_body_multi_thread!(branches, systems, expansion_order, leaf_index)
end

#####
##### create interaction lists
#####
function build_interaction_lists(branches, theta, farfield, nearfield)
    return build_interaction_lists(branches, branches, theta, farfield, nearfield)
end

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
    if center_spacing_squared * theta * theta >= summed_radii_squared && farfield # meet M2L criteria
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
    fmm!(tree, systems, tree, systems; theta, reset_source_tree=reset_tree, reset_target_tree=false, nearfield=nearfield, farfield=farfield, unsort_source_bodies=unsort_bodies, unsort_target_bodies=false)
end

function fmm!(target_tree::Tree, target_systems, source_tree::Tree, source_systems; theta=0.4, reset_source_tree=true, reset_target_tree=true, nearfield=true, farfield=true, unsort_source_bodies=true, unsort_target_bodies=true)
    # reset multipole/local expansions
    reset_target_tree && (reset_expansions!(source_tree))
    reset_source_tree && (reset_expansions!(source_tree))

    # create interaction lists
    m2l_list, direct_list = build_interaction_lists(target_tree.branches, source_tree.branches, theta, farfield, nearfield)

    # run FMM
    if Threads.nthreads() == 1
        println("nearfield")
        @time nearfield && (nearfield_single_thread!(target_systems, target_tree.branches, source_systems, source_tree.branches, direct_list))
        if farfield
            println("upward pass:")
            @time upward_pass_single_thread!(source_tree.branches, source_systems, source_tree.expansion_order)
            println("horizontal pass:")
            @time horizontal_pass_single_thread!(target_tree.branches, source_tree.branches, m2l_list, source_tree.expansion_order)
            println("downward pass:")
            @time downward_pass_single_thread!(target_tree.branches, target_systems, target_tree.expansion_order)
        end
    else # multithread
        println("nearfield")
        @time nearfield && (nearfield_multi_thread!(target_systems, target_tree.branches, source_systems, source_tree.branches, direct_list))
        if farfield
            upward_pass_multi_thread!(source_tree.branches, source_systems, source_tree.expansion_order, source_tree.levels_index, source_tree.leaf_index)
            println("horizontal pass:")
            @time horizontal_pass_multi_thread!(target_tree.branches, source_tree.branches, m2l_list, target_tree.expansion_order)
            downward_pass_multi_thread!(target_tree.branches, target_systems, target_tree.expansion_order, target_tree.levels_index, target_tree.leaf_index)
        end
    end

    # unsort bodies
    unsort_target_bodies && (unsort!(target_systems, target_tree))
    unsort_source_bodies && (unsort!(source_systems, source_tree))
end

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
