
function localized_particle_interaction!(systems::Tuple, radius; 
    leaf_size=default_leaf_size(systems), 
    is_target=SVector{length(systems), Bool}(true for _ in 1:length(systems)), 
    is_source=SVector{length(systems), Bool}(true for _ in 1:length(systems)))

    tree = Tree(systems; leaf_size, is_target, is_source)
    direct_list = build_localized_interaction_list(tree, radius)

    for i in eachindex(systems)
        !is_source[i] && continue
        source_system = systems[i]
        for j in eachindex(systems)
            !is_target[j] && continue
            target_system = systems[j]
            for interaction in direct_list
                t_leaf, s_branch = interaction
                target_idx = tree.branches[t_leaf].bodies_index[j]
                source_idx = tree.branches[s_branch].bodies_index[i]
                localized_direct!(target_system, target_idx, source_system, source_idx)
            end
        end
    end

    unsort!(systems, tree)
    return nothing
end

function build_localized_interaction_list(tree::Tree, radius)
    direct_list = Vector{SVector{2,Int32}}(undef,0)
    for i_leaf in tree.leaf_index
        if !tree.branches[i_leaf].target
            continue
        end
        build_localized_interaction_list!(direct_list, tree.branches, 1, i_leaf, tree, radius)
    end
    return direct_list
end

function build_localized_interaction_list!(list, branches, j_branch, i_leaf, tree, radius)
    current_branch = branches[j_branch]

    if is_outside_radius(tree.branches[i_leaf], current_branch, radius)
        return nothing
    end

    if current_branch.i_leaf >= 0
        push!(list, SVector{2,Int32}(i_leaf, j_branch))
        return nothing
    end

    if is_encapsulated(tree.branches[i_leaf], current_branch, radius)
        push!(list, SVector{2,Int32}(i_leaf, j_branch))
        return nothing
    end

    for k in current_branch.branch_index
        build_localized_interaction_list!(list, branches, k, i_leaf, tree, radius)
    end
    return nothing
end

# return true if branch is encapsulated by radius around leaf
function is_encapsulated(leaf, branch, radius)
    dx, dy, dz = maximum_edge_distance(leaf.target_center, leaf.target_box, branch.target_center, branch.target_box)
    return sqrt(dx^2 + dy^2 + dz^2) < radius
end

function is_outside_radius(leaf, branch, radius)
    dx, dy, dz = minimum_edge_distance(leaf.target_center, leaf.target_box, branch.target_center, branch.target_box)
    return sqrt(dx^2 + dy^2 + dz^2) > radius
end

function localized_direct!(target_system, target_index, source_system, source_index)
    @warn "localized_direct! not overloaded for type $(typeof(source_system)); interaction ignored"
    return nothing
end
