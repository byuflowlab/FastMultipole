using FastMultipole

include("../test/gravitational.jl")
include("simple_gravitational.jl")

expansion_order = 4
multipole_threshold = 0.4
n_bodies = 1_000_000
leaf_size = 50
sys = generate_gravitational(123, n_bodies; radius_factor=0)
tree = Tree(sys; expansion_order, leaf_size, shrink_recenter=false)

# upward pass
FastMultipole.upward_pass_singlethread!(tree.branches, sys, tree.expansion_order)

# horizontal pass, existing functions
farfield, nearfield, self_induced = true, true, true
m2l_list, direct_target_bodies, direct_source_bodies = FastMultipole.build_interaction_lists(tree.branches, tree.branches, tree.leaf_index, multipole_threshold, farfield, nearfield, self_induced)

#@time FastMultipole.horizontal_pass_singlethread!(tree.branches, tree.branches, m2l_list, tree.expansion_order)
@time FastMultipole.horizontal_pass_singlethread!(tree.branches, tree.branches, m2l_list, tree.expansion_order)
@time FastMultipole.horizontal_pass_singlethread!(tree.branches, tree.branches, m2l_list, tree.expansion_order)

@profview FastMultipole.horizontal_pass_singlethread!(tree.branches, tree.branches, m2l_list, tree.expansion_order)

# using Profile
# using PProf
# @profile FastMultipole.horizontal_pass_singlethread!(tree.branches, tree.branches, m2l_list, tree.expansion_order)
# Profile.clear()
# @profile FastMultipole.horizontal_pass_singlethread!(tree.branches, tree.branches, m2l_list, tree.expansion_order)
# pprof()
