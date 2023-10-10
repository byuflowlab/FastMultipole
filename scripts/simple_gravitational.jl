include("../test/gravitational.jl")
using BenchmarkTools

n_bodies = 30000
bodies = rand(8,n_bodies)
system = Gravitational(bodies)
# options = fmm.Options(13,700,10.0)
expansion_order, n_per_branch, theta = 14, 1000, 0.34
old_bodies = deepcopy(system.bodies)
println("Create Single Tree:")
@time tree = fmm.Tree(system, expansion_order, n_per_branch)
# note: old_bodies[tree.index_list[1]] = systems[1].bodies
println("Run FMM:")
# @btime fmm.fmm!($tree, $systems, $options; unsort_bodies=true)
@time fmm.fmm!(tree, system; theta=theta, reset_tree=true, nearfield=true, farfield=true, unsort_bodies=true)
println("done.")
system2 = Gravitational(bodies)
println("Run direct:")
@time fmm.direct!(system2, 1:n_bodies, system2, 1:n_bodies)
println("done.")
phi = system.potential[1,:]
phi2 = system2.potential[1,:]
@show maximum(abs.(phi2 - phi))

# run some tests
system3 = Gravitational(bodies)
fmm.fmm!(system3; n_per_branch=n_per_branch, theta=0.34, nearfield=true, farfield=true, unsort_bodies=true)
@time fmm.fmm!(system3; n_per_branch=n_per_branch, theta=0.34, nearfield=true, farfield=true, unsort_bodies=true)
system4 = Gravitational(bodies)
fmm.fmm!((system4,); n_per_branch=n_per_branch, theta=0.34, nearfield=true, farfield=true, unsort_bodies=true)
@time fmm.fmm!((system4,); n_per_branch=n_per_branch, theta=0.34, nearfield=true, farfield=true, unsort_bodies=true)
system5 = Gravitational(bodies)
fmm.fmm!(system5, system5; n_per_branch_source=n_per_branch, n_per_branch_target=n_per_branch, theta=0.34, nearfield=true, farfield=true, unsort_source_bodies=true, unsort_target_bodies=true)
@time fmm.fmm!(system5, system5; n_per_branch_source=n_per_branch, n_per_branch_target=n_per_branch, theta=0.34, nearfield=true, farfield=true, unsort_source_bodies=true, unsort_target_bodies=true)
system6 = Gravitational(bodies)
fmm.fmm!((system6,), system6; n_per_branch_source=n_per_branch, n_per_branch_target=n_per_branch, theta=0.34, nearfield=true, farfield=true, unsort_source_bodies=true, unsort_target_bodies=true)
@time fmm.fmm!((system6,), system6; n_per_branch_source=n_per_branch, n_per_branch_target=n_per_branch, theta=0.34, nearfield=true, farfield=true, unsort_source_bodies=true, unsort_target_bodies=true)
system7 = Gravitational(bodies)
fmm.fmm!((system7,), (system7,); n_per_branch_source=n_per_branch, n_per_branch_target=n_per_branch, theta=0.34, nearfield=true, farfield=true, unsort_source_bodies=true, unsort_target_bodies=true)
@time fmm.fmm!((system7,), (system7,); n_per_branch_source=n_per_branch, n_per_branch_target=n_per_branch, theta=0.34, nearfield=true, farfield=true, unsort_source_bodies=true, unsort_target_bodies=true)
println("done.")

# n_bodies = 30000
# bodies = rand(8,n_bodies)
# systems = (Gravitational(bodies),)
# # systems = (fmm.SortWrapper(Gravitational(bodies)),)
# options = fmm.Options(13,700,10.0)
# old_bodies = deepcopy(systems[1].bodies)
# @time tree = fmm.Tree(systems, options)
# # note: old_bodies[tree.index_list[1]] = systems[1].bodies
# println("Run FMM:")
# # @btime fmm.fmm!($tree, $systems, $options; unsort_bodies=true)
# @time fmm.fmm!(tree, systems, options; unsort_bodies=true)
# println("done.")
# systems2 = (Gravitational(bodies),)
# println("Run direct:")
# @time fmm.direct!(systems2[1], 1:n_bodies, systems2[1], 1:n_bodies)
# println("done.")
# phi = systems[1].potential[1,:]
# phi2 = systems2[1].potential[1,:]
# @show maximum(abs.(phi2 - phi))
