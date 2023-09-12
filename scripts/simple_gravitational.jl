include("../test/gravitational.jl")
using BenchmarkTools

n_bodies = 300000
bodies = rand(8,n_bodies)
systems = (Gravitational(bodies),)
options = fmm.Options(13,700,10.0)
old_bodies = deepcopy(systems[1].bodies)
tree = fmm.Tree(systems, options)
# note: old_bodies[tree.index_list[1]] = systems[1].bodies
println("Run FMM:")
@btime fmm.fmm!($tree, $systems, $options; unsort_bodies=true)
println("done.")
systems2 = (Gravitational(bodies),)
println("Run direct:")
@btime fmm.direct!(systems2[1], 1:n_bodies, systems2[1], 1:n_bodies)
println("done.")
phi = systems[1].potential[1,:]
phi2 = systems2[1].potential[1,:]
@show maximum(abs.(phi2 - phi))
