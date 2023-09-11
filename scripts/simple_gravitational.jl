include("../test/gravitational.jl")

n_bodies = 30000
bodies = rand(8,n_bodies)
systems = (Gravitational(bodies),)
options = fmm.Options(13,700,10.0)
tree = fmm.Tree(systems, options)
println("Run FMM:")
@time fmm.fmm!(tree, systems, options; unsort_bodies=true)
println("done.")
systems2 = (Gravitational(bodies),)
println("Run direct:")
@time fmm.direct!(systems2[1], 1:n_bodies, systems2[1], 1:n_bodies)
println("done.")
phi = systems[1].potential[1,:]
phi2 = systems2[1].potential[1,:]
@show maximum(abs.(phi2 - phi))
