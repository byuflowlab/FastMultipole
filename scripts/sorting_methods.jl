include("../test/gravitational.jl")

#####
##### sorting in-place
#####

# check tree
n_bodies = 10
bodies = rand(8,n_bodies)
systems = (Gravitational(bodies),)
old_bodies = deepcopy(systems[1].bodies)
options = fmm.Options(13,1,10.0)
tree = fmm.Tree(systems, options)
fmm.fmm!(tree,systems,options; unsort_bodies=true)
potential_fmm = systems[1].potential[1,:]

systems_direct = (Gravitational(bodies),)
fmm.direct!(systems_direct)
potential_direct = systems_direct[1].potential[1,:]

@show maximum(potential_direct - potential_fmm)

#####
##### sorting not in-place
#####
systems2 = (fmm.SortWrapper(Gravitational(bodies)),)
old_bodies2 = deepcopy(systems2[1].system.bodies)
options2 = fmm.Options(13,1,10.0)
tree2 = fmm.Tree(systems2, options2)
fmm.fmm!(tree2,systems2,options2; unsort_bodies=true)
potential_fmm2 = systems2[1].system.potential[1,:]

# check if branches match
for i in 1:length(tree.branches)
    @show tree.branches[i].center - tree2.branches[i].center
end

@show maximum(abs.(potential_fmm-potential_direct)) maximum(abs.(potential_fmm2-potential_direct))