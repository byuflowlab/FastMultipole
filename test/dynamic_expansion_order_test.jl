using Statistics

expansion_order, leaf_size, multipole_threshold = 10, 100, 0.5
n_bodies = 10000

shrink_recenter = true
seed = 123
validation_system = generate_gravitational(seed, n_bodies; radius_factor=0.0)
FastMultipole.direct!(validation_system)
validation_potential = validation_system.potential[1,:]

validation_system2 = generate_gravitational(seed, n_bodies; radius_factor=0.1)
FastMultipole.direct!(validation_system2)
validation_potential2 = validation_system2.potential[1,:]

@assert validation_potential == validation_potential2

ε_tol = 1e-3
# ε_tol = nothing
system = generate_gravitational(seed, n_bodies; radius_factor=0.0)
system2 = generate_gravitational(seed, n_bodies; radius_factor=0.1)

println("\n===== radius factor = 0.0 =====\n")

tree, m2l_list, direct_list, derivatives_switches = FastMultipole.fmm!(system; expansion_order, leaf_size, multipole_threshold, nearfield=true, farfield=true, unsort_bodies=true, shrink_recenter, ε_tol)

@show length(m2l_list) length(direct_list)
potential = system.potential[1,:]
@show maximum(abs.(potential - validation_potential))

println("\n===== radius factor = 0.1 =====\n")

tree2, m2l_list2, direct_list2, derivatives_switches2 = FastMultipole.fmm!(system2; expansion_order, leaf_size, multipole_threshold, nearfield=true, farfield=true, unsort_bodies=true, shrink_recenter, ε_tol)

@show length(m2l_list2) length(direct_list2)
potential2 = system2.potential[1,:]
@show maximum(abs.(potential2 - validation_potential))

#--- get true potential for nearfield = false ---#

# system2.potential .= 0.0
#tree3, m2l_list3, direct_list3, derivatives_switches3 = FastMultipole.fmm!(system2; expansion_order, leaf_size, multipole_threshold, nearfield=true, farfield=false, unsort_bodies=true, shrink_recenter, ε_tol)
# potential_nearfield_only = system2.potential[1,:]
# potential_farfield_only = validation_potential - potential_nearfield_only

# @test isapprox(maximum(abs.(potential - validation_potential)), 0.0; atol=1e-10)

#=
tree, m2l_list, direct_list, derivatives_switches = FastMultipole.fmm!(system2; expansion_order, leaf_size, multipole_threshold, nearfield=false, farfield=true, unsort_bodies=false, shrink_recenter, ε_tol)

# save potential from farfield only
potential_fmm = system2.potential[1,:]

# zero potential
system2.potential .= 0.0

lamb_helmholtz = Val(false)
εs_obs, εs_pred = FastMultipole.horizontal_pass_debug!(system2, tree.branches, system2, tree.branches, m2l_list, lamb_helmholtz, expansion_order, ε_tol)

potential_hp = system2.potential[1,:]

# zero potential
system2.potential .= 0.0
εs_obs, εs_pred = FastMultipole.horizontal_pass_debug!(system2, tree.branches, system2, tree.branches, m2l_list, lamb_helmholtz, expansion_order, ε_tol; store="multipole")

potential_multipole = system2.potential[1,:]

# zero potential
system2.potential .= 0.0
εs_obs, εs_pred = FastMultipole.horizontal_pass_debug!(system2, tree.branches, system2, tree.branches, m2l_list, lamb_helmholtz, expansion_order, ε_tol; store="direct")

potential_true = system2.potential[1,:]

function output_stuff(v)
    @show mean(v)
    @show std(v)
    @show minimum(v)
    @show maximum(v)
end

println()
println("observed / predicted:")
output_stuff(εs_obs ./ εs_pred)

@show maximum(εs_obs)

dp_true_mp = potential_true - potential_multipole
dp_true_fmm = potential_true - potential_fmm
dp_true_hp = potential_true - potential_hp

@show mean(dp_true_mp)
@show mean(dp_true_fmm)
@show mean(dp_true_hp)
=#
