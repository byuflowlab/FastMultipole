@testset "dynamic expansion order: absolute rotated coefficients, point source" begin

expansion_order, leaf_size_source, multipole_threshold = 20, SVector{1}(10), 0.5
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

ε = 1e-5
ε_tol = FastMultipole.RotatedCoefficientsAbsoluteVectorField(ε, false)
# ε_tol = nothing
system = generate_gravitational(seed, n_bodies; radius_factor=0.0)
system2 = generate_gravitational(seed, n_bodies; radius_factor=0.1)

# println("\n===== radius factor = 0.0 =====\n")

vector_field_null = system.potential[5:7,:]
optimized_args, cache, target_tree, source_tree, m2l_list, direct_list, derivatives_switches, error_success = FastMultipole.fmm!(system; expansion_order, leaf_size_source, multipole_threshold, nearfield=true, farfield=true, shrink_recenter, ε_tol)

vector_field_fmm = system.potential[5:7,:]
vector_field_err = [norm(system.potential[5:7,i] - validation_system.potential[5:7,i]) for i in 1:size(system.potential,2)]

@test ε * 0.1 < maximum(vector_field_err) < ε * 30

# println("\n===== radius factor = 0.1 =====\n")

FastMultipole.fmm!(system2; expansion_order, leaf_size_source, multipole_threshold, nearfield=true, farfield=true, shrink_recenter, ε_tol)

vector_field_err = [norm(system2.potential[5:7,i] - validation_system.potential[5:7,i]) for i in 1:size(system.potential,2)]

@test ε * 0.1 < maximum(vector_field_err) < ε * 30

end

@testset "dynamic expansion order: absolute rotated coefficients, point vortex" begin

expansion_order, leaf_size_source, multipole_threshold = 20, SVector{1}(10), 0.5
n_bodies = 10000

shrink_recenter = true
seed = 123
validation_system = generate_vortex(seed, n_bodies; radius_factor=0.0)
FastMultipole.direct!(validation_system)
validation_potential = validation_system.vector_field_stretching[1:3,:]

validation_system2 = generate_vortex(seed, n_bodies; radius_factor=0.1)
FastMultipole.direct!(validation_system2)
validation_potential2 = validation_system2.vector_field_stretching[1:3,:]

@assert validation_potential == validation_potential2

ε = 1e-5
ε_tol = FastMultipole.RotatedCoefficientsAbsoluteVectorField(ε, false)
lamb_helmholtz = true
# ε_tol = nothing
system = generate_vortex(seed, n_bodies; radius_factor=0.0)
system2 = generate_vortex(seed, n_bodies; radius_factor=0.1)

# println("\n===== radius factor = 0.0 =====\n")

FastMultipole.fmm!(system; expansion_order, leaf_size_source, multipole_threshold, lamb_helmholtz, shrink_recenter, ε_tol)

vector_field_err = [norm(system.vector_field_stretching[1:3,i] - validation_system.vector_field_stretching[1:3,i]) for i in 1:size(system.vector_field_stretching,2)]

@test ε * 0.1 < maximum(vector_field_err) < ε * 80

# println("\n===== radius factor = 0.1 =====\n")

FastMultipole.fmm!(system2; expansion_order, leaf_size_source, multipole_threshold, lamb_helmholtz, shrink_recenter, ε_tol)

vector_field_err = [norm(system2.vector_field_stretching[1:3,i] - validation_system.vector_field_stretching[1:3,i]) for i in 1:size(system.vector_field_stretching,2)]

@test ε * 0.1 < maximum(vector_field_err) < ε * 80

end

@testset "dynamic expansion order: absolute multipole power, point source" begin

expansion_order, leaf_size_source, multipole_threshold = 20, SVector{1}(10), 0.5
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

ε = 1e-5
ε_tol = FastMultipole.PowerAbsoluteVectorField(ε, false)
system = generate_gravitational(seed, n_bodies; radius_factor=0.0)
system2 = generate_gravitational(seed, n_bodies; radius_factor=0.1)

# println("\n===== radius factor = 0.0 =====\n")

FastMultipole.fmm!(system; expansion_order, leaf_size_source, multipole_threshold, nearfield=true, farfield=true, shrink_recenter, ε_tol)

vector_field_err = [norm(system.potential[5:7,i] - validation_system.potential[5:7,i]) for i in 1:size(system.potential,2)]

@test ε * 0.1 < maximum(vector_field_err) < ε * 10

# println("\n===== radius factor = 0.1 =====\n")

FastMultipole.fmm!(system2; expansion_order, leaf_size_source, multipole_threshold, nearfield=true, farfield=true, shrink_recenter, ε_tol)

vector_field_err = [norm(system2.potential[5:7,i] - validation_system.potential[5:7,i]) for i in 1:size(system.potential,2)]

@test ε * 0.1 < maximum(vector_field_err) < ε * 10

end

@testset "dynamic expansion order: absolute multipole power, point vortex" begin

expansion_order, leaf_size_source, multipole_threshold = 20, SVector{1}(10), 0.5
n_bodies = 10000
lamb_helmholtz = true

shrink_recenter = true
seed = 123
validation_system = generate_vortex(seed, n_bodies; radius_factor=0.0)
FastMultipole.direct!(validation_system)
validation_potential = validation_system.potential

validation_system2 = generate_vortex(seed, n_bodies; radius_factor=0.1)
FastMultipole.direct!(validation_system2)
validation_potential2 = validation_system2.potential

@assert validation_potential == validation_potential2

ε = 1e-5
ε_tol = FastMultipole.PowerAbsoluteVectorField(ε, false)
# ε_tol = nothing
system = generate_vortex(seed, n_bodies; radius_factor=0.0)
system2 = generate_vortex(seed, n_bodies; radius_factor=0.1)

# println("\n===== radius factor = 0.0 =====\n")

tree, m2l_list, direct_list, derivatives_switches = FastMultipole.fmm!(system; expansion_order, leaf_size_source, multipole_threshold, lamb_helmholtz, shrink_recenter, ε_tol)

vector_field_err = [norm(system.vector_field_stretching[1:3,i] - validation_system.vector_field_stretching[1:3,i]) for i in 1:size(system.potential,2)]

@test ε * 0.1 < maximum(vector_field_err) < ε * 10

# println("\n===== radius factor = 0.1 =====\n")

tree2, m2l_list2, direct_list2, derivatives_switches2 = FastMultipole.fmm!(system2; expansion_order, leaf_size_source, multipole_threshold, lamb_helmholtz, shrink_recenter, ε_tol)

vector_field_err = [norm(system2.vector_field_stretching[1:3,i] - validation_system.vector_field_stretching[1:3,i]) for i in 1:size(system.potential,2)]

@test ε * 0.1 < maximum(vector_field_err) < ε * 10

end




# @testset "dynamic expansion order: influence estimate, point source" begin

# expansion_order, leaf_size_source, multipole_threshold = 20, SVector{1}(100), 0.5
# n_bodies = 10000

# shrink_recenter = true
# seed = 123
# validation_system = generate_gravitational(seed, n_bodies; radius_factor=0.0)
# FastMultipole.direct!(validation_system)
# validation_potential = validation_system.potential[1,:]
# validation_vector = validation_system.potential[5:7,:]

# validation_system2 = generate_gravitational(seed, n_bodies; radius_factor=0.1)
# fmm!(validation_system2;
#     scalar_potential = true, vector_field = true, vector_gradient = false,
#     leaf_size_source = FastMultipole.to_vector(5, 1),
#     expansion_order = 3, multipole_threshold = 0.6,
#     ε_tol = nothing, shrink_recenter = true, nearfield_device = false,
#     update_target_systems = true,
#     silence_warnings = true
# )
# validation_potential2 = validation_system2.potential[1,:]
# validation_vector2 = validation_system2.potential[5:7,:]

# vector_field = [norm(validation_vector_field[1:3,i]) for i in 1:size(validation_vector_field,2)]
# vector_field_err = [norm(validation_vector_field[1:3,i] - validation_vector_field2[1:3,i]) for i in 1:size(validation_vector_field2,2)]
# relative_err = vector_field_err ./ vector_field

# @test maximum(relative_err) < 1e-1

# end







# @testset "dynamic expansion order: relative rotated coefficients, point source" begin

# expansion_order, leaf_size_source, multipole_threshold = 20, SVector{1}(100), 0.5
# n_bodies = 10000

# shrink_recenter = true
# seed = 123
# validation_system = generate_gravitational(seed, n_bodies; radius_factor=0.0)
# FastMultipole.direct!(validation_system)
# validation_potential = validation_system.potential[1,:]

# validation_system2 = generate_gravitational(seed, n_bodies; radius_factor=0.1)
# FastMultipole.direct!(validation_system2)
# validation_potential2 = validation_system2.potential[1,:]

# @assert validation_potential == validation_potential2

# ε = 1e-4
# ε_tol = FastMultipole.RotatedCoefficientsRelativeVectorField(ε, false)
# # ε_tol = nothing
# system = generate_gravitational(seed, n_bodies; radius_factor=0.0)
# system2 = generate_gravitational(seed, n_bodies; radius_factor=0.1)

# # println("\n===== radius factor = 0.0 =====\n")

# FastMultipole.DEBUG[] = true
# _, _, target_tree, _ = FastMultipole.fmm!(system; expansion_order, leaf_size_source, multipole_threshold, nearfield=true, farfield=true, shrink_recenter, ε_tol)
# FastMultipole.DEBUG[] = false

# vector_field = [norm(validation_system.potential[5:7,i]) for i in 1:size(validation_system.potential,2)]
# vector_field_err = [norm(system.potential[5:7,i] - validation_system.potential[5:7,i]) for i in 1:size(system.potential,2)]
# relative_err = vector_field_err ./ vector_field

# i_worst = findfirst((x) -> x == maximum(relative_err), relative_err)
# @show i_worst FastMultipole.unsorted_index_2_sorted_index(i_worst, 1, target_tree)

# @show mean(relative_err) std(relative_err) maximum(relative_err)
# @show mean(vector_field_err) std(vector_field_err) maximum(vector_field_err)
# @show mean(vector_field) std(vector_field) maximum(vector_field)
# @show sum(relative_err .> ε) / length(relative_err)

# @test ε * 0.1 < maximum(relative_err) < ε * 10

# # println("\n===== radius factor = 0.1 =====\n")

# FastMultipole.fmm!(system2; expansion_order, leaf_size_source, multipole_threshold, nearfield=true, farfield=true, shrink_recenter, ε_tol)

# vector_field = [norm(system2.potential[5:7,:]) for i in 1:size(system2.potential,2)]
# vector_field_err = [norm(system2.potential[5:7,i] - validation_system.potential[5:7,i]) for i in 1:size(system2.potential,2)]
# relative_err = vector_field_err ./ vector_field

# @test ε * 0.1 < maximum(relative_err) < ε * 10

# end

# @testset "dynamic expansion order: relative rotated coefficients, point vortex" begin

# expansion_order, leaf_size_source, multipole_threshold = 20, SVector{1}(100), 0.5
# n_bodies = 10000

# shrink_recenter = true
# seed = 123
# validation_system = generate_vortex(seed, n_bodies; radius_factor=0.0)
# FastMultipole.direct!(validation_system)
# validation_potential = validation_system.vector_field_stretching[1:3,:]

# validation_system2 = generate_vortex(seed, n_bodies; radius_factor=0.1)
# FastMultipole.direct!(validation_system2)
# validation_potential2 = validation_system2.vector_field_stretching[1:3,:]

# @assert validation_potential == validation_potential2

# ε = 1e-4
# ε_tol = FastMultipole.RotatedCoefficientsRelativeVectorField(ε, false)
# lamb_helmholtz = true
# # ε_tol = nothing
# system = generate_vortex(seed, n_bodies; radius_factor=0.0)
# system2 = generate_vortex(seed, n_bodies; radius_factor=0.1)

# # println("\n===== radius factor = 0.0 =====\n")

# FastMultipole.fmm!(system; expansion_order, leaf_size_source, multipole_threshold, lamb_helmholtz, shrink_recenter, ε_tol)

# vector_field = [norm(system.vector_field_stretching[1:3,i]) for i in 1:size(system.potential,2)]
# vector_field_err = [norm(system.vector_field_stretching[1:3,i] - validation_system.vector_field_stretching[1:3,i]) for i in 1:size(system.potential,2)]
# relative_err = vector_field_err ./ vector_field

# @test ε * 0.1 < maximum(vector_field_err) < ε * 60

# # println("\n===== radius factor = 0.1 =====\n")

# FastMultipole.fmm!(system2; expansion_order, leaf_size_source, multipole_threshold, lamb_helmholtz, shrink_recenter, ε_tol)

# vector_field = [norm(system2.vector_field_stretching[1:3,i]) for i in 1:size(system2.potential,2)]
# vector_field_err = [norm(system2.vector_field_stretching[1:3,i] - validation_system.vector_field_stretching[1:3,i]) for i in 1:size(system2.potential,2)]
# relative_err = vector_field_err ./ vector_field

# @test ε * 0.1 < maximum(vector_field_err) < ε * 60

# end

# @testset "dynamic expansion order: relative multipole power, point source" begin

# expansion_order, leaf_size_source, multipole_threshold = 20, SVector{1}(100), 0.5
# n_bodies = 10000

# shrink_recenter = true
# seed = 123
# validation_system = generate_gravitational(seed, n_bodies; radius_factor=0.0)
# FastMultipole.direct!(validation_system)
# validation_potential = validation_system.potential[1,:]

# validation_system2 = generate_gravitational(seed, n_bodies; radius_factor=0.1)
# FastMultipole.direct!(validation_system2)
# validation_potential2 = validation_system2.potential[1,:]

# @assert validation_potential == validation_potential2

# ε = 1e-4
# ε_tol = FastMultipole.PowerRelativeVectorField(ε, false)
# system = generate_gravitational(seed, n_bodies; radius_factor=0.0)
# system2 = generate_gravitational(seed, n_bodies; radius_factor=0.1)

# # println("\n===== radius factor = 0.0 =====\n")

# FastMultipole.fmm!(system; expansion_order, leaf_size_source, multipole_threshold, nearfield=true, farfield=true, shrink_recenter, ε_tol)

# vector_field = [norm(system.potential[5:7,:]) for i in 1:size(system.potential,2)]
# vector_field_err = [norm(system.potential[5:7,i] - validation_system.potential[5:7,i]) for i in 1:size(system.potential,2)]
# relative_err = vector_field_err ./ vector_field

# @test ε * 0.1 < maximum(vector_field_err) < ε * 10

# # println("\n===== radius factor = 0.1 =====\n")

# FastMultipole.fmm!(system2; expansion_order, leaf_size_source, multipole_threshold, nearfield=true, farfield=true, shrink_recenter, ε_tol)

# vector_field = [norm(system2.potential[5:7,:]) for i in 1:size(system2.potential,2)]
# vector_field_err = [norm(system2.potential[5:7,i] - validation_system.potential[5:7,i]) for i in 1:size(system2.potential,2)]
# relative_err = vector_field_err ./ vector_field

# @test ε * 0.1 < maximum(relative_err) < ε * 10

# end

# @testset "dynamic expansion order: relative multipole power, point vortex" begin

# expansion_order, leaf_size_source, multipole_threshold = 20, SVector{1}(100), 0.5
# n_bodies = 10000
# lamb_helmholtz = true

# shrink_recenter = true
# seed = 123
# validation_system = generate_vortex(seed, n_bodies; radius_factor=0.0)
# FastMultipole.direct!(validation_system)
# validation_potential = validation_system.potential

# validation_system2 = generate_vortex(seed, n_bodies; radius_factor=0.1)
# FastMultipole.direct!(validation_system2)
# validation_potential2 = validation_system2.potential

# @assert validation_potential == validation_potential2

# ε = 1e-4
# ε_tol = FastMultipole.PowerRelativeVectorField(ε, false)
# system = generate_vortex(seed, n_bodies; radius_factor=0.0)
# system2 = generate_vortex(seed, n_bodies; radius_factor=0.1)

# # println("\n===== radius factor = 0.0 =====\n")

# FastMultipole.fmm!(system; expansion_order, leaf_size_source, multipole_threshold, lamb_helmholtz, shrink_recenter, ε_tol)

# vector_field = [norm(system.vector_field_stretching[1:3,i]) for i in 1:size(system.potential,2)]
# vector_field_err = [norm(system.vector_field_stretching[1:3,i] - validation_system.vector_field_stretching[1:3,i]) for i in 1:size(system.potential,2)]
# relative_err = vector_field_err ./ vector_field

# @test ε * 0.1 < maximum(vector_field_err) < ε * 10

# # println("\n===== radius factor = 0.1 =====\n")

# FastMultipole.fmm!(system2; expansion_order, leaf_size_source, multipole_threshold, lamb_helmholtz, shrink_recenter, ε_tol)

# vector_field = [norm(system2.vector_field_stretching[1:3,i]) for i in 1:size(system2.potential,2)]
# vector_field_err = [norm(system2.vector_field_stretching[1:3,i] - validation_system.vector_field_stretching[1:3,i]) for i in 1:size(system2.potential,2)]
# relative_err = vector_field_err ./ vector_field

# @test ε * 0.1 < maximum(vector_field_err) < ε * 30

# end
