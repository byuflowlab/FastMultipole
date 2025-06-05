#--- define influence function ---#

function FastMultipole.influence!(influence, target_buffer, ::Gravitational, source_buffer)
    influence .= view(target_buffer, 4, :)
end

function FastMultipole.value_to_strength!(source_buffer, ::Gravitational, i_body, value)
    source_buffer[5, i_body] = value
end

@testset "Fast Gauss Seidel: self influence" begin

#--- create system ---#

n_bodies = 10
seed = 1234
system = generate_gravitational(seed, n_bodies)

direct!(system; scalar_potential=true, velocity=false)
phi_desired = system.potential[1, :]

#--- create FGS solver ---#

fgs = FastMultipole.FastGaussSeidel((system,), (system,); expansion_order=4, multipole_threshold=0.5, leaf_size=30)

#--- check influence matrix ---#

influence_matrix_check = zeros(Float64, n_bodies, n_bodies)
source_buffer = fgs.source_tree.buffers[1]
for i in 1:n_bodies
    for j in 1:n_bodies
        r = norm(FastMultipole.get_position(source_buffer, i) .- FastMultipole.get_position(source_buffer, j))
        if r > 0.0
            influence_matrix_check[i, j] = 1.0 / (r*4*pi)
        end
    end
end

influence_matrix_fgs, _ = FastMultipole.get_matrix_vector(fgs.self_matrices, 1)
for i in 1:n_bodies
    for j in 1:n_bodies
        @test influence_matrix_fgs[i, j] ≈ influence_matrix_check[i, j] atol=1e-6
    end
end

end

@testset "Fast Gauss Seidel: all influence matrices" begin

#--- larger system ---#

n_bodies = 1000
seed = 1234
system = generate_gravitational(seed, n_bodies)

direct!(system; scalar_potential=true, velocity=false)
phi_desired = system.potential[1, :]

#--- create FGS solver ---#

fgs = FastMultipole.FastGaussSeidel((system,), (system,); expansion_order=4, multipole_threshold=0.5, leaf_size=5) # try with leaf_size=3 for sources with no non-self influence

#--- check self influence matrices ---#

for (i, i_leaf) in enumerate(fgs.source_tree.leaf_index)

    bodies_index = fgs.source_tree.branches[i_leaf].bodies_index[1]
    n_bodies = length(bodies_index)
    influence_matrix_check = zeros(Float64, n_bodies, n_bodies)
    source_buffer = fgs.source_tree.buffers[1]
    for (i,i_body) in enumerate(bodies_index)
        for (j,j_body) in enumerate(bodies_index)
            r = norm(FastMultipole.get_position(source_buffer, i_body) .- FastMultipole.get_position(source_buffer, j_body))
            if r > 0.0
                influence_matrix_check[i, j] = 1.0 / (r*4*pi)
            end
        end
    end

    influence_matrix_fgs, _ = FastMultipole.get_matrix_vector(fgs.self_matrices, i)
    for (i,i_body) in enumerate(bodies_index)
        for (j,j_body) in enumerate(bodies_index)
            @test influence_matrix_fgs[i, j] ≈ influence_matrix_check[i, j] atol=1e-6
        end
    end

end

#--- check non-self influence matrices ---#

direct_list = fgs.direct_list
target_buffer = fgs.target_tree.buffers[1]
source_buffer = fgs.source_tree.buffers[1]
source_list = [dl[2] for dl in direct_list]

for i_leaf in 1:length(fgs.source_tree.leaf_index)
    this_source = fgs.source_tree.leaf_index[i_leaf]
    if this_source in source_list
        mat, rhs = FastMultipole.get_matrix_vector(fgs.nonself_matrices, i_leaf)

        # build influence matrix manually to check against FGS
        source_indices = fgs.source_tree.branches[this_source].bodies_index[1]
        target_indices = Int[]
        for (i_target, j_source) in direct_list
            if j_source == this_source
                target_bodies_index = fgs.target_tree.branches[i_target].bodies_index[1]
                target_indices = vcat(target_indices, collect(target_bodies_index))
            end
        end

        # create test influence matrix
        influence_matrix_check = zeros(Float64, length(target_indices), length(source_indices))
        for (i_target, i_target_body) in enumerate(target_indices)
            for (j_source, j_source_body) in enumerate(source_indices)
                r = norm(FastMultipole.get_position(source_buffer, j_source_body) .- FastMultipole.get_position(target_buffer, i_target_body))
                if r > 0.0
                    influence_matrix_check[i_target, j_source] = 1.0 / (r*4*pi)
                end
            end
        end

        # test
        for (i_target, i_target_body) in enumerate(target_indices)
            for (j_source, j_source_body) in enumerate(source_indices)
                @test mat[i_target, j_source] ≈ influence_matrix_check[i_target, j_source] atol=1e-6
            end
        end
    else
        @assert fgs.nonself_matrices.sizes[i_leaf][1] == 0 "this influence matrix should be empty"
    end
end

end