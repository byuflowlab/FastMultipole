import FLOWFMM as fmm
using SparseArrays
using Test

include("../test/gravitational.jl")

downloads_path = "/Users/randerson/Downloads"
notebook_path = "/Users/randerson/Dropbox/research/notebooks/"

positions = [
    0.0 0.20 1.0 0.94
    0.0 0.03 1.0 0.93
    0.0 0.14 1.0 0.84
]
strengths = [1.5, 2.1, 0.4, 1.8]
n_bodies = length(strengths)
strengths = vcat(strengths', zeros(3,length(strengths)))

element = Gravitational(vcat(positions, strengths))
elements_tuple = (element,)

save_vtk("test", element)

options = fmm.Options(2,2,8.0)
tree = fmm.Tree(elements_tuple, options)


S = spzeros(n_bodies, n_bodies)
old_strengths = deepcopy(element.bodies[4,:])
element.bodies[4,:] .= 1.0
for (i_leaf,leaf) in enumerate(tree.branches[tree.leaf_index])
    i_adder = tree.cumulative_count[i_leaf]
    for (i_element_source,element_source) in enumerate(elements_tuple)
        for i_source in 1:leaf.n_bodies[i_element_source]
            for (i_element_target,element_target) in enumerate(elements_tuple)
                for i_target in 1:leaf.n_bodies[i_element_target]
                    element.direct!(view(element_target.potential,:,leaf.first_body[i_element_target]+i_target-1), 
                    view(element.bodies,1:3,leaf.first_body[i_element_target]+i_target-1), 
                    view(element.bodies,:,leaf.first_body[i_element_source]+i_source-1)
                    )
                    S[i_target+i_adder, i_source+i_adder] = element.potential[1,leaf.first_body[i_element_target]+i_target-1]
                    element.potential[:,leaf.first_body[i_element_target]+i_target-1] .*= 0
                end
            end
        end
    end
end

@testset "S: P2P" begin

position_1 =  positions[:,1]
position_2 =  positions[:,2]
dx_vec = position_2 - position_1
dx = sqrt(dx_vec'*dx_vec)
check_potential = 1/dx

position_1 =  positions[:,3]
position_2 =  positions[:,4]
dx_vec = position_2 - position_1
dx = sqrt(dx_vec'*dx_vec)
check_potential_2 = 1/dx

@test isapprox(S[1,2], check_potential)
@test isapprox(S[2,1], check_potential)
@test isapprox(S[3,4], check_potential_2)
@test isapprox(S[4,3], check_potential_2)

end


n_max = options.expansion_order
n_coefficients = (n_max * (n_max + 1)) >> 1 + n_max + 1
harmonics = Vector{Complex{Float64}}(undef, (tree.expansion_order+1)^2)


V = spzeros(eltype(tree.branches[1].multipole_expansion[1]), n_coefficients * n_bodies, n_bodies)
for (i_leaf,leaf) in enumerate(tree.branches[tree.leaf_index])
    i_adder = tree.cumulative_count[i_leaf]
    for (i_element,element) in enumerate(elements_tuple)
        for i in 1:leaf.n_bodies[i_element]
            element.B2M!(tree, leaf, view(element.bodies,:,leaf.first_body[i_element]+i-1), 1, harmonics)
            V[(i+i_adder-1)*n_coefficients+1:(i+i_adder)*n_coefficients, i+i_adder] .= leaf.multipole_expansion[1][:]
            for j in 1:4; leaf.multipole_expansion[j] .*= 0; end # reset multipole expansion
        end
    end
end

fmm.fmm!(tree, elements_tuple, options)

@testset "V: P2M" begin

leaf_1 = tree.branches[tree.leaf_index[1]]
check_coefficients = leaf_1.multipole_expansion[1]
computed_coefficients = V[1:n_coefficients,1] + V[n_coefficients+1:2*n_coefficients,2]
for i_coeff in 1:length(check_coefficients)
    @test isapprox(check_coefficients[i_coeff], computed_coefficients[i_coeff]; atol=1e-10)
end

leaf_2 = tree.branches[tree.leaf_index[2]]
check_coefficients = leaf_2.multipole_expansion[1]
computed_coefficients = V[2*n_coefficients+1:3*n_coefficients,3] + V[3*n_coefficients+1:4*n_coefficients,4]
for i_coeff in 1:length(check_coefficients)
    @test isapprox(check_coefficients[i_coeff], computed_coefficients[i_coeff]; atol=1e-10)
end

end