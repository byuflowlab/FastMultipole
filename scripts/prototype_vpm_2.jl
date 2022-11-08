include("../test/vortex.jl")

# 2-D vortex ring
xs = [
    0 0
    0.5 -0.5
    0 0
]

Gammas = [
    0 0;
    0 0;
    0 0;
    1 -1.0
]

vortons = Vector{Vorton}(undef,size(xs)[2])
for i in 1:size(xs)[2]
    vortons[i] = Vorton(xs[:,i], Gammas[:,i])
end

# using direct method
fmm.direct!(vortons, P2P!)
fmm.update_velocity!(vortons)

@show isapprox(vortons[1].velocity[1], 1/4/pi; atol=1e-10)
@show isapprox(vortons[2].velocity[1], 1/4/pi; atol=1e-10)
@show isapprox(vortons[1].velocity[2], 0; atol=1e-10)
@show isapprox(vortons[2].velocity[2], 0; atol=1e-10)
@show isapprox(vortons[1].velocity[3], 0; atol=1e-10)
@show isapprox(vortons[2].velocity[3], 0; atol=1e-10)

vorton_potential_check = zeros(4,2)
vorton_potential_check[:,1] = deepcopy(vortons[1].potential)
vorton_potential_check[:,2] = deepcopy(vortons[2].potential)

vorton_velocity_check = zeros(3,2)
vorton_velocity_check[:,1] = deepcopy(vortons[1].velocity)
vorton_velocity_check[:,2] = deepcopy(vortons[2].velocity)

# reset vortons
for vorton in vortons
    vorton.velocity .*= 0
    vorton.potential .*= 0
    vorton.J_potential .*= 0
end

# manually build tree for testing
# Branch(n_branches, n_bodies, i_child, i_start, center, radius, multipole_expansion, local_expansion)
expansion_order = 9
x_branch_1 = [0.0,0,0]
branch_1 = fmm.Branch(2, 2, 2, 1, x_branch_1, 1/8, fmm.initialize_expansion(expansion_order), fmm.initialize_expansion(expansion_order))
x_branch_2 = xs[:,1] .+ [0.01, 0.02, -0.03]
branch_2 = fmm.Branch(-1, 1, -1, 1, x_branch_2, 1/8, fmm.initialize_expansion(expansion_order), fmm.initialize_expansion(expansion_order))
x_branch_3 = xs[:,2] .+ [0.02, -0.04, 0.01]
branch_3 = fmm.Branch(-1, 1, -1, 2, x_branch_3, 1/8, fmm.initialize_expansion(expansion_order), fmm.initialize_expansion(expansion_order))

# using FMM
# tree = fmm.Tree(indices, branches, [expansion_order], n_per_branch, B2M!, P2P!)
tree = fmm.Tree([1,2], [branch_1, branch_2, branch_3], [expansion_order], 1, B2M!, P2P!)

# perform FMM manually since we artificially constructed the tree
fmm.B2M!(tree, vortons, 2)
fmm.B2M!(tree, vortons, 3)
fmm.M2L!(tree, 2, 3)
fmm.M2L!(tree, 3, 2)
fmm.L2B!(tree, vortons, 2)
fmm.L2B!(tree, vortons, 3)

# theta = 2 # cutoff radius squared
# fmm.fmm!(tree, vortons, theta)

fmm.update_velocity!(vortons)

for i in 1:2
    for ind in 1:4
        @show isapprox(vortons[i].potential[ind], vorton_potential_check[ind,i]; atol=1e-12)
    end
    for dim in 1:3
        @show isapprox(vortons[i].velocity[dim], vorton_velocity_check[dim,i]; atol=1e-12)
    end
end