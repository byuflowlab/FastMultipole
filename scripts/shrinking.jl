include("../test/gravitational.jl")

function run_gravitational_fmm(c;n_bodies=1000)

    sz = (8,n_bodies)
    bodies = zeros(eltype(c),sz)
    bodies[1,:] .= 1.0
    bodies[3,:] .= 1.0
    bodies[4,:] .= 1.0
    # bodies[8,:] .= 1.0
    for i=1:n_bodies
        bodies[2,i] = c[1]*i/n_bodies
        bodies[5,i] = i/n_bodies*1.0
        bodies[6,i] = i/n_bodies*1.0
    end
    systems = (Gravitational(bodies),)
    expansion_order, n_per_branch, theta = 10, 10, 0.4
    tree = fmm.Tree(systems;expansion_order=expansion_order, n_per_branch=n_per_branch, shrink_recenter=true)
    fmm.fmm!(tree,systems;theta=theta, reset_tree=true, nearfield=true, farfield=true, unsort_bodies=true)
    @show sum(systems[1].potential)
    return [sum(systems[1].potential)], systems, tree
end

_, sys, tree = run_gravitational_fmm([1.0])

fmm.visualize("test_shrinking_nan", sys, tree; toggle_bodies=true)