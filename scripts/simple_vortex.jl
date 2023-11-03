include("../test/vortex.jl")
using Random
using WriteVTK

function generate_vortices(seed, n_bodies; radius_factor=0.1)
    Random.seed!(123)
    bodies = rand(7,n_bodies)
    bodies[4,:] ./= (n_bodies^(1/3)*2)
    bodies[4,:] .*= radius_factor
    system = VortexParticles(bodies)
    return system
end

function bm_fmm()
    expansion_order, n_per_branch, theta = 5, 5, 0.4
    n_bodies = 50
    system = generate_vortices(123, n_bodies)
    fmm.fmm!(system; expansion_order=expansion_order, n_per_branch=n_per_branch, theta=theta, nearfield=true, farfield=true, unsort_bodies=true, shrink_recenter=false)
    return system.potential[2:4,:], system.velocity_stretching[1:3,:]
end

function bm_fmm(expansion_order, n_per_branch, theta, n_bodies, shrink_recenter)
    system = generate_vortices(123, n_bodies)
    tree = fmm.fmm!(system; expansion_order=expansion_order, n_per_branch=n_per_branch, theta=theta, nearfield=true, farfield=true, unsort_bodies=true, shrink_recenter=shrink_recenter)
    return system.potential[2:4,:], system.velocity_stretching[1:3,:], tree
end

function bm_direct()
    n_bodies = 50
    system = generate_vortices(123, n_bodies)
    fmm.direct!(system, 1:n_bodies, system, 1:n_bodies)
    return system.potential[2:4,:], system.velocity_stretching[1:3,:]
end

function bm_direct(n_bodies)
    system = generate_vortices(123, n_bodies)
    fmm.direct!(system, 1:n_bodies, system, 1:n_bodies)
    return system.potential[2:4,:], system.velocity_stretching[1:3,:]
end

function bm_fmm_accuracy(expansion_order, n_per_branch, theta, n_bodies, shrink_recenter)
    system = generate_vortices(123, n_bodies)
    tree = fmm.fmm!(system; expansion_order=expansion_order, n_per_branch=n_per_branch, shrink_recenter=shrink_recenter, theta=theta, nearfield=true, farfield=true, unsort_bodies=true)
    system2 = generate_vortices(123, n_bodies)
    fmm.direct!(system2, 1:n_bodies, system2, 1:n_bodies)
    phi = system.potential[2:4,:]
    phi2 = system2.potential[2:4,:]
    v = system.velocity_stretching[1:3,:]
    v2 = system2.velocity_stretching[1:3,:]
    return maximum(abs.(phi2 - phi))/maximum(abs.(phi2)), maximum(abs.(v2 - v))/maximum(abs.(v2)), system, tree, system2
end

expansion_order, n_per_branch, theta, n_bodies, shrink_recenter = 20, 2, 0.4, 14, false
# err_potential, err_velocity, sys, tree, sys2 = bm_fmm_accuracy(expansion_order, n_per_branch, theta, n_bodies, shrink_recenter)
# println("relative potential error: $err_potential")
# println("relative velocity error: $err_velocity")

potential, v, tree = bm_fmm(expansion_order, n_per_branch, theta, n_bodies, shrink_recenter)
potential_direct, v_direct = bm_direct(n_bodies)

err_potential = zeros(size(potential,2))
for i in 1:length(err_potential)
    err_potential[i] = max(abs.(potential[:,i] - potential_direct[:,i])...)
end