using FastMultipole
using FastMultipole.StaticArrays
using LinearAlgebra
using PythonPlot
using BSON

include("../test/gravitational.jl")
include("../test/evaluate_multipole.jl")

function spherical_to_cartesian(ρ,θ,ϕ)
    z = ρ * cos(θ)
    x = ρ * sin(θ) * cos(ϕ)
    y = ρ * sin(θ) * sin(ϕ)
    return x, y, z
end

function build_system(center, ρs, θs, ϕs, qs)
    source_bodies = zeros(8,length(ρs))
    for (i,(ρ,θ,ϕ,q)) in enumerate(zip(ρs, θs, ϕs, qs))
        x, y, z = spherical_to_cartesian(ρ,θ,ϕ)
        source_bodies[1:3,i] .= center + SVector{3}(x,y,z)
        source_bodies[5,i] = q
    end
    system = Gravitational(source_bodies)
    return system
end

function generate_multipole(center, ρs, θs, ϕs, qs, expansion_order)
    # build system
    system = build_system(center, ρs, θs, ϕs, qs)

    # build branch
    branch = Branch(1:1, 0, 1:0, 0, 1, center, 0.0, expansion_order)

    # multipole coefficients
    body_to_multipole!(branch, system, branch.harmonics, Val(expansion_order))

    return branch, system
end

function evaluate_multipole(xt, branch::Branch, expansion_order)
    ϕ_m2b, v_m2b, g_m2b = evaluate_multipole(xt, branch.center, branch.multipole_expansion, DerivativesSwitch(), Val(expansion_order))
    return ϕ_m2b
end

function evaluate_direct(xt, system::Gravitational)
    u = 0.0
    for body in system.bodies
        xs = body.position
        q = body.strength
        u += q / norm(xt-xs)
    end
    return u / 4 / pi
end

function generate_local(center, multipole_branch::Branch, expansion_order)
    local_branch = Branch(1:1, 0, 1:0, 0, 1, center, 0.0, expansion_order)

    # preallocate containers
    lamb_helmholtz = Val(false)
    Hs_π2 = [1.0]
    FastMultipole.update_Hs_π2!(Hs_π2, Val(expansion_order))
    Ts = zeros(FastMultipole.length_Ts(expansion_order))
    eimϕs = zeros(2, expansion_order+1)
    weights_tmp_1 = initialize_expansion(expansion_order, eltype(Ts))
    weights_tmp_2 = initialize_expansion(expansion_order, eltype(Ts))

    # normalization
    ζs_mag = zeros(FastMultipole.length_ζs(expansion_order))
    FastMultipole.update_ζs_mag!(ζs_mag, 0, expansion_order)
    ηs_mag = zeros(FastMultipole.length_ηs(expansion_order))
    FastMultipole.update_ηs_mag!(ηs_mag, 0, expansion_order)

    # local coefficients
    FastMultipole.multipole_to_local!(local_branch, multipole_branch, weights_tmp_1, weights_tmp_2, Ts, eimϕs, ζs_mag, ηs_mag, Hs_π2, Val(expansion_order), lamb_helmholtz)

    return local_branch
end

function evaluate_local(xt, branch::Branch, expansion_order)
    Δx = xt - branch.center
    velocity_n_m = zeros(2,3,size(branch.multipole_expansion,3))
    lamb_helmholtz = Val(false)
    u, velocity, gradient = FastMultipole.evaluate_local(Δx, branch.harmonics, velocity_n_m, branch.local_expansion, Val(expansion_order), lamb_helmholtz, DerivativesSwitch())
    return u
end

function upper_bound_error(center_separation, am, al, expansion_order)
    c̃ = (center_separation - am) / al
    ε = 1/(c̃ - 1) * (1/c̃)^expansion_order
    return ε
end

function old_upper_bound_error(center_separation, am, al, expansion_order)
    c = center_separation / am - 1.0
    ε = 1/(c - 1) * (1/c)^expansion_order
    return ε*2
end

function porter_error(center_separation, am, al, expansion_order)
    d = center_separation
    multi_err = (d - am) / (d-al-am) * (am / (d-al))^(expansion_order + 1)
    local_err = (d - am) / (d-am-al) * (al / (d-am))^(expansion_order + 1)
    return multi_err, local_err
end

#--- RUN TESTS ---#



function test_error(distance::Number, epsilon)
    # define multipole (source)
    multipole_center = SVector{3}(0.0,0,0)
    ρs = [1.0]
    θs = [π/2]
    ϕs = [0.0]
    qs = [1.0]
    multipole, system = generate_multipole(multipole_center, ρs, θs, ϕs, qs, expansion_order)

    # define target
    xt = SVector{3}(1.0+distance,0.0,0.0)

    # define local expansion
    local_center = SVector{3}(7.0,0,0)
    local_branch = generate_local(local_center, multipole, expansion_order)

    # potential
    u_multipole = evaluate_multipole(xt, multipole, expansion_order)
    u_local = evaluate_local(xt, local_branch, expansion_order)
    u_direct = evaluate_direct(xt, system)

    # upper bound (expected)
    ub = upper_bound_error(norm(multipole.center - local_center), norm(system.bodies[1].position-multipole.center), norm(xt-local_branch.center), expansion_order)
    old_ub = old_upper_bound_error(norm(multipole.center - local_center), norm(system.bodies[1].position-multipole.center), norm(xt-local_branch.center), expansion_order)
    me, le = porter_error(norm(multipole.center - local_center), norm(system.bodies[1].position-multipole.center), norm(xt-local_branch.center), expansion_order)

    return abs((u_multipole - u_direct)/u_direct), abs((u_local - u_direct)/u_direct), ub, old_ub, me, le
end

function test_error(ds::Vector, expansion_order = 130)
    multipole_error = zeros(length(ds))
    local_error = zeros(length(ds))
    upper_bounds = zeros(length(ds))
    old_upper_bounds = zeros(length(ds))
    multi_upper_bounds = zeros(length(ds))
    local_upper_bounds = zeros(length(ds))
    for (i,d) in enumerate(ds)
        m, l, u, u_old, me, le = test_error(d, expansion_order)
        multipole_error[i] = m
        local_error[i] = l
        upper_bounds[i] = u
        old_upper_bounds[i] = u_old
        multi_upper_bounds[i] = me
        local_upper_bounds[i] = le
    end
    return multipole_error, local_error, upper_bounds, old_upper_bounds, multi_upper_bounds, local_upper_bounds
end

#=
ds = collect(range(-0.5,7,length=200))
p = 130
m, l, u = test_error(ds, p)

fig = figure("test error")
fig.clear()
fig.add_subplot(111, xlabel="distance", ylabel="relative error")
ax = fig.get_axes()[0]
ax.plot(ds, m)
ax.plot(ds, l)
#ax.scatter([1.0],[0.0], "*")
ax.set_yscale("log")
#ax.set_ylim([1e-16,100])
ax.legend(["multipole", "local"])

BSON.@save "error_convergence_p$p.bson" ds m l
savefig("error_convergence_p$p.png")
=#

ds = collect(range(-0.5,6,length=200))
p = 20
m, l, u, old_u, me, le = test_error(ds, p)

fig = figure("test upper bound")
fig.clear()
fig.add_subplot(111, xlabel="distance", ylabel="relative error")
ax = fig.get_axes()[0]
#ax.plot(ds, m)
ax.plot(ds, l)
#ax.plot(ds, l .- m)
#ax.plot(ds, me)
#ax.plot(ds, le)
#max_error = zeros(length(me))
#for (i,(m,l)) in enumerate(zip(me,le))
#    max_error[i] = max(m,l)
#end
#ax.plot(ds, max_error)
ax.plot(ds, old_u)
ax.plot(ds, me .+ le)
#ax.plot(ds, u, "--")
#ax.scatter([1.0],[0.0], "*")
ax.set_yscale("log")
#ax.set_ylim([1e-16,100])
#ax.legend(["multipole", "local", "upper bound", "old upper bound", "multipole only", "local only", "combined"])
ax.legend(["actual error", "literature", "novel formula"])

BSON.@save "new_upper_bound_p$p.bson" ds m l
savefig("new_upper_bound_p$p.png")
