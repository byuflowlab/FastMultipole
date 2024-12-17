using Statistics, Random, LegendrePolynomials
using PythonPlot
using LaTeXStrings
using LinearAlgebra
include("../test/gravitational.jl")
include("../test/bodytolocal.jl")
include("../test/evaluate_multipole.jl")

function output_stuff(v; verbose=false)
    if verbose
        @show mean(abs.(v)) maximum(abs.(v)) minimum(abs.(v)) std(abs.(v))
    end
    return mean(abs.(v)), maximum(abs.(v)), minimum(abs.(v)), std(abs.(v))
end

function spherical_to_cartesian(ρ,θ,ϕ)
    z = ρ * cos(θ)
    x = ρ * sin(θ) * cos(ϕ)
    y = ρ * sin(θ) * sin(ϕ)
    return x, y, z
end

function build_system(center, ρs, θs, ϕs, radii, qs)
    source_bodies = zeros(8,length(ρs))
    for (i,(ρ,θ,ϕ,r,q)) in enumerate(zip(ρs, θs, ϕs, radii, qs))
        x, y, z = spherical_to_cartesian(ρ,θ,ϕ)
        source_bodies[1:3,i] .= center + SVector{3}(x,y,z)
        source_bodies[4,i] = r
        source_bodies[5,i] = q
    end
    system = Gravitational(source_bodies)
    return system
end

function build_system_range(center, ρ_range, θ_range, ϕ_range, r_range, A)
    ρs = Float64[]
    θs = Float64[]
    ϕs = Float64[]
    rs = Float64[]
    qs = Float64[]
    q_density = A / (4/3*pi*ρ_range[end]^3)
    if length(ρ_range) > 1
        for iρ in 1:length(ρ_range)-1
            ρ = (ρ_range[iρ] + ρ_range[iρ+1])/2
            dV = 4/3*pi*(ρ_range[iρ+1]^3 - ρ_range[iρ]^3)
            total_q = dV * q_density
            q = total_q / (length(θ_range) * length(ϕ_range))
            total_circ = 0.0
            for θ in θ_range
                total_circ += 2*pi*ρ*sin(θ)
            end
            dq = q / total_circ
            for θ in θ_range
                circ = 2*pi*ρ*sin(θ)
                this_q = circ * dq / length(ϕ_range)
                for ϕ in ϕ_range
                    r = r_range[1] + rand() * (r_range[2]-r_range[1])
                    #q = q_range[1] + rand() * (q_range[2]-q_range[1])
                    push!(ρs, ρ)
                    push!(θs, θ)
                    push!(ϕs, ϕ)
                    push!(rs, r)
                    push!(qs, this_q)
                end

            end
        end
    else
        ρ = ρ_range[end]
        q = A / (length(θ_range) * length(ϕ_range))
        for θ in θ_range
            for ϕ in ϕ_range
                r = r_range[1] + rand() * (r_range[2]-r_range[1])
                #q = q_range[1] + rand() * (q_range[2]-q_range[1])
                push!(ρs, ρ)
                push!(θs, θ)
                push!(ϕs, ϕ)
                push!(rs, r)
                push!(qs, q)
            end
        end
    end
    return build_system(center, ρs, θs, ϕs, rs, qs)
end

function build_system_cartesian(center, xs, ys, zs, radii, qs)
    source_bodies = zeros(8,length(xs))
    for (i,(x,y,z,r,q)) in enumerate(zip(xs, ys, zs, radii, qs))
        source_bodies[1:3,i] .= center + SVector{3}(x,y,z)
        source_bodies[4,i] = r
        source_bodies[5,i] = q
    end
    system = Gravitational(source_bodies)
    return system
end

function build_system_cartesian_range(center, x_range, y_range, z_range, r_range, q_range)
    xs = Float64[]
    ys = Float64[]
    zs = Float64[]
    rs = Float64[]
    qs = Float64[]
    for x in x_range
        for y in y_range
            for z in z_range
                push!(xs, x)
                push!(ys, y)
                push!(zs, z)
                r = r_range[1] + rand() * (r_range[2]-r_range[1])
                q = q_range[1] + rand() * (q_range[2]-q_range[1])
                push!(rs, r)
                push!(qs, q)
            end
        end
    end
    return build_system_cartesian(center, xs, ys, zs, rs, qs)
end

function generate_branch(center, ρ_range, θ_range, ϕ_range, r_range, A, expansion_order, shrink)
    # build system
    system = build_system_range(center, ρ_range, θ_range, ϕ_range, r_range, A)

    # build branch
    branch = Branch(1:get_n_bodies(system), 0, 1:0, 0, 1, center, sqrt(3.0), sqrt(3.0), SVector{6}(-1.0,1.0,-1.0,1.0,-1.0,1.0), SVector{3}(1.0,1.0,1.0), expansion_order)

    # shrink branch
    branches = [branch]
    shrink && FastMultipole.shrink_recenter!(branches, [1:1], system)
    branch = branches[1]

    return branch, system
end

function generate_branch_cartesian(center, x_range, y_range, z_range, r_range, q_range, expansion_order, shrink)
    # build system
    system = build_system_cartesian_range(center, x_range, y_range, z_range, r_range, q_range)

    # build branch
    branch = Branch(1:get_n_bodies(system), 0, 1:0, 0, 1, center, sqrt(3.0), sqrt(3.0), SVector{6}(-1.0,1.0,-1.0,1.0,-1.0,1.0), SVector{3}(1.0,1.0,1.0), expansion_order)

    # shrink branch
    branches = [branch]
    shrink && FastMultipole.shrink_recenter!(branches, [1:1], system)
    branch = branches[1]

    return branch, system
end

function evaluate_multipole(xt, branch::Branch, expansion_order, lamb_helmholtz::Bool)
    Δx = xt - branch.source_center
    ϕ_m2b, v_m2b, g_m2b = evaluate_multipole(Δx, branch.harmonics, branch.multipole_expansion, Val(expansion_order), Val(lamb_helmholtz), DerivativesSwitch(true,true,false))
    return ϕ_m2b, v_m2b
end

function evaluate_direct(xt, system::Gravitational, multipole_branch::SingleBranch)
    u = 0.0
    for i in multipole_branch.bodies_index
        body = system.bodies[i]
        xs = body.position
        q = body.strength
        u += q / norm(xt-xs)
    end
    return u / 4 / pi
end

function multipole_to_local!(local_branch, multipole_branch, expansion_order, lamb_helmholtz)
    # reset expansion
    local_branch.local_expansion .= zero(eltype(local_branch.local_expansion))

    # preallocate containers
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
    FastMultipole.multipole_to_local!(local_branch, multipole_branch, weights_tmp_1, weights_tmp_2, Ts, eimϕs, ζs_mag, ηs_mag, Hs_π2, expansion_order, lamb_helmholtz)

    return nothing
end

function evaluate_local(xt, branch::Branch, expansion_order, lamb_helmholtz)
    Δx = xt - branch.target_center
    velocity_n_m = zeros(2,3,size(branch.multipole_expansion,3))
    u, velocity, gradient = FastMultipole.evaluate_local(Δx, branch.harmonics, velocity_n_m, branch.local_expansion, Val(expansion_order), Val(lamb_helmholtz), DerivativesSwitch())
    return u, velocity
end

function relative_error(u_true, u_test)
    return (u_test .- u_true) ./ u_true
end

function absolute_error(u_true::Vector{Float64}, u_test::Vector{Float64})
    return u_test .- u_true
end

function reset!(system::Gravitational)
    system.potential .= 0.0
end

function reset!(system::VortexParticles)
    system.potential .= 0.0
    system.velocity_stretching .= 0.0
end

function absolute_error(direct, mp)
    delta = similar(direct)
    for i in eachindex(direct)
        delta[i] = direct[i] - mp[i]
    end
    err = [norm(v) for v in delta]
    return err
end

function multipole_local_error(local_branch::SingleBranch, local_system, multipole_branch::SingleBranch, multipole_system, expansion_order::Int; body_type=Point{Source})
    # reset bodies
    reset!(local_system)

    # get direct potential
    # potential_direct = [evaluate_direct(local_system[i,Position()], multipole_system, multipole_branch) for i in local_branch.bodies_index]
    FastMultipole._direct!(local_system, local_branch.bodies_index, DerivativesSwitch(true,true,false), multipole_system, multipole_branch.bodies_index)
    potential_direct = [local_system[i,FastMultipole.ScalarPotential()] for i in local_branch.bodies_index]
    velocity_direct = [SVector{3}(local_system[i,FastMultipole.Velocity()]) for i in local_branch.bodies_index]

    # multipole expansion
    multipole_branch.multipole_expansion .= 0.0
    FastMultipole.body_to_multipole!(multipole_branch, multipole_system, multipole_branch.harmonics, Val(expansion_order))

    # multipole error
    lamb_helmholtz = body_type == Point{Vortex}
    potential_mp = [evaluate_multipole(local_system[i,Position()], multipole_branch, expansion_order, lamb_helmholtz)[1] for i in local_branch.bodies_index]
    velocity_mp = [evaluate_multipole(local_system[i,Position()], multipole_branch, expansion_order, lamb_helmholtz)[2] for i in local_branch.bodies_index]
    error_mp = absolute_error(potential_direct, potential_mp)
    error_mp_velocity = absolute_error(velocity_direct, velocity_mp)

    # local expansion
    if body_type == Point{Source}
        local_branch.local_expansion .= 0.0
        x_l = local_branch.target_center
        for i in multipole_branch.bodies_index
            Δx = multipole_system[i,Position()] - x_l
            body_to_local_point!(body_type, local_branch.local_expansion, local_branch.harmonics, Δx, multipole_system[i, Strength()], Val(expansion_order))
        end

        # local error
        lamb_helmholtz = false
        potential_l = [evaluate_local(local_system[i,Position()], local_branch, expansion_order, lamb_helmholtz)[1] for i in local_branch.bodies_index]
        velocity_l = [evaluate_local(local_system[i,Position()], local_branch, expansion_order, lamb_helmholtz)[2] for i in local_branch.bodies_index]
        error_l = absolute_error(potential_direct, potential_l)
        error_l_velocity = absolute_error(velocity_direct, velocity_l)
    end

    # overall error
    local_branch.local_expansion .= 0.0
    multipole_to_local!(local_branch, multipole_branch, expansion_order, Val(lamb_helmholtz))
    potential_o = [evaluate_local(local_system[i,Position()], local_branch, expansion_order, lamb_helmholtz)[1] for i in local_branch.bodies_index]
    velocity_o = [evaluate_local(local_system[i,Position()], local_branch, expansion_order, lamb_helmholtz)[2] for i in local_branch.bodies_index]
    error_o = absolute_error(potential_direct, potential_o)
    error_o_velocity = absolute_error(velocity_direct, velocity_o)


    if body_type !== Point{Source}
        error_l = abs.(error_o) .- abs.(error_mp)
        dV_mp = [velocity_direct[i] - velocity_mp[i] for i in eachindex(velocity_mp)]
        dV_o = [velocity_direct[i] - velocity_o[i] for i in eachindex(velocity_o)]
        dV = [dV_o[i] - dV_mp[i] for i in eachindex(dV_mp)]
        error_l_velocity = [norm(V) for V in dV]
        # error_l_velocity = abs.(error_o_velocity) .- abs.(error_mp_velocity)
    end

    return error_mp, error_l, error_o, error_mp_velocity, error_l_velocity, error_o_velocity
end

function mymax(vect)
    mp, mn = 0.0, 0.0
    for v in vect
        mp = max(v,mp)
        mn = min(v,mn)
    end
    if mp > -mn
        return mp
    else
        return mn
    end
end

function test_error(local_branch, local_system, multipole_branch, multipole_system, expansion_order::Int)
    # experimental error
    if typeof(multipole_system) <: Gravitational
        body_type = Point{Source}
    elseif typeof(multipole_system) <: VortexParticles
        body_type = Point{Vortex}
    end
    e_mp, e_l, e_o, ev_mp, ev_l, ev_o = multipole_local_error(local_branch, local_system, multipole_branch, multipole_system, expansion_order; body_type)

    # unequal spheres
    ub_mp_us = multipole_error(local_branch, multipole_branch, expansion_order, UnequalSpheres())
    ub_l_us = local_error(local_branch, multipole_branch, expansion_order, UnequalSpheres())

    # uniform unequal spheres
    ub_mp_uus = multipole_error(local_branch, multipole_branch, expansion_order, UniformUnequalSpheres())
    ub_l_uus = local_error(local_branch, multipole_branch, expansion_order, UniformUnequalSpheres())

    # unequal boxes
    ub_mp_ub = multipole_error(local_branch, multipole_branch, expansion_order, UnequalBoxes())
    ub_l_ub = local_error(local_branch, multipole_branch, expansion_order, UnequalBoxes())

    # uniform unequal boxes
    ub_mp_uub = multipole_error(local_branch, multipole_branch, expansion_order, UniformUnequalBoxes())
    ub_l_uub = local_error(local_branch, multipole_branch, expansion_order, UniformUnequalBoxes())

    # uniform cubes
    e_mp_uc = multipole_error(local_branch, multipole_branch, expansion_order, UniformCubes())
    e_l_uc = local_error(local_branch, multipole_branch, expansion_order, UniformCubes())

    # uniform cubes velocity
    ev_mp_uc = multipole_error(local_branch, multipole_branch, expansion_order, UniformCubesVelocity())
    ev_l_uc = local_error(local_branch, multipole_branch, expansion_order, UniformCubesVelocity())

    return mymax(e_mp), mymax(e_l), mymax(e_o), mymax(ev_mp), mymax(ev_l), mymax(ev_o), ub_mp_us, ub_l_us, ub_mp_uus, ub_l_uus, ub_mp_ub, ub_l_ub, ub_mp_uub, ub_l_uub, e_mp_uc, e_l_uc, ev_mp_uc, ev_l_uc
end

function test_rms_error(R, r_multipole, r_local, nx_source_bodies, nx_target_bodies, expansion_order)
    source_bodies = zeros(7,nx_source_bodies^3)
    i = 1
    for x in range(-r_multipole, r_multipole, nx_source_bodies)
        for y in range(-r_multipole, r_multipole, nx_source_bodies)
            for z in range(-r_multipole, r_multipole, nx_source_bodies)
                source_bodies[1:3,i] .= x, y, z
                i += 1
            end
        end
    end
    source_bodies[5,:] .= 1.0

    target_bodies = zeros(7,nx_target_bodies^3)
    i = 1
    for x in range(R[1]-r_local, R[1]+r_local, nx_target_bodies)
        for y in range(R[2]-r_local, R[2]+r_local, nx_target_bodies)
            for z in range(R[3]-r_local, R[3]+r_local, nx_target_bodies)
                target_bodies[1:3,i] .= x, y, z
                i += 1
            end
        end
    end

    bodies = hcat(source_bodies, target_bodies)
    system = Gravitational(bodies)
    root = Branch(1:size(bodies,2), 0, 2:3, 0, 1, R/2, 0.5*norm(R)+max(r_multipole,r_local), 0.5*norm(R)+max(r_multipole,r_local), SVector{6}(0.0,0,0,0,0,0), SVector{3}(0.0,0,0), expansion_order)
    source_branch = Branch(1:nx_source_bodies^3, 0, 1:0, 0, 1, SVector{3}(0.0,0,0), r_multipole*sqrt(3.0), r_multipole*sqrt(3.0), SVector{6}(-1.0,1.0,-1.0,1.0,-1.0,1.0)*r_multipole, SVector{3}(1.0,1.0,1.0)*r_multipole, expansion_order)
    target_branch = Branch(nx_source_bodies^3+1:nx_source_bodies^3+nx_target_bodies^3, 0, 1:0, 0, 1, R, r_local*sqrt(3.0), r_local*sqrt(3.0), SVector{6}(-1.0,1.0,-1.0,1.0,-1.0,1.0)*r_local, SVector{3}(1.0,1.0,1.0)*r_local, expansion_order)

    # experimental error
    e_mp, e_l, e_o = multipole_local_error(target_branch, system, source_branch, system, expansion_order)
    FastMultipole.evaluate_local!(system, target_branch, Val(expansion_order), Val(false), DerivativesSwitch())

    #--- visualize error with paraview ---#

    system_copy = deepcopy(system)
    system_copy.potential .= 0.0
    for i_source in source_branch.bodies_index
        x_source = system_copy.bodies[i_source].position
        q = system_copy.bodies[i_source].strength
        for i_target in target_branch.bodies_index
            x_target = system_copy.bodies[i_target].position
            system_copy.potential[1,i_target] += q/norm(x_source-x_target)/(4*pi)
        end
    end

    # local error
    local_potential = deepcopy(system.potential[1,:])
    local_error = local_potential .- system_copy.potential[1,:]
    system.potential[1,:] .-= system_copy.potential[1,:]
    for i in 1:get_n_bodies(system)
        if iszero(system.potential[1,i])
            system.potential[1,i] = 1e-16
        end
    end
    tree = Tree([root, source_branch, target_branch], [1:1, 2:3], [3], collect(1:size(bodies,2)), collect(1:size(bodies,2)), system_copy, Val(expansion_order), size(target_bodies,2))
    FastMultipole.visualize_bodies("local_error_check", system, tree)


    # multipole error
    system.potential .= 0.0
    evaluate_multipole!(system, 1+nx_source_bodies^3:get_n_bodies(system), source_branch.multipole_expansion, DerivativesSwitch(), Val(expansion_order), source_branch.source_center)
    multipole_potential = deepcopy(system.potential[1,:])
    multipole_error = multipole_potential .- system_copy.potential[1,:]
    system.potential[1,:] .-= system_copy.potential[1,:]
    for i in 1:get_n_bodies(system)
        if iszero(system.potential[1,i])
            system.potential[1,i] = 1e-16
        end
    end
    FastMultipole.visualize_bodies("multipole_error_check", system, tree)

    return system, system_copy, multipole_potential, local_potential, e_mp, e_l, e_o
end

function closest_corner(target_branch, source_branch)
    x0 = source_branch.source_center
    dx = Inf

    # check midpoints
    lx, ly, lz = target_branch.target_box
    points = [
                 target_branch.target_center + SVector{3}(lx,ly,lz),
                 target_branch.target_center + SVector{3}(-lx,ly,lz),
                 target_branch.target_center + SVector{3}(lx,-ly,lz),
                 target_branch.target_center + SVector{3}(lx,ly,-lz),
                 target_branch.target_center + SVector{3}(-lx,-ly,lz),
                 target_branch.target_center + SVector{3}(-lx,ly,-lz),
                 target_branch.target_center + SVector{3}(lx,-ly,-lz),
                 target_branch.target_center + SVector{3}(-lx,-ly,-lz)
                ]
    xt = SVector{3}(0.0,0.0,0.0)
    for pt in points
        if norm(pt - x0) < dx
            dx = norm(pt-x0)
            xt = pt
        end
    end

    return xt
end

function get_points(target_branch)
    lx, ly, lz = target_branch.target_box
    points = [
                 target_branch.target_center,
                 target_branch.target_center + SVector{3}(lx,0.0,0.0),
                 target_branch.target_center + SVector{3}(0.0,ly,0.0),
                 target_branch.target_center + SVector{3}(0.0,0.0,lz),
                 target_branch.target_center - SVector{3}(lx,0.0,0.0),
                 target_branch.target_center - SVector{3}(0.0,ly,0.0),
                 target_branch.target_center - SVector{3}(0.0,0.0,lz),
                 target_branch.target_center + SVector{3}(lx,ly,lz),
                 target_branch.target_center + SVector{3}(-lx,ly,lz),
                 target_branch.target_center + SVector{3}(lx,-ly,lz),
                 target_branch.target_center + SVector{3}(lx,ly,-lz),
                 target_branch.target_center + SVector{3}(-lx,-ly,lz),
                 target_branch.target_center + SVector{3}(-lx,ly,-lz),
                 target_branch.target_center + SVector{3}(lx,-ly,-lz),
                 target_branch.target_center + SVector{3}(-lx,-ly,-lz),
                 target_branch.target_center + SVector{3}(lx,0.0,0.0)*0.5,
                 target_branch.target_center + SVector{3}(0.0,ly,0.0)*0.5,
                 target_branch.target_center + SVector{3}(0.0,0.0,lz)*0.5,
                 target_branch.target_center - SVector{3}(lx,0.0,0.0)*0.5,
                 target_branch.target_center - SVector{3}(0.0,ly,0.0)*0.5,
                 target_branch.target_center - SVector{3}(0.0,0.0,lz)*0.5,
                 target_branch.target_center + SVector{3}(lx,ly,lz)*0.5,
                 target_branch.target_center + SVector{3}(-lx,ly,lz)*0.5,
                 target_branch.target_center + SVector{3}(lx,-ly,lz)*0.5,
                 target_branch.target_center + SVector{3}(lx,ly,-lz)*0.5,
                 target_branch.target_center + SVector{3}(-lx,-ly,lz)*0.5,
                 target_branch.target_center + SVector{3}(-lx,ly,-lz)*0.5,
                 target_branch.target_center + SVector{3}(lx,-ly,-lz)*0.5,
                 target_branch.target_center + SVector{3}(-lx,-ly,-lz)*0.5
                ]
    return points
end

function closest_point(target_branch, source_branch)
    x0 = source_branch.source_center
    dx = Inf

    # check midpoints
    lx, ly, lz = target_branch.target_box
    points = [
                 target_branch.target_center + SVector{3}(lx,0.0,0.0),
                 target_branch.target_center + SVector{3}(0.0,ly,0.0),
                 target_branch.target_center + SVector{3}(0.0,0.0,lz),
                 target_branch.target_center - SVector{3}(lx,0.0,0.0),
                 target_branch.target_center - SVector{3}(0.0,ly,0.0),
                 target_branch.target_center - SVector{3}(0.0,0.0,lz),
                 target_branch.target_center + SVector{3}(lx,ly,lz),
                 target_branch.target_center + SVector{3}(-lx,ly,lz),
                 target_branch.target_center + SVector{3}(lx,-ly,lz),
                 target_branch.target_center + SVector{3}(lx,ly,-lz),
                 target_branch.target_center + SVector{3}(-lx,-ly,lz),
                 target_branch.target_center + SVector{3}(-lx,ly,-lz),
                 target_branch.target_center + SVector{3}(lx,-ly,-lz),
                 target_branch.target_center + SVector{3}(-lx,-ly,-lz)
                ]
    xt = SVector{3}(0.0,0.0,0.0)
    for pt in points
        if norm(pt - x0) < dx
            dx = norm(pt-x0)
            xt = pt
        end
    end

    return xt
end

function f_gpt(r, s, n)
    term1 = (-2 * (s * sqrt(2) / 2)^(n + 3)) / (n + 3) + (3 / 2) * s * (s * sqrt(2) / 2)^(n + 2) / (n + 2)
    term2 = (-2 * (s / 2)^(n + 3)) / (n + 3) + (3 / 2) * s * (s / 2)^(n + 2) / (n + 2)
    term3 = ((s * 3 * sqrt(6) - 4) / (2 * s * (sqrt(3) - sqrt(2)))) * (s * sqrt(3) / 2)^(n + 3) / (n + 3) - ((3 * sqrt(2) - 4) / (s * (sqrt(3) - sqrt(2)))) * (s * sqrt(3) / 2)^(n + 4) / (n + 4)
    term4 = ((s * 3 * sqrt(6) - 4) / (2 * s * (sqrt(3) - sqrt(2)))) * (s * sqrt(2) / 2)^(n + 3) / (n + 3) - ((3 * sqrt(2) - 4) / (s * (sqrt(3) - sqrt(2)))) * (s * sqrt(2) / 2)^(n + 4) / (n + 4)

    return (1 / r^(n + 1)) * (term1 - term2 + term3 - term4)
end

function multipole_preintegration(rvec, n_max, s, nx)
    dx = s/nx
    x = y = z = -0.5*s + dx*0.5
    R2 = s*s*0.25 # squared radius of largest enclosed sphere
    r = norm(rvec)
    r̂x, r̂y, r̂z = rvec / r # unit vector pointing to target
    rinv = 1/r

    res = zeros(n_max) # preallocate result
    dV = dx * dx * dx # differential volume
    for ix in 1:nx
        y = -0.5*s + dx*0.5
        for iy in 1:nx
            z = -0.5*s + dx*0.5
            for iz in 1:nx
                ρ2 = x*x + y*y + z*z
                if ρ2 > R2 # only integrate outside the sphere
                    ρ = sqrt(ρ2)
                    cθ = (r̂x * x + r̂y * y + r̂z * z) / ρ # cosine theta
                    ρn = ρ # ρ^n
                    Pnm1 = 1.0 # Legendre polynomial of degree n-1
                    Pn = cθ # Legendre polynomial of degree n
                    #ρr = ρ * rinv
                    #ρr_n = ρr

                    # update integral for each n
                    for n in 1:n_max
                        # res[n] += ρr_n * Pn * rinv # actual integral for ρ/r
                        res[n] += ρn * Pn

                        # next Legendre polynomial
                        Pnp1 = ((2n+1) * cθ * Pn - n * Pnm1) / (n+1)

                        # recurse
                        Pnm1 = Pn
                        Pn = Pnp1
                        ρn *= ρ
                        #ρr_n *= ρr
                    end
                end
                z += dx # increment z
            end
            y += dx # increment y
        end
        x += dx # increment x
    end

    # finish integration
    res .*= dV

    return res
end

function multipole_full_integration(rvec, n_max, s, nx)
    dx = s/nx
    x = y = z = -0.5*s + dx*0.5
    R2 = s*s*0.25 # squared radius of largest enclosed sphere
    r = norm(rvec)
    r̂x, r̂y, r̂z = rvec / r # unit vector pointing to target
    rinv = 1/r

    res = zeros(n_max) # preallocate result
    dV = dx * dx * dx # differential volume
    for ix in 1:nx
        y = -0.5*s + dx*0.5
        for iy in 1:nx
            z = -0.5*s + dx*0.5
            for iz in 1:nx
                ρ2 = x*x + y*y + z*z
                if ρ2 > R2 # only integrate outside the sphere
                    ρ = sqrt(ρ2)
                    cθ = (r̂x * x + r̂y * y + r̂z * z) / ρ # cosine theta
                    ρn = ρ # ρ^n
                    Pnm1 = 1.0 # Legendre polynomial of degree n-1
                    Pn = cθ # Legendre polynomial of degree n
                    ρr = ρ * rinv
                    ρr_n = ρr

                    # update integral for each n
                    for n in 1:n_max
                        res[n] += ρr_n * Pn * rinv # actual integral for ρ/r
                        # res[n] += ρn * Pn

                        # next Legendre polynomial
                        Pnp1 = ((2n+1) * cθ * Pn - n * Pnm1) / (n+1)

                        # recurse
                        Pnm1 = Pn
                        Pn = Pnp1
                        ρn *= ρ
                        ρr_n *= ρr
                    end
                end
                z += dx # increment z
            end
            y += dx # increment y
        end
        x += dx # increment x
    end

    # finish integration
    res .*= dV

    return res
end

function multipole_preintegration(n_max, s, nx, nθ, nϕ)
    res = zeros(n_max, nθ, nϕ)
    for (iθ,θ) in enumerate(range(0, stop=pi, length=nθ))
    # for (icθ,cθ) in enumerate(range(1.0, stop=-1.0, length=nθ))
        for (iϕ,ϕ) in enumerate(range(0, stop=pi, length=nϕ))
            # sθ = sqrt(1-cθ*cθ)
            sθ, cθ = sincos(θ)
            sϕ, cϕ = sincos(ϕ)
            r = SVector{3}(sθ*cϕ, sθ*sϕ, cθ)
            # res[:,icθ,iϕ] .= multipole_preintegration(r, n_max, s, nx)
            res[:,iθ,iϕ] .= multipole_preintegration(r, n_max, s, nx)
        end
    end
    return res
end

function get_iΔθ(Δθr,nθ)
    pi05 = pi*0.5
    dΔθ = pi05/nθ
    Δθ = 0.0
    for i in 1:nθ
        if Δθ >= Δθr
            Δθ-Δθr > Δθr-Δθ+dΔθ && (return i-1)
            return i
        end
        Δθ += dΔθ
    end
    pi05-Δθr > Δθr-pi05+dΔθ && (return nθ-1)
    return nθ
end

function local_preintegration(dtheta,n; R=1.0)
    dtheta < 1e-6 && (dtheta += 1e-6)
    drho = R*sin(dtheta)
	nx = 30
	dx = drho*2/nx
	x0 = -drho + dx*0.5
	y0 = -drho + dx*0.5
	z0 = R - drho + dx*0.5
	cx, cy, cz = 0.0, 0.0, R
	val = 0.0
	for x in x0:dx:x0+2*drho
		for y in y0:dx:y0+2*drho
			for z in z0:dx:z0+2*drho
				if (x^2+y^2+(z-cz)^2) <= drho*drho
                    rho = sqrt(x*x+y*y+z*z)
					val += rho^(-n-1)
				end
			end
		end
	end
	val *= dx*dx*dx
	return val
end

function local_preintegration(n_max, θmax, nθ)
    res = zeros(n_max, nθ)
    for (iθ,Δθ) in enumerate(range(0, stop=pi/2, length=nθ))
        for n in 1:n_max
            res[n,iθ] = local_preintegration(Δθ,n)
        end
    end
    return res
end

function local_error_r(x_t, l_branch, m_branch, system, x0, y0, z0, dx, dy, dz, nx, ny, nz, lx, ly, lz, A, P; method="loop", debug=false)
    integral_local = 0.0

    if method == "numerical"
        rvec = x_t - l_branch.target_center
        r = norm(x_t - l_branch.target_center)
        for z in range(z0, step=dz, length=nz)
            for y in range(y0, step=dy, length=ny)
                for x in range(x0, step=dx, length=nx)
                    rho_vec = SVector{3}(x,y,z) + m_branch.source_center - l_branch.target_center
                    rho = norm(rho_vec)
                    cos_gamma = dot(rvec, rho_vec) / (r * rho)
                    for n in P+1:P+2
                        rho_inv = 1/rho
                        Pn = Plm(cos_gamma, n, 0)
                        integral_local += (r * rho_inv)^n * rho_inv * Pn
                    end
                end
            end
        end
        integral_local *= A * dx*dy*dz / (8*lx*ly*lz)

    elseif method == "sphere"
        rvec = x_t - l_branch.target_center
        r = norm(x_t - l_branch.target_center)
        cosines = Float64[]
        s_over_2 = mean([lx,ly,lz])
        # radius2 = (x0*x0 + y0*y0 + z0*z0) / 3 # 3 makes this the largest encompassed sphere
        # radius2 = 3/(4pi)^(0.666666666) * s2

        for z in range(z0, step=dz, length=nz)
            for y in range(y0, step=dy, length=ny)
                for x in range(x0, step=dx, length=nx)
                    rhop = SVector{3}(x,y,z)
                    if rhop'*rhop < s_over_2 * s_over_2 # radius2 # contained in the largest encompassed sphere
                        rho_vec = rhop + m_branch.source_center - l_branch.target_center
                        rho = norm(rho_vec)
                        cos_gamma = dot(rvec, rho_vec) / (r * rho)
                        push!(cosines, cos_gamma)
                        for n in P+1:P+2
                            rho_inv = 1/rho
                            integral_local += (r * rho_inv)^n * rho_inv * Plm(cos_gamma, n, 0)
                        end
                    end
                end
            end
        end
        integral_local *= A * dx*dy*dz / (8*lx*ly*lz)

    elseif method == "separate_radial_sphere"
        rvec = x_t - l_branch.target_center
        r = norm(x_t - l_branch.target_center)
        cosines = Float64[]
        s_over_2 = mean([lx,ly,lz])
        # radius2 = (x0*x0 + y0*y0 + z0*z0) / 3 # 3 makes this the largest encompassed sphere
        # radius2 = 3/(4pi)^(0.666666666) * s2

        for z in range(z0, step=dz, length=nz)
            for y in range(y0, step=dy, length=ny)
                for x in range(x0, step=dx, length=nx)
                    rhop = SVector{3}(x,y,z)
                    if rhop'*rhop < s_over_2 * s_over_2 # radius2 # contained in the largest encompassed sphere
                        rho_vec = rhop + m_branch.source_center - l_branch.target_center
                        rho = norm(rho_vec)
                        cos_gamma = dot(rvec, rho_vec) / (r * rho)
                        push!(cosines, cos_gamma)
                        for n in P+1:P+1
                            rho_inv = 1/rho
                            integral_local += (rho_inv)^(n+1) * Plm(cos_gamma, n, 0)
                        end
                    end
                end
            end
        end
        integral_local *= dx*dy*dz / (8*lx*ly*lz)
        integral_local *= A * r^(P+1)

    elseif method == "separate_angular_sphere"
        rvec = x_t - l_branch.target_center
        r = norm(x_t - l_branch.target_center)
        cosines = Float64[]
        s2 = x0*x0 + y0*y0 + z0*z0
        radius2 = (x0*x0 + y0*y0 + z0*z0) / 3 # 3 makes this the largest encompassed sphere
        # radius2 = 3/(4pi)^(0.666666666) * s2
        this_cos = dot(rvec, m_branch.source_center-l_branch.target_center) / (r*norm(m_branch.source_center-l_branch.target_center))
        this_P = Plm(this_cos, P+1,0)
        for z in range(z0, step=dz, length=nz)
            for y in range(y0, step=dy, length=ny)
                for x in range(x0, step=dx, length=nx)
                    rhop = SVector{3}(x,y,z)
                    if rhop'*rhop < radius2 # contained in the largest encompassed sphere
                        rho_vec = rhop + m_branch.source_center - l_branch.target_center
                        rho = norm(rho_vec)
                        cos_gamma = dot(rvec, rho_vec) / (r * rho)
                        push!(cosines, cos_gamma)
                        for n in P+1:P+2
                            this_integral = 0.0
                            rho_inv = 1/rho
                            integral_local += (r * rho_inv)^n * rho_inv #* this_P#Plm(cos_gamma, n, 0)
                            this_integral += (r * rho_inv)^n * rho_inv #* this_P#Plm(cos_gamma, n, 0)
                        end
                    end
                end
            end
        end

        integral_local *= A * dx*dy*dz / (8*lx*ly*lz) * this_P

    elseif method == "precomputed"
        rvec = x_t - l_branch.target_center
        r = norm(x_t - l_branch.target_center)
        cosines = Float64[]
        # s2 = x0*x0 + y0*y0 + z0*z0
        # radius2 = (x0*x0 + y0*y0 + z0*z0) / 3 # 3 makes this the largest encompassed sphere
        radius2 = mean(l_branch.target_box)^2
        this_cos = dot(rvec, m_branch.source_center-l_branch.target_center) / (r*norm(m_branch.source_center-l_branch.target_center))
        this_P = Plm(this_cos, P+1,0)
        R = norm(m_branch.source_center - l_branch.target_center)
        Δθ = asin(sqrt(radius2)/R)
        iΔθ = get_iΔθ(Δθ, size(FastMultipole.LOCAL_INTEGRALS,2))
        Rinv = 1/R
        R_np2 = R
        rn = r
        for n in 1:P
            R_np2 *= Rinv
            rn *= r
        end

        for n in P+1:P+2
            pint = FastMultipole.LOCAL_INTEGRALS[n,5,5] * R_np2
            integral_local += pint * rn
            R_np2 *= Rinv
            rn *= r
        end
        integral_local *= A / (8*lx*ly*lz) * this_P

    elseif method == "separate_angular"
        rvec = x_t - l_branch.target_center
        r = norm(x_t - l_branch.target_center)
        cosines = Float64[]
        for z in range(z0, step=dz, length=nz)
            for y in range(y0, step=dy, length=ny)
                for x in range(x0, step=dx, length=nx)
                    rho_vec = SVector{3}(x,y,z) + m_branch.source_center - l_branch.target_center
                    rho = norm(rho_vec)
                    cos_gamma = dot(rvec, rho_vec) / (r * rho)
                    push!(cosines, cos_gamma)
                    for n in P+1:20
                        this_cos = dot(rvec, m_branch.source_center-l_branch.target_center) / (r*norm(m_branch.source_center-l_branch.target_center))
                        this_P = Plm(this_cos, n,0)
                        rho_inv = 1/rho
                        integral_local += (r * rho_inv)^n * rho_inv * this_P
                    end
                end
            end
        end
        integral_local *= A * dx*dy*dz / (8*lx*ly*lz)

    elseif method == "separate_angular"
        rvec = x_t - l_branch.target_center
        r = norm(x_t - l_branch.target_center)
        cosines = Float64[]
        for z in range(z0, step=dz, length=nz)
            for y in range(y0, step=dy, length=ny)
                for x in range(x0, step=dx, length=nx)
                    rho_vec = SVector{3}(x,y,z) + m_branch.source_center - l_branch.target_center
                    rho = norm(rho_vec)
                    cos_gamma = dot(rvec, rho_vec) / (r * rho)
                    push!(cosines, cos_gamma)
                    for n in P+1:20
                        this_cos = dot(rvec, m_branch.source_center-l_branch.target_center) / (r*norm(m_branch.source_center-l_branch.target_center))
                        this_P = Plm(this_cos, n,0)
                        rho_inv = 1/rho
                        integral_local += (r * rho_inv)^n * rho_inv * this_P
                    end
                end
            end
        end
        integral_local *= A * dx*dy*dz / (8*lx*ly*lz)

    elseif method == "loop"
        rvec = x_t - l_branch.target_center
        r = norm(rvec)
        err = 0.0
        for i_body in m_branch.bodies_index
            x_s = system[i_body,FastMultipole.Position()]
            q = system[i_body,FastMultipole.Strength()]
            rho_vec = x_s - l_branch.target_center
            rho = norm(rho_vec)
            rho_inv = 1/rho
            cos_gamma = dot(rvec, rho_vec) / (r*rho)
            for n in P+1:20
                if typeof(system) <: VortexParticles
                    err += dot(q,rvec)/(r*r) * n * (r*rho_inv)^n * rho_inv * Plm(cos_gamma, n, 0)
                else
                    err += q * (r*rho_inv)^n * rho_inv * Plm(cos_gamma, n, 0)
                end
            end
        end

        integral_local = err

    elseif method == "loop2"

        rvec = x_t - l_branch.target_center
        r, θt, ϕt = FastMultipole.cartesian_to_spherical(rvec)
        cos_theta_t = cos(θt)
        err = 0.0 + 0im
        for i_body in m_branch.bodies_index
            x_s = system[i_body,FastMultipole.Position()]
            q = system[i_body,FastMultipole.Strength()]
            rho_vec = x_s - l_branch.target_center
            rho, θs, ϕs = FastMultipole.cartesian_to_spherical(rho_vec)
            rho_inv = 1/rho
            cos_theta_s = cos(θs)
            for n in P+1:20
                for m in -n:n
                    if typeof(system) <: VortexParticles
                        q = dot(q,rvec)/(r*r) * n
                    end
                    err += q * (-1)^(n+m) * Float64(factorial(big(n-abs(m)))/factorial(big(n+abs(m)))) * (r*rho_inv)^n * rho_inv * Plm(cos_theta_s, n, abs(m)) * Plm(cos_theta_t, n, abs(m)) * exp(im*m*(phi_t-phi_s))
                end
            end
        end

        integral_local = err

    end

    return integral_local
end

function get_iθ(θr,nθ)
    Δθ = pi/nθ
    θ = 0.0
    for i in 1:nθ
        if θ >= θr
            θ-θr > θr-θ+Δθ && (return i-1)
            return i
        end
        θ += Δθ
    end
    pi-θr > θr-pi+Δθ && (return nθ-1)
    return nθ
end

function get_iϕ(ϕr,nϕ)
    Δϕ = pi/nϕ
    # ϕr < 0.0 && (ϕr = -ϕr)
    ϕr = abs(ϕr)
    ϕ = 0.0
    for i in 1:nϕ
        if ϕ >= ϕr
            ϕ-ϕr > ϕr-ϕ+Δϕ && (return i-1)
            return i
        end
        ϕ += Δϕ
    end
    pi-ϕr > ϕr-pi+Δϕ && (return nϕ-1)
    return nϕ
end

function integrate_multipole(x_t, m_branch, n)
    s = mean(m_branch.target_box) * 2
    dx = s / 200
    x0 = m_branch.source_center[1] - s*0.5 + dx * 0.5
    y0 = m_branch.source_center[2] - s*0.5 + dx * 0.5
    z0 = m_branch.source_center[3] - s*0.5 + dx * 0.5
    r_vec = x_t - m_branch.source_center
    r = norm(r_vec)
    r_inv = 1/r
    integral = 0.0
    for x in range(x0, step=dx, length=200)
        for y in range(y0, step=dx, length=200)
            for z in range(z0, step=dx, length=200)
                rho_vec = SVector{3}(x,y,z) - m_branch.source_center
                rho = norm(rho_vec)
                cos_gamma = dot(rho_vec, r_vec) / (rho*r)
                integral += rho^n * r_inv^(n+1) * Plm(cos_gamma,n,0)
            end
        end
    end
    Q = m_branch.multipole_expansion[1,1,1]
    integral *= dx*dx*dx * Q/(s*s*s)
    return integral
end

function multipole_error_r(rvec, r, r_inv, m_branch, system, x0, y0, z0, dx, dy, dz, nx, ny, nz, lx, ly, lz, A, P; method="loop", debug=false)

    integral = 0.0
    if method == "loop"
        err = 0.0
        for i_body in m_branch.bodies_index
            x_s = system[i_body,FastMultipole.Position()]
            q = system[i_body,FastMultipole.Strength()]
            rho_vec = x_s - m_branch.source_center
            rho = norm(rho_vec)
            cos_gamma = dot(rvec, rho_vec) / (r*rho)
            for n in P+1:20
                if typeof(system) <: VortexParticles
                    err += (n+1) * dot(q,rvec) * r_inv * r_inv * (rho*r_inv)^n * r_inv * Plm(cos_gamma, n, 0)
                else
                    err += q * (rho*r_inv)^n * r_inv * Plm(cos_gamma, n, 0)
                end
            end
        end
        integral = err

    elseif method == "numerical"
        V = 0.0
        for z in range(z0, step=dz, length=nz)
            for y in range(y0, step=dy, length=ny)
                for x in range(x0, step=dx, length=nx)
                    rho = sqrt(x^2+y^2+z^2)
                    cos_gamma = dot(rvec, SVector{3}(x,y,z)) / (r * rho)
                    for n in P+1:20
                        integral += (rho * r_inv)^n * r_inv * Plm(cos_gamma, n, 0)
                    end
                end
            end
        end
        integral *= A * dx*dy*dz / (8*lx*ly*lz)

    elseif method == "numerical_corners"

        s = mean(m_branch.target_box)*2
        V = 0.0
        for z in range(z0, step=dz, length=nz)
            for y in range(y0, step=dy, length=ny)
                for x in range(x0, step=dx, length=nx)
                    rho = sqrt(x^2+y^2+z^2)
                    if rho >= s*0.5
                        cos_gamma = dot(rvec, SVector{3}(x,y,z)) / (r * rho)
                        for n in P+1:100
                            integral += (rho * r_inv)^n * r_inv * Plm(cos_gamma, n, 0)
                        end
                        V += dx*dy*dz
                    end
                end
            end
        end
        integral *= A * dx*dy*dz / (8*lx*ly*lz)

    elseif method == "numerical_squared"

        for z in range(z0, step=dz, length=nz)
            for y in range(y0, step=dy, length=ny)
                for x in range(x0, step=dx, length=nx)
                    rho = sqrt(x^2+y^2+z^2)
                    cos_gamma = dot(rvec, SVector{3}(x,y,z)) / (r * rho)
                    for n in P+1:100
                        integral += (rho * r_inv)^(2*n) * r_inv * r_inv * Plm(cos_gamma, n, 0)^2
                    end
                end
            end
        end
        integral *= A * A * dx*dy*dz / (8*lx*ly*lz)

    elseif method == "remove_legendre"
        s = mean(m_branch.target_box) * 2
        for n in P+1:100
            #=
            t1 = (s*0.5)^(n+3) / (n+3)
            t2 = f_gpt(r,s,n)
            integral += 1/r^(n+1)*sqrt(2/(2*n+1)) * (t1+t2)
            =#
            t1 = 1/r^(n+1)
            t2 = 0.0
            nint = 5000
            dθ = pi/nint
            for theta in range(dθ/2, step=dθ, length=nint)
                t2 += dθ * Plm(cos(theta),n,0) * sin(theta)
            end
            # t2 = sqrt(2/(2*n+1))
            t3 = 0.0
            rstop = s*0.5
            dρ = rstop/50
            for rho in range(dρ*0.5, stop=rstop-dρ*0.5, length=50)
                t3 += rho^n * dρ
            end
            rstart = s*0.5
            rstop = s*sqrt(2)*0.5
            dρ = (rstart-rstop)/50
            for rho in range(dρ*0.5+rstart, stop=rstop-dρ*0.5, length=50)
                t3 += rho^n * dρ * (-2+1.5*s/rho)
            end
            integral += t1*t2*t3
        end
        integral *= 2*pi*A/(8*lx*ly*lz)

    elseif method == "precomputed"
        s = mean(m_branch.target_box) * 2
        integral = 0.0

        # get polar/azimuth angles for accessing error database
        _, θr, ϕr = FastMultipole.cartesian_to_spherical(rvec)
        _, nθ, nϕ = size(FastMultipole.MULTIPOLE_INTEGRALS)
        iθ = get_iθ(θr,nθ)
        iϕ = get_iϕ(ϕr,nϕ)

        # get s^(P+4)/r^(P+1)
        s_rinv = s / r
        scalar = s * s * s_rinv
        for n in 1:P
            scalar *= s_rinv
        end

        n_max, ndiv_x = P+1, 200
        isodd(n_max) && (n_max += 1)
        # mpint = multipole_preintegration(rvec, n_max, s, ndiv_x)
        # mpint_full = multipole_full_integration(rvec, n_max, s, nx)

        # compute error
        test_integral = 0.0
        for n in P+1:n_max
            integral += FastMultipole.MULTIPOLE_INTEGRALS[n,iθ,iϕ] * scalar
            scalar *= s_rinv
            # test_integral += integrate_multipole(rvec, m_branch, n)
        end

        # integral *= A/(8*lx*ly*lz)
        integral *= A/(s*s*s)
        # println()

        iszero(m_branch.source_box) && (integral = 0.0)

    end
    return integral
end

function rectangle_error(l_branch, l_system, m_branch, m_system, P; m_method="loop", l_method="loop")
    # multipole error (integrate over source)
    nx = ny = nz = 10
    lx, ly, lz = m_branch.target_box
    dx = 2*lx/nx
    dy = 2*ly/ny
    dz = 2*lz/nz
    x0 = -lx+dx*0.5
    y0 = -ly+dy*0.5
    z0 = -lz+dz*0.5

    x_t = closest_point(l_branch, m_branch)
    rvec = x_t - m_branch.source_center
    r = norm(rvec)
    r_inv = 1/r
    A = m_branch.multipole_expansion[1,1,1]

    integral = 0.0
    pts = get_points(l_branch) # absolute coordinates of all corners and centers
    i_pt = 0
    for (i,pt) in enumerate(pts)
        r_t2 = pt - m_branch.source_center
        r_t2n = norm(r_t2)
        r_t2inv = 1/r_t2n
        this_integral = multipole_error_r(r_t2, r_t2n, r_t2inv, m_branch, m_system, x0, y0, z0, dx, dy, dz, nx, ny, nz, lx, ly, lz, A, P; method=m_method)
        integral = max(integral, abs(this_integral))
        if integral == abs(this_integral)
            i_pt = i
        end
    end

    # local error (integrate over source)
    # x_t = closest_corner(l_branch, m_branch)
    # @assert x_t in pts "x_t not found in pts"
    max_local = 0.0
    this_pt = pts[1]
    for x_t2 in pts
        integral_local = local_error_r(x_t2, l_branch, m_branch, m_system, x0, y0, z0, dx, dy, dz, nx, ny, nz, lx, ly, lz, A, P; method=l_method)
        max_local = max(abs(integral_local), max_local)
        if max_local == abs(integral_local)
            this_pt = x_t2
        end
    end

    m_err = integral
    l_err = max_local

    return m_err/4/pi, l_err/4/pi
end

function test_error_from_m2l_list(system; expansion_order=5, multipole_threshold=0.5, leaf_size=30, shrink_recenter=true)
    tree = Tree(system; expansion_order, leaf_size, shrink_recenter)
    m2l_list, direct_list = build_interaction_lists(tree.branches, tree.branches, tree.leaf_index, multipole_threshold, true, true, true, UnequalSpheres(), expansion_order)
    leaf_sizes = Int64[]
    for b in tree.branches
        push!(leaf_sizes, length(b.bodies_index))
    end
    errs_mp, errs_l, errs_o = Float64[], Float64[], Float64[]
    errs_mp_v, errs_l_v, errs_o_v = Float64[], Float64[], Float64[]
    ubs_mp_us, ubs_l_us = Float64[], Float64[]
    ubs_mp_uus, ubs_l_uus = Float64[], Float64[]
    ubs_mp_ub, ubs_l_ub = Float64[], Float64[]
    ubs_mp_uub, ubs_l_uub = Float64[], Float64[]
    errs_mp_uc, errs_l_uc = Float64[], Float64[]
    errs_mp_uc_v, errs_l_uc_v = Float64[], Float64[]
    errs_mp_r, errs_l_r = Float64[], Float64[]
    i_m2l = 1
    for (i_target, i_source, P) in m2l_list
        # println("\n========== i_m2l = $i_m2l ==========\n")
        if i_m2l == 5
            FastMultipole.DEBUG[] = true
        end
        local_branch = tree.branches[i_target]
        multipole_branch = tree.branches[i_source]
        e_mp, e_l, e_o, e_mp_v, e_l_v, e_o_v, ub_mp_us, ub_l_us, ub_mp_uus, ub_l_uus, ub_mp_ub, ub_l_ub, ub_mp_uub, ub_l_uub, e_mp_uc, e_l_uc, ev_mp_uc, ev_l_uc = test_error(local_branch, system, multipole_branch, system, P)
        # settings
        e_mp_r, e_l_r = rectangle_error(local_branch, system, multipole_branch, system, P; m_method="loop", l_method="loop")

        if i_m2l == 5
            FastMultipole.DEBUG[] = false
        end
        i_m2l += 1
        push!(errs_mp, e_mp)
        push!(errs_l, e_l)
        push!(errs_o, e_o)
        push!(errs_mp_v, e_mp_v)
        push!(errs_l_v, e_l_v)
        push!(errs_o_v, e_o_v)
        push!(ubs_mp_us, ub_mp_us)
        push!(ubs_l_us, ub_l_us)
        push!(ubs_mp_uus, ub_mp_uus)
        push!(ubs_l_uus, ub_l_uus)
        push!(ubs_mp_ub, ub_mp_ub)
        push!(ubs_l_ub, ub_l_ub)
        push!(ubs_mp_uub, ub_mp_uub)
        push!(ubs_l_uub, ub_l_uub)
        push!(errs_mp_uc, e_mp_uc)
        push!(errs_l_uc, e_l_uc)
        push!(errs_mp_uc_v, ev_mp_uc)
        push!(errs_l_uc_v, ev_l_uc)
        push!(errs_mp_r, e_mp_r)
        push!(errs_l_r, e_l_r)
    end

    return errs_mp, errs_l, errs_o, errs_mp_v, errs_l_v, errs_o_v, ubs_mp_us, ubs_l_us, ubs_mp_uus, ubs_l_uus, ubs_mp_ub, ubs_l_ub, ubs_mp_uub, ubs_l_uub, errs_mp_uc, errs_l_uc, errs_mp_uc_v, errs_l_uc_v, errs_mp_r, errs_l_r
end



#--- run tests ---#

n_max, nθ, nϕ = 20, 80, 10
nx = 100
s = 1.0
# const MULTIPOLE_INTEGRALS = multipole_preintegration(n_max, s, nx, nθ, nϕ)
θ_max, nθ = pi/2, 50
# const LOCAL_INTEGRALS = local_preintegration(n_max, θ_max, nθ)

#function expected_legendre(n_max, cos0, dcos

# result
n_bodies = 2000
Random.seed!(123)
bodies = rand(7,n_bodies)
bodies[4,:] .= 0.0
# bodies[5,:] .-= 0.5 # coulombic interaction (positive and negative charges)
system = Gravitational(bodies)

# vortex system
# system = generate_vortex(123, n_bodies)

expansion_order, multipole_threshold, leaf_size = 1, 0.5, 260

#=
leaf_size, multipole_threshold = 25,0.5
FastMultipole.direct!(system)
velocity_direct = deepcopy(system.velocity_stretching)
reset!(system)
FastMultipole.fmm!(system; expansion_order=10, multipole_threshold, leaf_size, lamb_helmholtz=true)
velocity_fmm = deepcopy(system.velocity_stretching)
println("\nP=10:")
@show mean(abs.(velocity_direct[1:3,:] .- velocity_fmm[1:3,:]))
reset!(system)
FastMultipole.fmm!(system; expansion_order=5, multipole_threshold, leaf_size, lamb_helmholtz=true)
velocity_fmm = deepcopy(system.velocity_stretching)
println("\nP=5:")
@show mean(abs.(velocity_direct[1:3,:] .- velocity_fmm[1:3,:]))
reset!(system)
FastMultipole.fmm!(system; expansion_order=1, multipole_threshold, leaf_size, lamb_helmholtz=true)
velocity_fmm = deepcopy(system.velocity_stretching)
println("\nP=1:")
@show mean(abs.(velocity_direct[1:3,:] .- velocity_fmm[1:3,:]))
=#

# results
means_local = Float64[]
means_multipole = Float64[]
Ps = 1:17
for expansion_order in Ps
    errs_mp, errs_l, errs_o, errs_mp_v, errs_l_v, errs_o_v, ubs_mp_us, ubs_l_us, ubs_mp_uus, ubs_l_uus, ubs_mp_ub, ubs_l_ub, ubs_mp_uub, ubs_l_uub, errs_mp_uc, errs_l_uc, errs_mp_uc_v, errs_l_uc_v, errs_mp_r, errs_l_r = test_error_from_m2l_list(system; expansion_order, multipole_threshold, leaf_size, shrink_recenter=true) # leaf_size=260
    push!(means_local, output_stuff(errs_l ./ errs_l_r)[1])
    push!(means_multipole, output_stuff(errs_mp ./ errs_mp_r)[1])
end

fig = figure("error with expansion order")
fig.add_subplot(111,xlabel="expansion order", ylabel="mean potential error")
ax = fig.get_axes()[0]
ax.plot(Ps, means_multipole)
ax.plot(Ps, means_local)
ax.legend(["multipole error", "local error"])


#=
@testset "error: multipole" begin

#--- test multipole error ---#

v_mp = errs_mp ./ errs_mp_r
me_mp, ma_mp, mi_mp, st_mp = output_stuff(v_mp)
v_mp2 = errs_mp ./ errs_mp_uc
me_mp2, ma_mp2, mi_mp2, st_mp2 = output_stuff(v_mp2)

@test isapprox(me_mp, me_mp2; rtol=0.5)
@test isapprox(ma_mp, ma_mp2; rtol=0.5)
@test isapprox(mi_mp, mi_mp2; rtol=0.5)
@test isapprox(st_mp, st_mp2; rtol=0.5)

end

@testset "error: local" begin

#--- test local error ---#

v_l = errs_l ./ errs_l_r
me_l, ma_l, mi_l, st_l = output_stuff(v_l)
v_l2 = errs_l ./ errs_l_uc
me_l2, ma_l2, mi_l2, st_l2 = output_stuff(v_l2)

@test isapprox(me_l, me_l2; rtol=0.5)
@test isapprox(ma_l, ma_l2; rtol=0.5)
@test isapprox(mi_l, mi_l2; rtol=0.5)
@test isapprox(st_l, st_l2; rtol=0.5)

end
=#

# explore multipole preintegration
rvec = SVector{3}(5*sqrt(3),0,0)
n_max = 10
s = 1.0
nx = 100
res = multipole_preintegration(rvec, n_max, s, nx)

# explore scaling as we change s
res2 = multipole_preintegration(rvec, n_max, s*3, nx)

# try different polar angles
new_rvec = SVector{3}(5.0,0,5.0)
new_rvec = new_rvec ./ norm(new_rvec) * 5 * sqrt(3)
res3 = multipole_preintegration(new_rvec, n_max, s, nx)

# try different azimuthal angles
new_rvec = SVector{3}(5.0,5.0,0.0)
new_rvec = new_rvec ./ norm(new_rvec) * 5 * sqrt(3)
res4 = multipole_preintegration(new_rvec, n_max, s, nx)

# check convergence
# new_nx = 1000
# res5 = multipole_preintegration(rvec, n_max, s, new_nx)
# @show res5 .- res

# new_nx = 2000
# res6 = multipole_preintegration(rvec, n_max, s, new_nx)
# @show res6 ./ res

function vary_polar_angle(res, itheta, iϕ)
    fig = figure("vary polar angle")
    fig.clear()
    fig.add_subplot(111,xlabel="n",ylabel="integral")
    ax = fig.get_axes()[0]
    iθ = 16
    for iθ in itheta
        ax.plot(1:size(res,1), abs.(res[:,iθ,iϕ]))
    end
    dθ = pi / (size(res,2)-1)
    ax.legend([L"\theta="*"$((round(i*dθ*180/pi, digits=1)))"*L"^\circ" for i in 0:size(res,2)-1 if i+1 in itheta])
    ax.set_ylim([1e-7,1e-1])
    ax.set_yscale("log")

    fig2 = figure("vary polar angle 2")
    fig2.clear()
    fig2.add_subplot(111,xlabel="theta",ylabel="integral")
    ax2 = fig2.get_axes()[0]
    for n in 1:8#size(res,1)
        ax2.plot(range(0,stop=pi,length=size(res,2)), res[n,:,iϕ])# ./ maximum(abs.(res[n,:,iϕ])))
    end
    ax2.legend(["n=$i" for i in 1:size(res,1)])
    # ax2.set_yscale("log")

    return fig, fig2
end

itheta = collect(1:10)
iϕ = 5
fig = vary_polar_angle(FastMultipole.MULTIPOLE_INTEGRALS, itheta, iϕ)

function vary_azimuth(res, iθ, iphi)
    fig = figure("vary azimuth")
    fig.clear()
    fig.add_subplot(111,xlabel="n",ylabel="integral")
    ax = fig.get_axes()[0]
    for iϕ in iphi
        ax.scatter(1:size(res,1), abs.(res[:,iθ,iϕ]))
    end
    dϕ = 2*pi / (size(res,3)-1)
    ax.legend([L"\theta="*"$((round(i*dϕ*180/pi, digits=1)))"*L"^\circ" for i in 0:size(res,3)-1 if i+1 in iphi])
    ax.set_ylim([1e-7,1e-1])
    ax.set_yscale("log")
    return fig
end

iphi = collect(1:size(FastMultipole.MULTIPOLE_INTEGRALS,3))
iθ = 4
fig_az = vary_azimuth(FastMultipole.MULTIPOLE_INTEGRALS, iθ, iphi)

#=

# x axis adjacent
#R = SVector{3}(2.0,0,0)
#r_multipole, r_local = 0.5, 0.5
#nx_source_bodies, nx_target_bodies = 10, 10
#expansion_order = 3
#system, direct, local_potential = test_rms_error(R, r_multipole, r_local, nx_source_bodies, nx_target_bodies, expansion_order)

# z axis adjacent
R = SVector{3}(2.0,0,0)
r_multipole, r_local = 0.5, 0.5
nx_source_bodies, nx_target_bodies = 10, 10
expansion_order = 5
system, direct, multipole_potential, local_potential, e_mp, e_l, e_o = test_rms_error(R, r_multipole, r_local, nx_source_bodies, nx_target_bodies, expansion_order)
e_o_rms_rel = sqrt(Base.sum((e_o) .^ 2)/1000)
s_error_rel = sqrt(solvason_error(expansion_order))
s_error_abs = direct.potential[1,:] .* s_error_rel
s_diff = abs.(e_o) ./ abs.(s_error_abs[1001:2000])
@show maximum(s_diff) mean(s_diff) median(s_diff) minimum(s_diff)

s = sqrt(1.5^2+0.5^2*2)

r_err = r_averaged_error(expansion_order, s)

=#
