using Statistics, Random, LegendrePolynomials
using PythonPlot
using LaTeXStrings
using LinearAlgebra
using KernelDensity
using DelimitedFiles
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

    return mymax(e_mp), mymax(e_l), mymax(e_o), mymax(ev_mp), mymax(ev_l), mymax(ev_o)
end

function test_error_from_m2l_list(system; error_method=UniformCubes(), expansion_order=Dynamic(20,1e-6), multipole_threshold=0.5, leaf_size=30, shrink_recenter=true, check_dipole=Val(true))
    tree = Tree(system; expansion_order=FastMultipole.get_Pmax(expansion_order), leaf_size, shrink_recenter)
    m2l_list, direct_list = build_interaction_lists(tree.branches, tree.branches, tree.leaf_index, multipole_threshold, true, true, true, UnequalSpheres(), expansion_order, check_dipole)
    leaf_sizes = Int64[]
    for b in tree.branches
        push!(leaf_sizes, length(b.bodies_index))
    end
    errs_mp, errs_l, errs_o = Float64[], Float64[], Float64[]
    errs_mp_v, errs_l_v, errs_o_v = Float64[], Float64[], Float64[]
    i_m2l = 1
    print_every = Int(ceil(length(m2l_list) / 100))
    @show length(m2l_list)

    for (i_target, i_source, P) in m2l_list

        # println("\n========== i_m2l = $i_m2l ==========")
        if i_m2l % print_every == 0
            println("\tProgress: $(i_m2l/length(m2l_list)*100) %")
        end
        if i_m2l == -1
            FastMultipole.DEBUG[] = true
        end
        local_branch = tree.branches[i_target]
        multipole_branch = tree.branches[i_source]
        Δx, Δy, Δz = local_branch.target_center - multipole_branch.source_center

        # run this once so the multipole expansion is generated
        e_mp, e_l, e_o, e_mp_v, e_l_v, e_o_v = test_error(local_branch, system, multipole_branch, system, 1)

        # choose expansion order
        P = FastMultipole.get_P(Δx, Δy, Δz, local_branch, multipole_branch, expansion_order, error_method, check_dipole)

        # compute error
        e_mp, e_l, e_o, e_mp_v, e_l_v, e_o_v = test_error(local_branch, system, multipole_branch, system, P)
        if i_m2l == -1
            FastMultipole.DEBUG[] = false
        end

        i_m2l += 1

        # save outputs
        push!(errs_mp, e_mp)
        push!(errs_l, e_l)
        push!(errs_o, e_o)
        push!(errs_mp_v, e_mp_v)
        push!(errs_l_v, e_l_v)
        push!(errs_o_v, e_o_v)
    end

    return errs_mp, errs_l, errs_o, errs_mp_v, errs_l_v, errs_o_v
end



#--- run tests ---#

function run_tests(coulombic; fast=false)
    # results

    if fast
        n_bodies = 2000
    else
        n_bodies = 10000
    end
    Random.seed!(123)
    bodies = rand(7,n_bodies)
    bodies[4,:] .= 0.0
    if coulombic
        bodies[5,:] .-= 0.5 # coulombic interaction (positive and negative charges)
    end
    system = Gravitational(bodies)

    # vortex system
    # system = generate_vortex(123, n_bodies)

    for ε_tol in (1e-3, 1e-6, 1e-9)

        expansion_order = Dynamic(20,ε_tol)
        if fast
            multipole_threshold, leaf_size = 0.5, 260
        else
            multipole_threshold, leaf_size = 0.5, 100
        end

        # UniformCubes, check_dipole=true
        error_method = UniformCubes()
        check_dipole = Val(true)
        errs_mp, errs_l, errs_o, errs_mp_v, errs_l_v, errs_o_v = test_error_from_m2l_list(system; error_method, expansion_order, multipole_threshold, leaf_size, shrink_recenter=true, check_dipole)

        # UniformCubes, check_dipole=false
        error_method = UniformCubes()
        check_dipole = Val(false)
        errs_mp_ndp, errs_l_ndp, errs_o_ndp, errs_mp_v_ndp, errs_l_v_ndp, errs_o_v_ndp = test_error_from_m2l_list(system; error_method, expansion_order, multipole_threshold, leaf_size, shrink_recenter=true, check_dipole)

        # UniformCubesVelocity, check_dipole=true
        error_method = UniformCubesVelocity()
        check_dipole = Val(true)
        errs_mp_v, errs_l_v, errs_o_v, errs_mp_v_v, errs_l_v_v, errs_o_v_v = test_error_from_m2l_list(system; error_method, expansion_order, multipole_threshold, leaf_size, shrink_recenter=true, check_dipole)

        # UniformCubesVelocity, check_dipole=false
        error_method = UniformCubesVelocity()
        check_dipole = Val(false)
        errs_mp_v_ndp, errs_l_v_ndp, errs_o_v_ndp, errs_mp_v_v_ndp, errs_l_v_v_ndp, errs_o_v_v_ndp = test_error_from_m2l_list(system; error_method, expansion_order, multipole_threshold, leaf_size, shrink_recenter=true, check_dipole)

        # UnequalBoxes, check_dipole=true
        error_method = UnequalBoxes()
        check_dipole = Val(true)
        errs_mp_ub, errs_l_ub, errs_o_ub, errs_mp_v_ub, errs_l_v_ub, errs_o_v_ub = test_error_from_m2l_list(system; error_method, expansion_order, multipole_threshold, leaf_size, shrink_recenter=true, check_dipole)

        # UnequalBoxes, check_dipole=false
        error_method = UnequalBoxes()
        check_dipole = Val(false)
        errs_mp_ub_ndp, errs_l_ub_ndp, errs_o_ub_ndp, errs_mp_v_ub_ndp, errs_l_v_ub_ndp, errs_o_v_ub_ndp = test_error_from_m2l_list(system; error_method, expansion_order, multipole_threshold, leaf_size, shrink_recenter=true, check_dipole)


        #------- get data for paper -------#

        # checking multipole plus local error equals total error (close...)
        v_error = (abs.(errs_mp) .+ abs.(errs_l)) ./ abs.(errs_o)
        output_stuff(v_error; verbose=true)

        # log scale comparison
        v = log10.(abs.(errs_o ./ ε_tol))
        v_ndp = log10.(abs.(errs_o_ndp ./ ε_tol))
        v_v = log10.(abs.(errs_o_v_v ./ ε_tol))
        v_v_ndp = log10.(abs.(errs_o_v_v_ndp ./ ε_tol))
        v_ub = log10.(abs.(errs_o_ub ./ ε_tol))
        v_ub_ndp = log10.(abs.(errs_o_ub_ndp ./ ε_tol))

        # kde partitioning
        npoints = 2^7
        u = kde(v; npoints)
        u_ndp = kde(v_ndp; npoints)
        u_v = kde(v_v; npoints)
        u_v_ndp = kde(v_v_ndp; npoints)
        u_ub = kde(v_ub; npoints)
        u_ub_ndp = kde(v_ub_ndp; npoints)

        # make preliminary plots
        fig = figure("dynamic expansion order: tol=1e$(Int(round(log10(ε_tol);sigdigits=1)))")
        fig.clear()
        fig.add_subplot(111, xlabel=L"\log_{10} ε/ε_{tol}", ylabel="density")
        ax = fig.get_axes()[0]
        ax.plot(u.x, u.density)
        ax.plot(u_ndp.x, u_ndp.density, "--")
        ax.plot(u_v.x, u_v.density)
        ax.plot(u_v_ndp.x, u_v_ndp.density, "--")
        ax.plot(u_ub.x, u_ub.density)
        ax.plot(u_ub_ndp.x, u_ub_ndp.density, "--")
        # ax.set_xticks([0.0,1.0,2.0,3.0], labels=["1", "10", "100", "1000"])
        # ax.set_xlim([-0.5,4.0])
        legend(["potential", "potential, no dipole", "velocity", "velocity, no dipole", "setting "*L"P_n=1", "setting "*L"P_n=1"*", no dipole"])
        savefig("figs/dynamic_p_n$(n_bodies)_eps$(Int(round(log10(ε_tol);sigdigits=1)))_coulombic_$(coulombic).png")

        # save data
        dynamic_P_performance = zeros(npoints,4)
        dynamic_P_performance[:,1] .= u.x
        dynamic_P_performance[:,2] .= u.density
        header = Matrix{String}(undef,1,4)
        for i in 1:4
            header[i] = ["potential x", "potential density", "velocity magnitude x", "velocity magnitude density"][i]
        end
        writedlm("data/dynamic_p_n$(n_bodies)_eps$(Int(round(log10(0.001);sigdigits=1)))_coulombic_$(coulombic).png", vcat(header,dynamic_P_performance), ','; header=true)

    end

end

run_tests(true; fast=false)
run_tests(false; fast=false)
