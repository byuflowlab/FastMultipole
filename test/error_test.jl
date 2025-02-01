using Statistics, Random, LegendrePolynomials
using PythonPlot
using LaTeXStrings
using LinearAlgebra
using DelimitedFiles
using FastMultipole

include("../test/gravitational.jl")
include("../test/bodytolocal.jl")
include("../test/evaluate_multipole.jl")

function predicted_errors(tree, m2l_list, system, error_method, expansion_order, lamb_helmholtz::Val)
    # preallocate containers
    multipole_errs_hat = zeros(length(m2l_list))
    local_errs_hat = zeros(length(m2l_list))

    # predict errors for each m2l interaction
    for (i,(i_target, i_source)) in enumerate(m2l_list)
        # extract branches
        multipole_branch = tree.branches[i_source]
        local_branch = tree.branches[i_target]
        multipole_errs_hat[i] = multipole_error(local_branch, multipole_branch, expansion_order, error_method, lamb_helmholtz)
        local_errs_hat[i] = local_error(local_branch, multipole_branch, expansion_order, error_method, lamb_helmholtz)
    end

    return multipole_errs_hat, local_errs_hat
end

function expansion_errors(tree, m2l_list, system, expansion_order, lamb_helmholtz::Val)

    # preallocate containers
    multipole_errs_mean = zeros(length(m2l_list))
    local_errs_mean = zeros(length(m2l_list))
    overall_errs_mean = zeros(length(m2l_list))
    multipole_errs_max = zeros(length(m2l_list))
    local_errs_max = zeros(length(m2l_list))
    overall_errs_max = zeros(length(m2l_list))

    # create derivatives switch
    derivatives_switch = DerivativesSwitch(true, true, false, system)
    if typeof(system) <: Gravitational
        body_type = Point{Source}
    elseif typeof(system) <: VortexParticles
        body_type = Point{Vortex}
    else
        throw("body type not defined for expansion_errors")
    end

    for (i,(i_target, i_source)) in enumerate(m2l_list)
        # extract branches
        local_branch = tree.branches[i_target]
        multipole_branch = tree.branches[i_source]

        # zero velocity influence
        for i in local_branch.bodies_index
            system[i,Velocity()] = zero(system[i,Velocity()])
        end

        # direct velocity influence
        FastMultipole._direct!(system, local_branch.bodies_index, derivatives_switch, system, multipole_branch.bodies_index)
        v_direct = [system[i,Velocity()] for i in local_branch.bodies_index]

        # zero velocity influence
        for i in local_branch.bodies_index
            system[i,Velocity()] = zero(system[i,Velocity()])
        end

        # multipole velocity error
        evaluate_multipole!(system, local_branch, multipole_branch, multipole_branch.harmonics, expansion_order, lamb_helmholtz, derivatives_switch)
        v_multipole = [system[i,Velocity()] for i in local_branch.bodies_index]
        multipole_error = norm.(v_direct - v_multipole)
        @assert length(multipole_error) == length(local_branch.bodies_index)
        multipole_errs_mean[i] = mean(multipole_error)
        multipole_errs_max[i] = maximum(multipole_error)

        # zero velocity influence
        for i in local_branch.bodies_index
            system[i,Velocity()] = zero(system[i,Velocity()])
        end

        # local velocity error
        local_branch.local_expansion .= zero(eltype(local_branch.local_expansion))
        for j in multipole_branch.bodies_index
            Δx = system[j,Position()] - local_branch.target_center
            strength = system[j,Strength()]
            body_to_local_point!(body_type, local_branch.local_expansion, local_branch.harmonics, Δx, strength, expansion_order)
        end
        FastMultipole.evaluate_local!(system, local_branch, expansion_order, lamb_helmholtz, derivatives_switch)
        v_local = [system[i,Velocity()] for i in local_branch.bodies_index]
        local_error = norm.(v_direct - v_local)
        @assert length(local_error) == length(local_branch.bodies_index)
        local_errs_mean[i] = mean(local_error)
        local_errs_max[i] = maximum(local_error)

        # zero velocity influence
        for i in local_branch.bodies_index
            system[i,Velocity()] = zero(system[i,Velocity()])
        end

        # overall velocity error
        local_branch.local_expansion .= zero(eltype(local_branch.local_expansion))
        FastMultipole.multipole_to_local!(local_branch, multipole_branch, expansion_order, lamb_helmholtz, nothing)
        FastMultipole.evaluate_local!(system, local_branch, expansion_order, lamb_helmholtz, derivatives_switch)
        v_overall = [system[i,Velocity()] for i in local_branch.bodies_index]
        overall_error = norm.(v_direct - v_overall)
        @assert length(overall_error) == length(local_branch.bodies_index)
        overall_errs_mean[i] = mean(overall_error)
        overall_errs_max[i] = maximum(overall_error)

    end

    return multipole_errs_mean, local_errs_mean, overall_errs_mean, multipole_errs_max, local_errs_max, overall_errs_max
end


function output_stuff(v; verbose=true)
    if verbose
        @show mean(abs.(v)) maximum(abs.(v)) minimum(abs.(v)) std(abs.(v))
    end
    return mean(abs.(v)), maximum(abs.(v)), minimum(abs.(v)), std(abs.(v))
end

#--- Gravitational System ---#

# create system
seed = 123
n_bodies = 1000
system = generate_gravitational(seed, n_bodies; radius_factor=0.0)

# generate multipole expansions and m2l-list
leaf_size, expansion_order, multipole_threshold, lamb_helmholtz = 10, 5, 0.5, false
tree = Tree(system; expansion_order=expansion_order+1, leaf_size, shrink_recenter=true)
m2l_list, direct_list, derivatives_switch = fmm!(tree, system; multipole_threshold, expansion_order=expansion_order+1, unsort_bodies=false, lamb_helmholtz)

# compute expansion errors
multipole_errs_mean, local_errs_mean, overall_errs_mean, multipole_errs_max, local_errs_max, overall_errs_max = expansion_errors(tree, m2l_list, system, expansion_order, Val(lamb_helmholtz))

# predicted errors
error_method = FastMultipole.RotatedCoefficients()
multipole_errs_hat, local_errs_hat = predicted_errors(tree, m2l_list, system, error_method, expansion_order, Val(lamb_helmholtz))

# predicted over actual
v_mp = multipole_errs_hat ./ multipole_errs_max
v_l = local_errs_hat ./ local_errs_max
v_o = (multipole_errs_hat .+ local_errs_hat) ./ overall_errs_max

# plot results

