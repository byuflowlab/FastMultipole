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

function expansion_errors(tree::FastMultipole.Tree{TF,<:Any,<:Any}, m2l_list, system, expansion_order, lamb_helmholtz::Val) where TF

    # preallocate containers
    multipole_errs_mean = zeros(length(m2l_list))
    local_errs_mean = zeros(length(m2l_list))
    overall_errs_mean = zeros(length(m2l_list))
    multipole_errs_max = zeros(length(m2l_list))
    local_errs_max = zeros(length(m2l_list))
    overall_errs_max = zeros(length(m2l_list))

    # create derivatives switch
    derivatives_switch = DerivativesSwitch(true, true, false, (system,))
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
        for i in local_branch.bodies_index[1]
            system[i,Velocity()] = zero(SVector{3,TF})
        end

        # direct velocity influence
        FastMultipole._direct!(system, local_branch.bodies_index[1], derivatives_switch, system, multipole_branch.bodies_index[1])
        v_direct = [system[i,Velocity()] for i in local_branch.bodies_index[1]]

        # zero velocity influence
        for i in local_branch.bodies_index[1]
            system[i,Velocity()] = zero(system[i,Velocity()])
        end

        # multipole velocity error
        evaluate_multipole!((system,), local_branch, multipole_branch, multipole_branch.harmonics, expansion_order, lamb_helmholtz, derivatives_switch)
        v_multipole = [system[i,Velocity()] for i in local_branch.bodies_index[1]]
        multipole_error = norm.(v_direct - v_multipole)
        @assert length(multipole_error) == length(local_branch.bodies_index[1])
        multipole_errs_mean[i] = mean(multipole_error)
        multipole_errs_max[i] = maximum(multipole_error)

        # zero velocity influence
        for i in local_branch.bodies_index[1]
            system[i,Velocity()] = zero(system[i,Velocity()])
        end

        # local velocity error
        local_branch.local_expansion .= zero(eltype(local_branch.local_expansion))
        for j in multipole_branch.bodies_index[1]
            Δx = system[j,Position()] - local_branch.target_center
            strength = system[j,Strength()]
            body_to_local_point!(body_type, local_branch.local_expansion, local_branch.harmonics, Δx, strength, expansion_order)
        end
        FastMultipole.evaluate_local!((system,), local_branch, expansion_order, lamb_helmholtz, derivatives_switch, SVector{1}(true))
        v_local = [system[i,Velocity()] for i in local_branch.bodies_index[1]]
        local_error = norm.(v_direct - v_local)
        @assert length(local_error) == length(local_branch.bodies_index[1])
        local_errs_mean[i] = mean(local_error)
        local_errs_max[i] = maximum(local_error)

        # zero velocity influence
        for i in local_branch.bodies_index[1]
            system[i,Velocity()] = zero(system[i,Velocity()])
        end

        # overall velocity error
        local_branch.local_expansion .= zero(eltype(local_branch.local_expansion))
        FastMultipole.multipole_to_local!(local_branch, multipole_branch, expansion_order, lamb_helmholtz, nothing)
        FastMultipole.evaluate_local!((system,), local_branch, expansion_order, lamb_helmholtz, derivatives_switch, SVector{1}(true))
        v_overall = [system[i,Velocity()] for i in local_branch.bodies_index[1]]
        overall_error = norm.(v_direct - v_overall)
        @assert length(overall_error) == length(local_branch.bodies_index[1])
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

# expansion order
expansion_order = 1

# create system
seed = 123
n_bodies = 1000
system = generate_gravitational(seed, n_bodies; radius_factor=0.0)

# generate multipole expansions and m2l-list
leaf_size, multipole_threshold, lamb_helmholtz = SVector{1}(10), 0.5, false
tree = Tree((system,); expansion_order=expansion_order+2, leaf_size, shrink_recenter=true)
is_target=SVector{1}(true)
is_source=SVector{1}(true)
tree, m2l_list, direct_list, derivatives_switches = fmm!((system,), tree, is_target, is_source; multipole_threshold, expansion_order=expansion_order+1, unsort_bodies=false, lamb_helmholtz)

# compute expansion errors
multipole_errs_mean, local_errs_mean, overall_errs_mean, multipole_errs_max, local_errs_max, overall_errs_max = expansion_errors(tree, m2l_list, system, expansion_order, Val(lamb_helmholtz))

# predicted errors
error_method = FastMultipole.RotatedCoefficients()
multipole_errs_hat, local_errs_hat = predicted_errors(tree, m2l_list, system, error_method, expansion_order, Val(lamb_helmholtz))

# actual errors
reset!(system)
fmm!(system; expansion_order, leaf_size, multipole_threshold)
fmm_velocity = [SVector{3}(system[i,Velocity()]) for i in 1:get_n_bodies(system)]
reset!(system)
direct!(system)
direct_velocity = [SVector{3}(system[i,Velocity()]) for i in 1:get_n_bodies(system)]
error_velocity = [norm(v1 - v2) for (v1,v2) in zip(fmm_velocity, direct_velocity)]
mean_error = mean(error_velocity)
std_error = std(error_velocity)
max_error = maximum(error_velocity)
min_error = minimum(error_velocity)


# print outputs
println("\n#------- Checking Error -------#")
println("\tlength of m2l_list:    $(length(m2l_list))")
println("\tlength of direct_list: $(length(direct_list))")
println("\tmean error (manual):   $(mean_error)")
println("\tstd error (manual):   $(std_error)")
println("\tmax error (manual):   $(max_error)")
println("\tmin error (manual):   $(min_error)")
println()

# predicted over actual
@assert sum(isnan.(multipole_errs_hat)) == 0
@assert sum(isnan.(multipole_errs_max)) == 0
@assert sum(isnan.(local_errs_hat)) == 0
@assert sum(isnan.(local_errs_max)) == 0
@assert sum(isnan.(overall_errs_max)) == 0

function clean_nans!(v)
    for i in eachindex(v)
        if abs(v[i]) < 1e-12
            v[i] = 1e-12
        end
    end
end

clean_nans!(multipole_errs_hat)
clean_nans!(multipole_errs_max)
clean_nans!(local_errs_hat)
clean_nans!(local_errs_max)
clean_nans!(overall_errs_max)

v_mp = multipole_errs_hat ./ multipole_errs_max
v_l = local_errs_hat ./ local_errs_max
v_o = (multipole_errs_hat .+ local_errs_hat) ./ overall_errs_max

@assert sum(isnan.(v_mp)) == 0
@assert sum(isnan.(v_l)) == 0
@assert sum(isnan.(v_o)) == 0

# plot results
fig = figure("histogram")
fig.clear()
fig.add_subplot(111, xlabel="prediction over actual", ylabel="density")
ax = fig.get_axes()[0]
ax.hist(log10.(abs.(v_mp)), bins=35, label="multipole")
ax.hist(log10.(abs.(v_l)), bins=35, label="local")
ax.hist(log10.(abs.(v_o)), bins=35, label="total")
