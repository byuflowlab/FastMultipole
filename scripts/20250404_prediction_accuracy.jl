#=
this file is used to generate Figure 5 in the paper (examples_verify_prediction_source.pdf and examples_verify_prediction_vortex.pdf)
=#

using FastMultipole
using Statistics
using Random
using BSON

include("../test/gravitational.jl")
include("../test/vortex.jl")

function get_velocity(system::Gravitational)
    return system.potential[5:7,:]
end

function get_velocity(system::VortexParticles)
    return system.velocity_stretching[1:3,:]
end

function velocity_err!(v_fmm, v_direct, bodies_index)
    # difference
    for d in 1:3
        for (i,ib) in enumerate(bodies_index)
            v_fmm[d, i] -= v_direct[d+4, ib]
        end
    end

    # square
    for d in 1:3
        for i in eachindex(bodies_index)
            v_fmm[d, i] *= v_fmm[d, i]
        end
    end

    # sum to d=1
    for d in 2:3
        for i in eachindex(bodies_index)
            v_fmm[1, i] += v_fmm[d, i]
        end
    end
    
    # sqrt
    for i in eachindex(bodies_index)
        v_fmm[1, i] = sqrt(v_fmm[1, i])
    end
end
    

function test_accuracy(target_systems::Tuple, source_systems::Tuple, expansion_orders, multipole_threshold, lamb_helmholtz, shrink_recenter, leaf_size)
    
    # outputs
    max_errs_list = Vector{Float64}[]
    ε_mp_hat_list = Vector{Float64}[]
    ε_l_hat_list = Vector{Float64}[]
    ε_hat_list = Vector{Float64}[]

    # wrap val
    lamb_helmholtz_bool = lamb_helmholtz
    lamb_helmholtz = Val(lamb_helmholtz_bool)

    # loop over expansion orders
    for expansion_order in expansion_orders
        println("\n#----- Expansion order: $expansion_order -----#")

        # get leaf size
        optargs, cache, _ = fmm!(target_systems, source_systems; lamb_helmholtz=lamb_helmholtz_bool, expansion_order, multipole_threshold, nearfield=false, farfield=false, self_induced=false)
        # leaf_size_source = optargs.leaf_size_source
        # leaf_size_source = SVector{length(source_systems)}(leaf_size for _ in eachindex(source_systems))
        leaf_size_source = leaf_size
        leaf_size_target = FastMultipole.to_vector(minimum(leaf_size_source), length(target_systems))

        # build trees
        target_tree = Tree(target_systems, true; buffers=cache.target_buffers, small_buffers=cache.target_small_buffers, expansion_order=expansion_order+1, leaf_size=leaf_size_target, shrink_recenter)
        source_tree = Tree(source_systems, false; buffers=cache.source_buffers, small_buffers=cache.source_small_buffers, expansion_order=expansion_order+1, leaf_size=leaf_size_source, shrink_recenter)

        # get m2l list
        m2l_list, _ = FastMultipole.build_interaction_lists(target_tree.branches, source_tree.branches, leaf_size_source, multipole_threshold, true, false, false, SelfTuning())

        # upward pass
        FastMultipole.upward_pass_singlethread!(source_tree, source_systems, expansion_order+1, lamb_helmholtz)

        # preallocate containers to be reused
        weights_tmp_1 = initialize_expansion(expansion_order+1, Float64)
        weights_tmp_2 = initialize_expansion(expansion_order+1, Float64)
        Ts = zeros(Float64, FastMultipole.length_Ts(expansion_order+1))
        eimϕs = zeros(Float64, 2, expansion_order + 2)
        harmonics = initialize_harmonics(expansion_order)
        velocity_n_m = FastMultipole.initialize_velocity_n_m(expansion_order, Float64)
        derivatives_switches = DerivativesSwitch(false, true, false, target_systems)

        # target velocity container
        target_velocity = Tuple(zeros(Float64, 3, get_n_bodies(target_systems[i])) for i in 1:length(target_systems))

        # outputs
        max_errs = zeros(length(m2l_list))
        ε_mp_hat = zeros(length(m2l_list))
        ε_l_hat = zeros(length(m2l_list))

        # manual horizontal pass
        println("\n\t#--- Manual Horizontal Pass ---#\n")
        for (i,(i_target, j_source)) in enumerate(m2l_list)
            if i % 1000 == 0
                println("\t\ti = $i / $(length(m2l_list))")
            end

            #--- expansions ---#
            
            # perform M2L
            target_branch = target_tree.branches[i_target]
            target_expansion = view(target_tree.expansions, :, :, :, i_target)
            target_expansion .= 0.0
            source_branch = source_tree.branches[j_source]
            source_expansion = view(source_tree.expansions, :, :, :, j_source)
            FastMultipole.multipole_to_local!(target_expansion, target_branch, source_expansion, source_branch, weights_tmp_1, weights_tmp_2, Ts, eimϕs, FastMultipole.ζs_mag, FastMultipole.ηs_mag, FastMultipole.Hs_π2, expansion_order, lamb_helmholtz, nothing)
            
            # evaluate local expansion
            FastMultipole.reset!(cache.target_buffers)
            for (i_system, system) in enumerate(cache.target_buffers)
                FastMultipole.evaluate_local!(system, target_branch.bodies_index[i_system], harmonics, velocity_n_m, target_expansion, target_branch.target_center, expansion_order, lamb_helmholtz, derivatives_switches[i_system])
            end
            
            # save velocity
            for i_target_system in 1:length(cache.target_buffers)
                target_velocity[i_target_system][:, 1:length(target_branch.bodies_index[i_target_system])] .= cache.target_buffers[i_target_system][5:7, target_branch.bodies_index[i_target_system]]
            end

            #--- predict error ---#

            ε_mp, ε_l = FastMultipole.predict_error(target_branch, source_expansion, source_branch, weights_tmp_1, weights_tmp_2, Ts, eimϕs, FastMultipole.ζs_mag, FastMultipole.Hs_π2, expansion_order, lamb_helmholtz)
            ε_mp_hat[i] = ε_mp
            ε_l_hat[i] = ε_l

            #--- direct ---#
            
            # reset systems
            FastMultipole.reset!(cache.target_buffers)

            # loop over source systems
            for i_source_system in eachindex(source_systems)
                source_system = source_systems[i_source_system]
                source_buffer = cache.source_buffers[i_source_system]
        
                # loop over target systems
                for (i_target_system, target_system) in enumerate(cache.target_buffers)

                    # extract derivatives switch
                    derivatives_switch = derivatives_switches[i_target_system]

                    # identify sources
                    source_index = source_tree.branches[j_source].bodies_index[i_source_system]

                    # identify targets
                    target_index = target_tree.branches[i_target].bodies_index[i_target_system]

                    # compute interaction
                    direct!(target_system, target_index, derivatives_switch, source_system, source_buffer, source_index)

                end
            end

            # calculate error
            for i_target_system in 1:length(cache.target_buffers)
                velocity_err!(target_velocity[i_target_system], cache.target_buffers[i_target_system], target_branch.bodies_index[i_target_system])
                max_errs[i] = max(max_errs[i], maximum(view(target_velocity[i_target_system],1,1:length(target_branch.bodies_index[i_target_system]))))
            end
        end

        # save errors
        push!(max_errs_list, max_errs)
        push!(ε_mp_hat_list, ε_mp_hat)
        push!(ε_l_hat_list, ε_l_hat)
        push!(ε_hat_list, ε_mp_hat .+ ε_l_hat)
    end
    return max_errs_list, ε_mp_hat_list, ε_l_hat_list, ε_hat_list
end

#--- create systems ---#

println("#----- Creating systems -----#")

n_bodies = 50_000
source = generate_gravitational(123, n_bodies; strength_scale=1.0/3798.955768976926)
vortex = generate_vortex(123, n_bodies; strength_scale=1.0/10566.33461495282)

#--- check error prediction ---#

expansion_orders = 1:15
multipole_threshold = 0.5
shrink_recenter = true
leaf_size = SVector{1}(50)

lamb_helmholtz = false
max_errs_list, ε_mp_hat_list, ε_l_hat_list, ε_hat_list = test_accuracy((source,), (source,), expansion_orders, multipole_threshold, lamb_helmholtz, shrink_recenter, leaf_size)
BSON.@save "20250404_prediction_accuracy_source.bson" max_errs_list, ε_mp_hat_list, ε_l_hat_list, ε_hat_list, expansion_orders, multipole_threshold, leaf_size

lamb_helmholtz = true
max_errs_list, ε_mp_hat_list, ε_l_hat_list, ε_hat_list = test_accuracy((vortex,), (vortex,), expansion_orders, multipole_threshold, lamb_helmholtz, shrink_recenter, leaf_size)
BSON.@save "20250404_prediction_accuracy_vortex.bson" max_errs_list, ε_mp_hat_list, ε_l_hat_list, ε_hat_list, expansion_orders, multipole_threshold, leaf_size