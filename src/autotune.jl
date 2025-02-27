#------- choose max expansion order -------#

function max_expansion_order(ε_abs, max_velocity, min_r, multipole_threshold)

    # convert ε_abs to ε_rel
    ε_rel = ε_abs / (max_velocity * 10)

    # convert from multipole_threshold to c
    c = 2.0 / multipole_threshold - 1.0

    # choose max expansion order (see Pringle)
    ρ_over_r = 1 / c

    for P in 1:20

    end

    return P
end

tune_fmm!(system; optargs...) = tune_fmm!(system, system; optargs...)

function tune_fmm!(target_systems, source_systems; optargs...)
    # promote arguments to Tuples
    target_systems = to_tuple(target_systems)
    source_systems = to_tuple(source_systems)

    return tune_fmm!(target_systems, source_systems; optargs...)
end

function tune_fmm!(target_systems::Tuple, source_systems::Tuple;
    expansion_order=5,
    ε_abs=nothing,
    shrink_recenter=true, lamb_helmholtz=true,
    max_expansion_order=20,
    multipole_thresholds=range(0.3, stop=0.7, step=0.1),
    leaf_size_source=default_leaf_size(source_systems),
    scalar_potential=true, velocity=true, velocity_gradient=false,
    verbose=true
)

    #--- predict maximum velocity ---#

    if verbose
        println("\n#======= Begin FastMultipole.autotune!() =======#")
        # println("\n\tpredicting max velocity...")
    end

    # reset!(target_systems)
    # fmm!(target_systems, source_systems;
    #      expansion_order=1, leaf_size_source=SVector{length(source_systems)}()
    #      lamb_helmholtz, ε_abs, tune=true,
    #      scalar_potential, velocity, velocity_gradient)

    # max_velocity = zero(eltype(target_systems[1]))
    # for system in target_systems
    #     for i in 1:get_n_bodies(system)
    #         max_velocity = max(max_velocity, norm(system[i, Velocity()]))
    #     end
    # end

    #--- preallocate buffers ---#

    optargs, _ = fmm!(target_systems, source_systems; expansion_order=1)
    source_buffers = optargs.source_buffers
    target_buffers = optargs.target_buffers
    source_small_buffers = optargs.source_small_buffers
    target_small_buffers = optargs.target_small_buffers

    #--- loop over multipole_threshold ---#

    # preallocated storage for when we find the best case
    leaf_size_sources = Vector{SVector{length(source_systems),Int}}(undef, length(multipole_thresholds))
    expansion_orders = Vector{Int}(undef, length(multipole_thresholds))
    ts_fmm = @MVector zeros(length(multipole_thresholds))

    for (i_mt, multipole_threshold) in enumerate(multipole_thresholds)
        if verbose
            println("\n\ttuning multipole_threshold = $multipole_threshold...")
        end

        #--- determine max expansion_order ---#

        if !isnothing(ε_abs)
            expansion_order = max_expansion_order
            # expansion_order = max_expansion_order(ε_abs, max_velocity, multipole_threshold)
        end

        #--- tune expansion order and leaf size ---#

        # get leaf size
        _ = @elapsed optargs, _ =
            fmm!(target_systems, source_systems;
                source_buffers, target_buffers,
                source_small_buffers, target_small_buffers,
                multipole_threshold, leaf_size_source,
                expansion_order, ε_abs, lamb_helmholtz,
                scalar_potential, velocity, velocity_gradient,
                tune=true,
            )

            leaf_size_source = optargs.leaf_size_source

        # get expansion order now that leaf_size has been chosen
        _ = @elapsed optargs, _ =
            fmm!(target_systems, source_systems;
                source_buffers, target_buffers,
                source_small_buffers, target_small_buffers,
                multipole_threshold, leaf_size_source,
                expansion_order, ε_abs, lamb_helmholtz,
                scalar_potential, velocity, velocity_gradient,
                tune=true,
            )
            expansion_order = optargs.expansion_order

        #--- benchmark with tuned parameters ---#

        t_fmm = @elapsed fmm!(target_systems, source_systems;
                source_buffers, target_buffers,
                source_small_buffers, target_small_buffers,
                multipole_threshold, leaf_size_source,
                expansion_order, ε_abs, lamb_helmholtz,
                scalar_potential, velocity, velocity_gradient,
                tune=true,
            )

        leaf_size_sources[i_mt] = leaf_size_source
        expansion_orders[i_mt] = expansion_order
        ts_fmm[i_mt] = error_success ? t_fmm : Inf # if error was not successfully constrained, penalize this case so it is not chosen

        if verbose
            tuned_params = (leaf_size_source => leaf_size_source,
                    expansion_order => expansion_order,
                    multipole_threshold => multipole_threshold)
            println("\t\tparameters: ", tuned_params)
            println("\t\tcost = $t_fmm seconds")
            println("\t\tlength(m2l_list): ", length(m2l_list))
            println("\t\tlength(direct_list): ", length(direct_list))
        end

    end

    #--- compile tuned parameters ---#

    _, i = findmin(ts_fmm)
    leaf_size_source = leaf_size_sources[i]
    expansion_order = expansion_orders[i]
    multipole_threshold = multipole_thresholds[i]

    tuned_params = (
                    source_buffers = source_buffers,
                    target_buffers = target_buffers,
                    source_small_buffers = source_small_buffers,
                    target_small_buffers = target_small_buffers,
                    leaf_size_source = leaf_size_source,
                    expansion_order = expansion_order,
                    multipole_threshold = multipole_threshold
                   )

    #--- return ---#

    if verbose
        println("\n\tfinished!")
        println("\n\tparameters: ", tuned_params)
        println("\n#===============================================#\n")
    end


    return tuned_params
end

