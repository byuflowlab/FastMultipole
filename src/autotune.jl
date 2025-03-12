#------- choose max expansion order -------#

function leaf_size_converged(optargs, leaf_size_source)
    converged = true
    for i in eachindex(optargs.leaf_size_source)
        converged = abs(optargs.leaf_size_source[i] - leaf_size_source[i]) < 0.1 * leaf_size_source[i]
    end
    return converged
end

tune_fmm!(system; kwargs...) = tune_fmm!(system, system; kwargs...)

function tune_fmm!(target_systems, source_systems; kwargs...)
    # promote arguments to Tuples
    target_systems = to_tuple(target_systems)
    source_systems = to_tuple(source_systems)

    return tune_fmm!(target_systems, source_systems; kwargs...)
end

function tune_fmm!(target_systems::Tuple, source_systems::Tuple;
    ε_abs=nothing,
    expansion_order=4, leaf_size_source=default_leaf_size(source_systems),
    max_expansion_order=20, max_iter=10,
    multipole_thresholds=range(0.3, stop=0.8, step=0.1),
    verbose=true, kwargs...
)

    if verbose
        println("\n#======= Begin FastMultipole.tune_fmm() =======#")
    end

    #--- save best parameters ---#

    t_fmm_best = Inf
    expansion_order_best = max_expansion_order
    leaf_size_source_best = leaf_size_source
    multipole_threshold_best = multipole_thresholds[1]

    #--- preallocate cache ---#

    t_fmm = @elapsed _, cache, _ = fmm!(target_systems, source_systems;
                       expansion_order=1, leaf_size_source,
                       nearfield=false, farfield=false, self_induced=false,
                       tune=true, update_target_systems=false
                      )

    #--- error tolerance selected ---#

    for multipole_threshold in multipole_thresholds

        if verbose
            println("\nmultipole_threshold = $multipole_threshold...")
        end

        # initial fmm! call with max_expansion_order to get expansion order
        t_fmm = @elapsed optargs, _, _, _, m2l_list, _, _, error_success = fmm!(target_systems, source_systems;
                                                                                expansion_order=isnothing(ε_abs) ? expansion_order : max_expansion_order,
                                                                                leaf_size_source, multipole_threshold,
                                                                                ε_abs, kwargs..., cache...,
                                                                                tune=true, update_target_systems=false,
                                                                               )

        # in case error is not satisfied
        if !error_success
            println("\terror tolerance not satisfied for max expansion order P=$max_expansion_order;")
            println("\tskipping this multipole_threshold...")
            continue
        end

        # track the best parameters for this multipole_threshold
        this_t_fmm_best = Inf
        this_expansion_order_best = max_expansion_order
        this_leaf_size_source_best = leaf_size_source

        # iterate to (loose) convergence
        i = 1
        for _ in 1:max_iter

            # in case m2l list is empty (direct calculation is probably best)
            if length(m2l_list) == 0 # likely won't get much better
                println("\tM2L list is empty; \n\tending iterations for this multipole_threshold...")
                if t_fmm < t_fmm_best
                    this_t_fmm_best = t_fmm_best = t_fmm
                    this_leaf_size_source_best = leaf_size_source_best = get_n_bodies_vec(source_systems)
                    this_expansion_order_best = expansion_order_best = 1
                    multipole_threshold_best = multipole_threshold
                end
                break # next multipole_threshold
            end

            # save the expansion order
            expansion_order = optargs.expansion_order

            # predict optimal leaf size
            t_fmm = @elapsed optargs, cache, _, _, m2l_list, _, _, error_success = fmm!(target_systems, source_systems;
                                                                                        expansion_order,
                                                                                        leaf_size_source, multipole_threshold,
                                                                                        ε_abs, kwargs..., cache...,
                                                                                        tune=true
                                                                                       )

            # save leaf size
            leaf_size_source = optargs.leaf_size_source

            # benchmark and check for convergence
            t_fmm = @elapsed optargs, cache, _, _, m2l_list, _, _, error_success = fmm!(target_systems, source_systems;
                                                                                        expansion_order,
                                                                                        leaf_size_source, multipole_threshold,
                                                                                        ε_abs, kwargs..., cache...,
                                                                                        tune=true
                                                                                       )
            if error_success # (loosely) converged

                # check if this is our best yet
                if t_fmm < t_fmm_best
                    this_t_fmm_best = t_fmm_best = t_fmm
                    this_leaf_size_source_best = leaf_size_source_best = leaf_size_source
                    this_expansion_order_best = expansion_order_best = expansion_order
                    multipole_threshold_best = multipole_threshold
                end

                break # move to the next multipole_threshold
            end

            i += 1
        end

        if verbose
            println("\n\tBest Parameters: ")
            println("\t\tleaf_size_source:    ", this_leaf_size_source_best)
            println("\t\texpansion_order:     ", this_expansion_order_best)
            println("\t\tmultipole_threshold: ", multipole_threshold)
            println("\t\tcost:                $this_t_fmm_best seconds")
            println("\t\titerations:          ", i)
            println("\t\toptargs:             ", optargs)
        end

    end

    if verbose
        println("\nFinished autotune!")
        println("\nParameters: ")
        println("\tleaf_size_source:    ", leaf_size_source_best)
        println("\texpansion_order:     ", expansion_order_best)
        println("\tmultipole_threshold: ", multipole_threshold_best)
        println("\tcost:                $t_fmm_best seconds")
        println("\n#===============================================#\n")
    end

    tuned_params = (
                    leaf_size_source = leaf_size_source_best,
                    expansion_order = expansion_order_best,
                    multipole_threshold = multipole_threshold_best,
                   )

    return tuned_params, cache
end

