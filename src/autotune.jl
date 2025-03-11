#------- choose max expansion order -------#

#=
tune_fmm_simple!(system; kwargs...) = tune_fmm_simple!(system, system; kwargs...)

function tune_fmm_simple!(target_systems, source_systems; kwargs...)
    # promote arguments to Tuples
    target_systems = to_tuple(target_systems)
    source_systems = to_tuple(source_systems)

    return tune_fmm_simple!(target_systems, source_systems; kwargs...)
end

function tune_fmm_simple!(target_systems::Tuple, source_systems::Tuple; max_iter=10, max_expansion_order=20,
        ε_abs, expansion_order=5, leaf_size_source=default_leaf_size(source_systems), multipole_thresholds=range(0.3, stop=0.8, step=0.1),
        verbose, kwargs...,
        ΔP=5
    )

    # variables to persist between iterations
    t_fmm_best = Inf
    leaf_size_source_best = leaf_size_source
    expansion_order_best = expansion_order
    multipole_threshold_best = multipole_thresholds[1]

    # get cache
    optargs, cache, _ = fmm!(target_systems, source_systems; ε_abs, expansion_order=1, kwargs...)

    # tune
    for multipole_threshold in multipole_thresholds
        ΔP = ΔP0
        for i in 1:max_iter
            t_fmm = @elapsed optargs, _, _, _, _, _, _, error_success = fmm!(target_systems, source_systems; ε_abs, expansion_order, multipole_threshold, leaf_size_source, kwargs...)

            # check if this is the best guess
            if t_fmm < t_fmm_best && error_success
                expansion_order_best = expansion_order
                leaf_size_source_best = leaf_size_source
                multipole_threshold_best = multipole_threshold
            end

            # update guesses
            leaf_size_source = optargs.leaf_size_source
            if expansion_order > optargs.expansion_order
                expansion_order = optargs.expansion_order
            else
                expansion_order = expansion_order + ΔP
            end
            if error_success
                ΔP = 1
            end
        end
        ΔP0 = 2
    end
end
=#

function leaf_size_converged(optargs, leaf_size_source)
    converged = true
    for i in eachindex(optargs.leaf_size_source)
        converged = abs(optargs.leaf_size_source[i] - leaf_size_source[i]) < 0.1 * leaf_size_source[i]
    end
    return converged
end

function autotune!(target_systems, source_systems, max_iter, max_expansion_order;
        ε_abs, expansion_order, leaf_size_source, multipole_threshold,
        verbose, kwargs, cache,
        α_LS=0.7, ΔP0=4, # α_P=0.6,
    )

    # initialize variables to persist between iterations
    t_fmm_best = Inf
    expansion_order_best = max_expansion_order
    leaf_size_source_best = leaf_size_source
    ΔP = ΔP0

    # begin tuning
    for i in 1:max_iter

        #--- run fmm! ---#

        t_fmm = @elapsed optargs, _, _, _, m2l_list, direct_list, _, error_success =
            fmm!(target_systems, source_systems;
                 ε_abs,
                 multipole_threshold, leaf_size_source, expansion_order,
                 tune=true, update_target_systems=false,
                 kwargs..., cache...
                )

        #--- save best parameters in case we terminate early ---#

        if t_fmm < t_fmm_best && error_success
            t_fmm_best = t_fmm
            expansion_order_best = expansion_order
            leaf_size_source_best = leaf_size_source
        end

        #--- check for empty M2L list ---#

        if length(m2l_list) == 0 # probably faster to use direct only
            if verbose
                println("\n\tempty m2l_list detected: it's probably faster to evaluate interactions without FMM")
            end
            return leaf_size_source, expansion_order, t_fmm
        end

        #--- error tolerance is satisfied ---#

        if error_success

            # once error tolerance is satisfied, we probably don't want to change expansion order as quickly
            ΔP = 1

            if optargs.expansion_order == expansion_order # this has been true 2 iterations in a row and is probably converged

                # check leaf size convergence
                if leaf_size_converged(optargs, leaf_size_source)

                    if verbose
                        println("\n\tConverged Parameters: ")
                        println("\t\tleaf_size_source:    ", leaf_size_source)
                        println("\t\texpansion_order:     ", expansion_order)
                        println("\t\tmultipole_threshold: ", multipole_threshold)
                        println("\t\tlength(m2l_list):    ", length(m2l_list))
                        println("\t\tlength(direct_list): ", length(direct_list))
                        println("\t\tcost:                $t_fmm seconds")
                    end

                    return leaf_size_source, expansion_order, t_fmm
                else
                    # set leaf size for next loop
                    # leaf_size_source = optargs.leaf_size_source

                    # try relaxation
                    leaf_size_source = Int.(ceil.( (1-α_LS) * leaf_size_source + α_LS * optargs.leaf_size_source))
                end

            else # expansion order is not converged

                # set expansion order for next loop
                expansion_order = optargs.expansion_order

                # try fine tuning
                # Δ = optargs.expansion_order > expansion_order ? ΔP : max(-ΔP0, optargs.expansion_order - expansion_order)
                # Δ = optargs.expansion_order > expansion_order ? ΔP : -ΔP
                # expansion_order += Δ
                # try relaxation
                # expansion_order = max(1, min(max_expansion_order, Int(ceil((1-α_P) * expansion_order + α_P * optargs.expansion_order))))

                # set leaf size for next loop (ACTUALLY, DON'T DO THIS TO AVOID CONVERGENCE ISSUES)
                # leaf_size_source = optargs.leaf_size_source

            end

        #--- error tolerance is not satisfied ---#

        else

            # max expansion order has been reached, leaf size is converged, and error is still not satisfied
            if expansion_order == max_expansion_order # && leaf_size_converged(optargs, leaf_size_source)

                if verbose
                    println("\n\terror tolerance not reached for max_expansion_order=$max_expansion_order;")
                    println("\tskipping this multipole_threshold\n")
                end

                t_fmm = Inf
                return leaf_size_source, expansion_order, t_fmm

            else

                # increment expansion order
                expansion_order = min(max_expansion_order, expansion_order + ΔP)

                # update leaf size
                # leaf_size_source = optargs.leaf_size_source

                # try relaxation
                leaf_size_source = Int.(ceil.( (1-α_LS) * leaf_size_source + α_LS * optargs.leaf_size_source))

            end
        end

        if verbose
            println("\tfinished iteration $i")
            println("\t\tParameters:")
            println("\t\tleaf_size_source:    ", leaf_size_source)
            println("\t\texpansion_order:     ", expansion_order)
            println("\t\tmultipole_threshold: ", multipole_threshold)
            println("\t\tcost:                $t_fmm seconds")
            println("\t\terror_success:       ", error_success)
        end
    end

    if verbose
        println("\n\tMax iterations reached! Returning best parameters:")
        println("\t\tleaf_size_source:    ", leaf_size_source_best)
        println("\t\texpansion_order:     ", expansion_order_best)
        println("\t\tmultipole_threshold: ", multipole_threshold_best)
        println("\t\tcost:                $t_fmm_best seconds")
    end

    return leaf_size_source_best, expansion_order_best, t_fmm_best
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

    #--- predict maximum velocity ---#

    if verbose
        println("\n#======= Begin FastMultipole.autotune!() =======#")
    end

    #--- preallocate buffers ---#

    _, cache, _ = fmm!(target_systems, source_systems;
                      expansion_order=1, update_target_systems=false, tune=true,
                      kwargs...
                     )

    #--- loop over multipole_threshold ---#

    # preallocated storage for when we find the best case
    leaf_size_sources = Vector{SVector{length(source_systems),Int}}(undef, length(multipole_thresholds))
    expansion_orders = Vector{Int}(undef, length(multipole_thresholds))
    ts_fmm = @MVector zeros(length(multipole_thresholds))

    #--- other variables ---#
    leaf_size_source = default_leaf_size(source_systems)

    for (i_mt, multipole_threshold) in enumerate(multipole_thresholds)
        if verbose
            println("\nTuning multipole_threshold = $multipole_threshold...\n")
        end

        leaf_size_source, expansion_order, t_fmm = autotune!(target_systems, source_systems, max_iter, max_expansion_order;
                                                             ε_abs, expansion_order, leaf_size_source, multipole_threshold,
                                                             verbose, kwargs, cache
                                                            )

        # update history
        leaf_size_sources[i_mt] = leaf_size_source
        expansion_orders[i_mt] = expansion_order
        ts_fmm[i_mt] = t_fmm

    end

    #--- compile tuned parameters ---#

    t_opt, i = findmin(ts_fmm)
    leaf_size_source = leaf_size_sources[i]
    expansion_order = expansion_orders[i]
    multipole_threshold = multipole_thresholds[i]

    tuned_params = (
                    leaf_size_source = leaf_size_source,
                    expansion_order = expansion_order,
                    multipole_threshold = multipole_threshold,
                   )

    #--- return ---#

    if verbose
        println("\nFinished autotune!")
        println("\nparameters: ")
        println("\tleaf_size_source:    ", tuned_params.leaf_size_source)
        println("\texpansion_order:     ", tuned_params.expansion_order)
        println("\tmultipole_threshold: ", tuned_params.multipole_threshold)
        println("\tcost:                $t_opt seconds")
        println("\n#===============================================#\n")
    end

    return tuned_params, cache, t_opt
end

