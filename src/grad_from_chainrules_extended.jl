macro grad_from_chainrules_extended(iip_inds, fcall)
    Meta.isexpr(fcall, :call) && length(fcall.args) >= 2 || # meta stuff I do not want to touch
        error("`@grad_from_chainrules_iip` has to be applied to a function signature")
    f = esc(fcall.args[1])
    xs = map(fcall.args[2:end]) do x
        if Meta.isexpr(x, :(::))
            if length(x.args) == 1 # ::T without var name
                return :($(gensym())::$(esc(x.args[1])))
            else # x::T
                @assert length(x.args) == 2
                return :($(x.args[1])::$(esc(x.args[2])))
            end
        else
            return x
        end
    end
    args_l, args_r, args_track, args_fixed, arg_types, kwargs = ReverseDiff._make_fwd_args(f, xs) # should be fine as is
    return quote
        #typeof($iip_inds) == Tuple || 
        #    error("in-place inputs must be referenced in a tuple! If there are no in-place inputs, use an empty tuple.")
        $f($(args_l...)) = ReverseDiff.track($(args_r...))
        function ReverseDiff.track($(args_track...))
            args = ($(args_fixed...),)
            tp = ReverseDiff.tape(args...)

            # get in-place and out-of-place output_values
            if length($iip_inds) > 0
                input_values = map(ReverseDiff.value, args)
                output_values_ip = map(ReverseDiff.value, args)
                #=for a in args[$iip_inds...] # any inputs modified in-place need to be tracked. If they aren't, weird behavior occurs.
                    if !(ReverseDiff.hastape(a))
                        ReverseDiff.track!(a,tp)
                    end
                end=#
                output_values, back = ChainRulesCore.rrule($f, output_values_ip...; $kwargs...)
                ReverseDiff.value!.(args[$iip_inds...],output_values_ip[$iip_inds...])
            else
                output_values, back = ChainRulesCore.rrule($f, map(ReverseDiff.value, args)...; $kwargs...)
                input_values = nothing
                output_values_ip = nothing
            end
            # broadcast tracking if multiple outputs.
            if output_values !== nothing
                output = length(output_values) > 1 ? broadcast((_output)->ReverseDiff.track(_output, tp),output_values) : ReverseDiff.track(output_values, tp)
            else
                output = nothing
            end
            if input_values !== nothing
                input = length(input_values) > 1 ? broadcast((_input)->ReverseDiff.track(_input, tp),input_values) : ReverseDiff.track(input_values, tp)
            else
                input = nothing
            end
            closure(cls_args...; cls_kwargs...) = ChainRulesCore.rrule($f, map(ReverseDiff.value, cls_args)...; cls_kwargs...)
            ReverseDiff.record!(tp, ReverseDiff.SpecialInstruction, $f, args, output, (back, closure, $kwargs, output_values_ip),)
            return output
        end

        @noinline function ReverseDiff.special_reverse_exec!(instruction::ReverseDiff.SpecialInstruction{typeof($f), <:Tuple{$(arg_types...)}})
            input = instruction.input
            back = instruction.cache[1]
            output = instruction.output
            modified_input = instruction.cache[4]
            if output !== nothing
                if modified_input !== nothing
                    back_derivs = (ReverseDiff.deriv.(output)..., ReverseDiff.deriv.(input[$iip_inds...])...)
                else
                    back_derivs = ReverseDiff.deriv.(output)
                end
            else
                if modified_input !== nothing
                    back_derivs = ReverseDiff.deriv.(input[$iip_inds...])
                else
                    error("no output for function $f")
                end
            end
            if length($iip_inds) > 0
                back_output = back(back_derivs) # for some reason the in-place derivatives end up in the input rather than the output. this accounts for that.
                ReverseDiff.unseed!(input[$iip_inds...]) # we no longer need the incoming derivatives, so they are cleared to avoid corrupting the new values
            else
                back_output = back(back_derivs)
            end
            input_derivs = back_output[2:end]
            @assert input_derivs isa Tuple
            ReverseDiff._add_to_deriv!.(input, input_derivs)
            ReverseDiff.unseed!(output)
            return nothing
        end

        @noinline function ReverseDiff.special_forward_exec!(instruction::ReverseDiff.SpecialInstruction{typeof($f), <:Tuple{$(arg_types...)}})
            output, input = instruction.output, instruction.input
            ReverseDiff.pull_value!.(input)
            pullback = instruction.cache[2]
            kwargs = instruction.cache[3]
            out_value = pullback(input...; kwargs...)[1]
            if out_value !== nothing
                ReverseDiff.value!.(output, out_value)
            end
            if instruction.cache[4] !== nothing
                ReverseDiff.value!.(input[$iip_inds...], instruction.cache[4][$iip_inds...])
            end
            return nothing
        end
    end
end