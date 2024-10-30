macro grad_from_chainrules_multiple_returns(fcall)
    Meta.isexpr(fcall, :call) && length(fcall.args) >= 2 || # meta stuff I do not want to touch
        error("`@grad_from_chainrules` has to be applied to a function signature")
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
        $f($(args_l...)) = ReverseDiff.track($(args_r...))
        function ReverseDiff.track($(args_track...))
            args = ($(args_fixed...),)
            tp = ReverseDiff.tape(args...)
            output_value, back = ChainRulesCore.rrule($f, map(ReverseDiff.value, args)...; $kwargs...)
            #output = ReverseDiff.track(output_value, tp)
            output = map(_o->ReverseDiff.track(_o,tp),output_value)
            closure(cls_args...; cls_kwargs...) = ChainRulesCore.rrule($f, map(ReverseDiff.value, cls_args)...; cls_kwargs...)
            ReverseDiff.record!(
                tp,
                ReverseDiff.SpecialInstruction,
                $f,
                args,
                output,
                (back, closure, $kwargs),
            )
            return output
        end

        @noinline function ReverseDiff.special_reverse_exec!(instruction::ReverseDiff.SpecialInstruction{typeof($f), <:Tuple{$(arg_types...)}})
            output = instruction.output
            input = instruction.input
            back = instruction.cache[1]
            back_output = back(ReverseDiff.deriv.(output))
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
            ReverseDiff.value!(output, out_value)
            return nothing
        end
    end
end

# for functions where the first argument is pre-allocated memory where the argument doesn't matter. 
macro grad_from_chainrules_cached(fcall)
    Meta.isexpr(fcall, :call) && length(fcall.args) >= 2 || # meta stuff I do not want to touch
        error("`@grad_from_chainrules` has to be applied to a function signature")
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
        $f($(args_l...)) = ReverseDiff.track($(args_r...))
        function ReverseDiff.track($(args_track...))
            args = ($(args_fixed...),)
            tp = ReverseDiff.tape(args...)
            output_value, back = ChainRulesCore.rrule($f, map(ReverseDiff.value, args)...; $kwargs...)
            #output = ReverseDiff.track(output_value, tp)
            output = map(_o->ReverseDiff.track(_o,tp),output_value)
            closure(cls_args...; cls_kwargs...) = ChainRulesCore.rrule($f, map(ReverseDiff.value, cls_args)...; cls_kwargs...)
            ReverseDiff.record!(
                tp,
                ReverseDiff.SpecialInstruction,
                $f,
                args,
                output,
                (back, closure, $kwargs),
            )
            return output
        end

        @noinline function ReverseDiff.special_reverse_exec!(instruction::ReverseDiff.SpecialInstruction{typeof($f), <:Tuple{$(arg_types...)}})
            output = instruction.output
            input = instruction.input
            back = instruction.cache[1]
            back_output = back(ReverseDiff.deriv.(output))
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
            ReverseDiff.value!(output, out_value)
            return nothing
        end
    end
end