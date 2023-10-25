# This file contains a set of methods that allow custom rrules with complex inputs/outputs to function.
#    Note that the inputs/outputs of the overall function still need to be purely real.

# Complex{TrackedReal} needs a bunch of new functions defined to act like a tracked complex number. These functions are inherited from
#    ReverseDiff's tracked.jl file.

# If you get weird ReverseDiff type-related errors (especially if TrackedReals/TrackedArrays end up somewhere unintended), it is possible that
#    an access rule is bugged or missing from this list.

# getters
@inline ReverseDiff.istracked(::Complex{<:ReverseDiff.TrackedReal}) = true
@inline ReverseDiff.value(z::Complex{<:ReverseDiff.TrackedReal}) = Complex{Float64}(real(z).value,imag(z).value)
@inline ReverseDiff.deriv(z::Complex{<:ReverseDiff.TrackedReal}) = Complex{Float64}(real(z).deriv,imag(z).deriv)
@inline ReverseDiff.valtype(::Complex{ReverseDiff.TrackedReal{V}}) where {V} = V
@inline ReverseDiff.valtype(::Type{Complex{ReverseDiff.TrackedReal{V,D,O}}}) where {V,D,O} = V
@inline ReverseDiff.valtype(::Complex{ReverseDiff.TrackedReal{V,D,O}}) where {V,D,O} = V
@inline ReverseDiff.derivtype(::Complex{ReverseDiff.TrackedReal{V,D}}) where {V,D} = D
@inline ReverseDiff.derivtype(::Type{Complex{ReverseDiff.TrackedReal{V,D,O}}}) where {V,D,O} = D
@inline ReverseDiff.derivtype(::Complex{ReverseDiff.TrackedReal{V,D,O}}) where {V,D,O} = D
@inline ReverseDiff.origintype(::Complex{ReverseDiff.TrackedReal{V,D,O}}) where {V,D,O} = O
@inline ReverseDiff.origintype(::Complex{Type{ReverseDiff.TrackedReal{V,D,O}}}) where {V,D,O} = O
@inline ReverseDiff.hasorigin(z::Complex{<:ReverseDiff.TrackedReal}) = real(z).index === imag(z).index && real(z).index !== ReverseDiff.NULL_INDEX
@inline ReverseDiff.hastape(z::Complex{<:ReverseDiff.TrackedReal}) = ReverseDiff.tape(z) !== ReverseDiff.NULL_TAPE
@inline ReverseDiff.tape(z::Complex{<:ReverseDiff.TrackedReal}) = real(z).tape

# additional getters for Complex{TrackedReal} in an array
@inline ReverseDiff.istracked(::AbstractArray{Complex{T}}) where {T} = T <: ReverseDiff.TrackedReal || !(isconcretetype(T))
@inline ReverseDiff.value(x::AbstractArray{C}) where {C <: Complex} = ReverseDiff.istracked(x) ? map(ReverseDiff.value, x) : x

# setters
@inline ReverseDiff.value!(t::Complex{<:ReverseDiff.TrackedReal}, v::Real) = (real(t).value = v; nothing)
@inline ReverseDiff.value!(t::Complex{<:ReverseDiff.TrackedReal}, z::Complex) = (t.real.value = real(z); t.imag.value = imag(z); nothing) # not 100% sure this won't break.
@inline ReverseDiff.deriv!(t::Complex{<:ReverseDiff.TrackedReal}, v::Real) = (real(t).deriv = v; nothing)
@inline ReverseDiff.deriv!(t::Complex{<:ReverseDiff.TrackedReal}, z::Complex) = (t.real.deriv = real(z); t.imag.deriv = imag(z); nothing) # not 100% sure this won't break.

ReverseDiff.pull_value!(t::Complex{<:ReverseDiff.TrackedReal}) = (ReverseDiff.hasorigin(t) && ReverseDiff.value!(t, ReverseDiff.value(t.origin)[t.index]); nothing)
ReverseDiff.pull_deriv!(t::Complex{<:ReverseDiff.TrackedReal}) = (ReverseDiff.hasorigin(t) && ReverseDiff.deriv!(t, ReverseDiff.deriv(t.origin)[t.index]); nothing)
ReverseDiff.push_deriv!(t::Complex{<:ReverseDiff.TrackedReal}) = (ReverseDiff.hasorigin(t) && (t.origin.deriv[t.index] = ReverseDiff.deriv(t)); nothing)

# seed/unseed... I don't think these are ever called, though.
function ReverseDiff.seed!(t::Complex{<:ReverseDiff.TrackedReal},i)
    if i == 1
        (t.real.deriv[i] = one(ReverseDiff.derivtype(t)); nothing)
    elseif i == 2
        (t.imag.deriv[i] = one(ReverseDiff.derivtype(t)); nothing)
    else
        error("index i was $(i), which was not 1 or 2!")
    end
end

ReverseDiff.seed!(z::Complex,i) = i == 1 ? seed!(z.real) : seed!(z.imag)
ReverseDiff.seed!(z::Complex) = error("there should be an index here...")

function ReverseDiff.unseed!(t::Complex{<:ReverseDiff.TrackedReal})
    real(t).deriv = zero(ReverseDiff.derivtype(t)); ReverseDiff.push_deriv!(t)
    imag(t).deriv = zero(ReverseDiff.derivtype(t)); ReverseDiff.push_deriv!(t)
    return nothing
end

ReverseDiff.capture(t::Complex) = ReverseDiff.istracked(t) ?  map!(ReverseDiff.capture, similar(t), t) : copy(t)
ReverseDiff.capture(t::Complex{<:ReverseDiff.TrackedReal}) = ifelse(ReverseDiff.hastape(t), t, ReverseDiff.value(t))
function Base.similar(z::Complex{TR}) where TR <: ReverseDiff.TrackedReal
    Complex{TR}(zero(realtype(z)),zero(realtype(z)))
end

# track/track! functions. These are pretty central to unpacking larger structures and propegating derivatives.
function ReverseDiff.track!(t::Complex{<:ReverseDiff.TrackedReal{D,D,Nothing}}, x::Complex, tp::ReverseDiff.InstructionTape) where D
    t.real = ReverseDiff.track(real(x),D,tp)
    t.imag = ReverseDiff.track(imag(x),D,tp)
    return t
end
realtype(::ComplexF64) = Float64 # because eltype(ComplexF64) == ComplexF64 rather than Float64
ReverseDiff.track(x::Complex, tp::ReverseDiff.InstructionTape = ReverseDiff.InstructionTape()) = ReverseDiff.track(x, realtype(x), tp)
ReverseDiff.track(x::Complex, ::Type{D}, tp::ReverseDiff.InstructionTape = ReverseDiff.InstructionTape()) where D = Complex(ReverseDiff.TrackedReal(real(x), zero(D), tp),ReverseDiff.TrackedReal(imag(x), zero(D), tp))


# these catch errors but should never be called.
#=function ReverseDiff.track(a::Float64,b::Bool)
    @show a b
    error(" ")
end=#
#=
function Base.convert(::Type{T1}, t::T2) where {T1<:Complex{ReverseDiff.TrackedReal},T2<:ReverseDiff.TrackedReal}
    error(" ")
    V1, D1, O1 = ReverseDiff.valtype(T1), ReverseDiff.derivtype(T1), ReverseDiff.origintype(T1)
    tp = ReverseDiff.tape(t)
    out = ReverseDiff.TrackedReal{V1,D1,O1}(_convert(V1, ReverseDiff.value(t)), _convert(D1, ReverseDiff.deriv(t)), tp)
    ReverseDiff.record!(tp, ReverseDiff.SpecialInstruction, convert, t, out)
    return out
end
function Base.convert(::Type{T},z::C) where {T <: Complex{ReverseDiff.TrackedReal}, C}
    error(" ")
end
function Base.convert(::Type{C},z::T) where {T <: ReverseDiff.TrackedReal, C}
    error(" ")
end=#


# these next functions make sure derivatives get propagated. Unimplemented ones still return an error so that I know if one of them is needed.
function ReverseDiff.increment_deriv!(t::Complex,x)
    error("t: $t\tx: $x")
end
function ReverseDiff.increment_deriv!(t::Complex,x,i)
    error("t: $t\tx: $x\ti: $i")
end

function ReverseDiff._add_to_deriv!(a::Complex,b)
    if ReverseDiff.istracked(a)
        real(a).deriv += real(b)
        imag(a).deriv += imag(b)
    end
    return nothing
end
#=
# not sure if these next two are needed, but they shouldn't break anything.
function ReverseDiff._add_to_deriv!(d1::AbstractArray{Complex{<:ReverseDiff.TrackedReal}}, d2::ReverseDiff.AbstractThunk)
    ReverseDiff.increment_deriv!(d1, ReverseDiff.unthunk(d2))
end
function ReverseDiff._add_to_deriv!(d1::AbstractArray{Complex{<:ReverseDiff.TrackedReal}}, d2)
    ReverseDiff.increment_deriv!(d1, d2)
end

@inline ReverseDiff.increment_deriv!(t::ReverseDiff.TrackedArray, x::Complex, i) = (t.deriv[i] += x; nothing)
@inline ReverseDiff.increment_deriv!(t::AbstractArray, x::Complex, i) = ReverseDiff.increment_deriv!(t[i], x)
=#
################ test case; will be moved to the actual unit tests once I've cleaned it up.

cmult(A,B) = A*B
function ChainRulesCore.rrule(::typeof(cmult),A,B)

    function C_pullback(C̄)
        return NoTangent(), B'*C̄, A'*C̄
    end
    return cmult(A,B), C_pullback

end
ReverseDiff.@grad_from_chainrules cmult(A::Complex{<:ReverseDiff.TrackedReal},B::Complex{<:ReverseDiff.TrackedReal})
ReverseDiff.@grad_from_chainrules cmult(A::Complex{<:ReverseDiff.TrackedReal},B)
ReverseDiff.@grad_from_chainrules cmult(A,B::Complex{<:ReverseDiff.TrackedReal})

function test_complex_multiplication()

    c = [1.0]
    function f(x)
        z1 = x[1] + cmult(x[1],1.0*im)
        z2 = 1.0 + 2.0*im
        z = cmult(z1,z2)
        return abs(z) + x[1]
    end
    fc = f(c)
    dfdc = ReverseDiff.gradient(f,c)
    return fc,dfdc

end

