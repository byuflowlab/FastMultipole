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
@inline ReverseDiff.origintype(::Type{Complex{ReverseDiff.TrackedReal{V,D,O}}}) where {V,D,O} = O
@inline ReverseDiff.hasorigin(z::Complex{<:ReverseDiff.TrackedReal}) = real(z).index === imag(z).index && real(z).index !== ReverseDiff.NULL_INDEX
@inline ReverseDiff.hastape(z::Complex{<:ReverseDiff.TrackedReal}) = ReverseDiff.tape(z) !== ReverseDiff.NULL_TAPE
@inline function ReverseDiff.tape(z::Complex{<:ReverseDiff.TrackedReal})
    if ReverseDiff.istracked(real(z))
        ReverseDiff.hastape(real(z)) && return ReverseDiff.tape(real(z))
    end
    if ReverseDiff.istracked(imag(z))
        ReverseDiff.hastape(imag(z)) && return ReverseDiff.tape(imag(z))
    end
    return ReverseDiff.NULL_TAPE
end

# additional getters for Complex{TrackedReal} in an array
@inline ReverseDiff.istracked(::AbstractArray{<:Complex{T}}) where {T} = T <: ReverseDiff.TrackedReal || !(isconcretetype(T))
@inline ReverseDiff.value(t::AbstractArray{<:Complex}) = ReverseDiff.istracked(t) ? map(ReverseDiff.value, t) : t
@inline ReverseDiff.deriv(t::AbstractArray{<:Complex}) = ReverseDiff.istracked(t) ? map(ReverseDiff.deriv, t) : t

# setters
@inline ReverseDiff.value!(t::Complex{<:ReverseDiff.TrackedReal}, v::Real) = (t.re.value = v; t.im.value = zero(typeof(v)); nothing)
#@inline ReverseDiff.value!(t::Complex{<:ReverseDiff.TrackedReal}, z::Complex) = (real(t).value = real(z); imag(t).value = imag(z); nothing) # not 100% sure this won't break.
@inline function ReverseDiff.value!(t::Complex{T}, z::Complex) where {T <: ReverseDiff.TrackedReal}
    #=temp = Complex{T}(real(t),imag(t))
    real(temp).value = real(z)
    imag(temp).value = imag(z)
    #real(temp).deriv = real(t).deriv
    #imag(temp).deriv = imag(t).deriv
    @show t temp z
    t = temp
    @show t temp z
    error(" ")=#
    t.re.value = real(z)
    t.im.value = imag(z)
    nothing # not 100% sure this won't break.
end
@inline ReverseDiff.deriv!(t::Complex{<:ReverseDiff.TrackedReal}, v::Real) = (t.re.deriv = v; t.im.deriv = zero(typeof(v)); nothing)
@inline ReverseDiff.deriv!(t::Complex{<:ReverseDiff.TrackedReal}, z::Complex) = (t.re.deriv = real(z); t.im.deriv = imag(z); nothing) # not 100% sure this won't break.

# additional setters for Complex{TrackedReal} in an array
#@inline ReverseDiff.value!(t::Complex{<:ReverseDiff.TrackedReal})
@inline function ReverseDiff.value!(t::AbstractArray{<:Complex},v::AbstractArray{<:Complex})
    if ReverseDiff.istracked(t)
        for i=1:length(t)
            t[i].re.value = real(v[i])
            t[i].im.value = imag(v[i])
            #(copyto!(value(t), v); nothing)
        end
    end
    #(copyto!(ReverseDiff.value(t), v); nothing)
    return nothing

end

@inline function ReverseDiff.deriv!(t::AbstractArray{<:Complex},v::AbstractArray{<:Complex})

    for i=1:length(t)
        t[i].re.value = real(v[i])
        t[i].im.value = imag(v[i])
    end

end

@inline function ReverseDiff.value!(t::Vector{<:Complex}, v::Vector{ComplexF64})
    if ReverseDiff.istracked(t)
        for i=1:length(t)
            t[i].re.value = real(v[i])
            t[i].im.value = imag(v[i])
            #(copyto!(value(t), v); nothing)
        end
    end
    return nothing
end
@inline function ReverseDiff.deriv!(t::Vector{<:Complex}, v::Vector{ComplexF64})
    for i=1:length(t)
        t[i].re.deriv = real(v[i])
        t[i].im.deriv = imag(v[i])
        #(copyto!(value(t), v); nothing)
    end
    return nothing
end


_get_origin(t::T) where T = T <: ReverseDiff.TrackedReal ? t.origin : T <: Complex ? t.re.origin : nothing
_get_index(t::T) where T = T <: ReverseDiff.TrackedReal ? t.index : T <: Complex ? t.re.index : nothing

ReverseDiff.pull_value!(t::Complex{<:ReverseDiff.TrackedReal}) = (ReverseDiff.hasorigin(t) && ReverseDiff.value!(t, ReverseDiff.value(_get_origin(t))[_get_index(t)]); nothing)
ReverseDiff.pull_deriv!(t::Complex{<:ReverseDiff.TrackedReal}) = (ReverseDiff.hasorigin(t) && ReverseDiff.deriv!(t, ReverseDiff.deriv(_get_origin(t))[_get_index(t)]); nothing)
ReverseDiff.push_deriv!(t::Complex{<:ReverseDiff.TrackedReal}) = (ReverseDiff.hasorigin(t) && (_get_origin(t).deriv[_get_index(t)] = ReverseDiff.deriv(t)); nothing)

#ReverseDiff.pull_value!(t::Complex{<:ReverseDiff.TrackedReal}) = (ReverseDiff.hasorigin(t) && ReverseDiff.value!(t, ReverseDiff.value(t.origin)[t.index]); nothing)
#ReverseDiff.pull_deriv!(t::Complex{<:ReverseDiff.TrackedReal}) = (ReverseDiff.hasorigin(t) && ReverseDiff.deriv!(t, ReverseDiff.deriv(t.origin)[t.index]); nothing)
#ReverseDiff.push_deriv!(t::Complex{<:ReverseDiff.TrackedReal}) = (ReverseDiff.hasorigin(t) && (t.origin.deriv[t.index] = ReverseDiff.deriv(t)); nothing)


# seed/unseed... I don't think these are ever called, though.
function ReverseDiff.seed!(t::Complex{<:ReverseDiff.TrackedReal},i)
    error(" ")
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
function ReverseDiff.seed!(t::Complex{<:ReverseDiff.TrackedReal})
    error(" ")
end
function ReverseDiff.seed!(t::AbstractArray{<:Complex{<:ReverseDiff.TrackedReal}})
    error(" ")
end
function ReverseDiff.seed!(t::AbstractArray{<:Complex{<:ReverseDiff.TrackedReal}},i)
    error(" ")
end

function ReverseDiff.unseed!(t::Complex{<:ReverseDiff.TrackedReal})
    t.re.deriv = zero(ReverseDiff.derivtype(t)); ReverseDiff.push_deriv!(t.re)
    t.im.deriv = zero(ReverseDiff.derivtype(t)); ReverseDiff.push_deriv!(t.im)
    return nothing
end

#ReverseDiff.capture(t::Complex) = ReverseDiff.istracked(t) ?  map!(ReverseDiff.capture, similar(t), t) : copy(t)
#ReverseDiff.capture(t::Complex) = t
#ReverseDiff.capture(t::Complex{<:ReverseDiff.TrackedReal}) = ifelse(ReverseDiff.hastape(t), t, ReverseDiff.value(t))
function Base.similar(z::Complex{TR}) where TR <: ReverseDiff.TrackedReal
    error(" ")
    zero(Complex{TR})
    #Complex{TR}(zero(Float64(z)),zero(Float64(z)))
end

# track/track! functions. These are pretty central to unpacking larger structures and propagating derivatives.

ReverseDiff.track(x::Complex, tp::ReverseDiff.InstructionTape = ReverseDiff.InstructionTape()) = ReverseDiff.track(x, typeof(x), tp)
ReverseDiff.track(x::AbstractArray{<:Complex}, tp::ReverseDiff.InstructionTape = ReverseDiff.InstructionTape()) = ReverseDiff.track(x, eltype(x), tp)
ReverseDiff.track(x::Complex, ::Type{Complex{D}}, tp::ReverseDiff.InstructionTape = ReverseDiff.InstructionTape()) where {D} = Complex{ReverseDiff.TrackedReal(real(x), zero(D), tp), ReverseDiff.TrackedReal(imag(x), zero(D), tp)}
ReverseDiff.track(x::Complex, ::Type{D}, tp::ReverseDiff.InstructionTape = ReverseDiff.InstructionTape()) where {D} = Complex{ReverseDiff.TrackedReal(real(x), zero(D), tp), ReverseDiff.TrackedReal(imag(x), zero(D), tp)}
ReverseDiff.track(x::AbstractArray{<:Complex}, ::Type{Complex{D}}, tp::ReverseDiff.InstructionTape = ReverseDiff.InstructionTape()) where {D} = ReverseDiff.TrackedArray(x, fill!(similar(x,Complex{D}), zero(Complex{D})), tp)
ReverseDiff.track(x::AbstractArray{<:Complex}, ::Type{D}, tp::ReverseDiff.InstructionTape = ReverseDiff.InstructionTape()) where {D} = ReverseDiff.TrackedArray(x, fill!(similar(x,D), zero(Complex{D})), tp)

#ReverseDiff.track!(t::ReverseDiff.TrackedArray{Complex{<:ReverseDiff.TrackedReal}}, x::AbstractArray) = (ReverseDiff.value!(t,x); ReverseDiff.unseed!(t); t)
#ReverseDiff.track!(t::Complex{<:ReverseDiff.TrackedReal}, x::Real) = (ReverseDiff.value!(t,x); ReverseDiff.unseed!(t); t)
#ReverseDiff.track!(t::Complex{<:ReverseDiff.TrackedReal}, x::Complex) = (ReverseDiff.value!(t,x); ReverseDiff.unseed!(t); t)
#function ReverseDiff.track!(t::AbstractArray{Complex{ReverseDiff.TrackedReal{D,D,Nothing}}}, x::AbstractArray, tp::ReverseDiff.InstructionTape) where D
#    for i in eachindex(t)
#        t[i] = track(x[i], D, tp)
#    end
#    return t
#end

function Base.convert(::Type{Complex{T1}}, t::T2) where {T1<:ReverseDiff.TrackedReal, T2<:ReverseDiff.TrackedReal}
    V1, D1, O1 = ReverseDiff.valtype(T1), ReverseDiff.derivtype(T1), ReverseDiff.origintype(T1)
    tp = ReverseDiff.tape(t)
    #out = Complex(ReverseDiff.TrackedReal{V1,D1,O1}(zero(V1), zero(D1), tp), ReverseDiff.TrackedReal{V1,D1,O1}(zero(V1), zero(D1), tp))
    #out = Complex(ReverseDiff.TrackedReal{V1,D1,O1}(ReverseDiff._convert(V1, ReverseDiff.value(t)), ReverseDiff._convert(D1, ReverseDiff.deriv(t)), tp), ReverseDiff.TrackedReal{V1,D1,O1}(zero(V1), zero(D1), tp))
    #out = ReverseDiff.TrackedReal{V1,D1,O1}(ReverseDiff._convert(V1, ReverseDiff.value(t)), ReverseDiff._convert(D1, ReverseDiff.deriv(t)), tp)
    #@show V1 D1
    imag_tracked_zero = ReverseDiff.TrackedReal{V1,D1,O1}(convert(V1,0.0),convert(D1,0.0),tp,t.index,t.origin)
    out = Complex(convert(T1,t),convert(T1,imag_tracked_zero))
    #out.re.index = 1
    #out.im.value = zero(V1)
    #out.im.deriv = zero(D1)
    ReverseDiff.record!(tp, ReverseDiff.SpecialInstruction, convert, t, out)
    #ReverseDiff.record!(tp, ReverseDiff.SpecialInstruction, convert, imag_tracked_zero, out.im)
    #@show t out.im real(t) imag(t) convert(T1,imag(t)) imag_tracked_zero
    #println("Converting $t to $T1 (real to complex)\nV1: $V1\tD1: $D1\tO1: $O1")
    #@show t out
    #println("\n")
    return out
end

function Base.convert(::Type{Complex{T1}},t::Complex{T2}) where {T1 <: ReverseDiff.TrackedReal, T2 <: ReverseDiff.TrackedReal}
    #V1,D1,O1 = ReverseDiff.valtype(T1), ReverseDiff.derivtype(T1), ReverseDiff.origintype(T1)
    tp = ReverseDiff.tape(t)
    out = Complex(ReverseDiff._convert(T1,t.re),ReverseDiff._convert(T1,t.im))
    #ReverseDiff.record!(tp, ReverseDiff.SpecialInstruction, convert, t, out)
    ReverseDiff.record!(tp, ReverseDiff.SpecialInstruction, convert, t.re, out.re)
    ReverseDiff.record!(tp, ReverseDiff.SpecialInstruction, convert, t.im, out.im)
    #println("Converting $t to $T1 (complex to complex)\nV1: $V1\tD1: $D1\tO1: $O1")
    #@show t out
    #println("\n")
    return out

end

function Base.convert(::Type{Complex{T1}},t::C) where {T1 <: ReverseDiff.TrackedReal, C <: Complex}
    error(" ")
end


function Base.convert(::Type{Complex{T1}},t::R) where {T1 <: ReverseDiff.TrackedReal, R <: Real}
    error(" ")
end

Base.convert(::Type{Complex{T}}, t::Complex{T}) where {T <: ReverseDiff.TrackedReal} = t
Base.promote_rule(::Type{Complex{ReverseDiff.TrackedReal{V1,D1,O1}}}, ::Type{Complex{ReverseDiff.TrackedReal{V2,D2,O2}}}) where {V1,V2,D1,D2,O1,O2} = Complex{ReverseDiff.TrackedReal{promote_type(V1,V2),promote_type(D1,D2),nothing},ReverseDiff.TrackedReal{promote_type(V1,V2),promote_type(D1,D2),nothing}}

Base.one(::Type{Complex{ReverseDiff.TrackedReal{V,D,O}}}) where {V,D,O} = error(" ")
Base.zero(::Type{Complex{ReverseDiff.TrackedReal{V,D,O}}}) where {V,D,O} = Complex(ReverseDiff.TrackedReal{V,D,O}(zero(V)),ReverseDiff.TrackedReal{V,D,O}(zero(V)))

function ReverseDiff.capture(t::Complex{<:ReverseDiff.TrackedReal})
    
    #@show t
    if ReverseDiff.hastape(t)
        return Complex(ReverseDiff.capture(t.re),ReverseDiff.capture(t.im))
    else
        return t
    end

end

# these next functions make sure derivatives get propagated. Unimplemented ones still return an error so that I know if one of them is needed.
function ReverseDiff.increment_deriv!(t::Complex,x)
    real(t).deriv += real(x)
    imag(t).deriv += imag(x)
    #error("t: $t\tx: $x")
end
function ReverseDiff.increment_deriv!(t::Complex,x,i)
    error("t: $t\tx: $x\ti: $i")
end

@inline function ReverseDiff.increment_deriv!(t::ReverseDiff.TrackedReal,x::Complex)

    t.deriv += real(x)
    #error("$t\n$x")

end

function ReverseDiff._add_to_deriv!(a::Complex,b)
    if ReverseDiff.istracked(a)
        real(a).deriv += real(b)
        imag(a).deriv += imag(b)
    end
    return nothing
end

ReverseDiff.increment_deriv!(a::AbstractArray{<:Complex{<:ReverseDiff.TrackedReal}},b::Real,i) = ReverseDiff.increment_deriv!(a[i],b)
ReverseDiff.increment_deriv!(a::AbstractArray{<:Complex{<:ReverseDiff.TrackedReal}},b::Complex,i) = ReverseDiff.increment_deriv!(a[i],b)

function ReverseDiff._add_to_deriv!(a::AbstractArray{<:Complex{<:ReverseDiff.TrackedReal}},b)

    if ReverseDiff.istracked(a)
        ReverseDiff.increment_deriv!(a,b)
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