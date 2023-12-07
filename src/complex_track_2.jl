# Idea: I really want to re-use already-written interfaces. Ideally, I would hand all ReverseDiff operations to the TrackedReal/TrackedArray interface and I would hand all Complex operations
#    to the Complex interface. In C++ language, I want an object that inherits from both TrackedArray and Complex.

struct TrackedComplexArray{V,D,N,VA,DA} <: AbstractArray{Complex{ReverseDiff.TrackedReal{V,D,ReverseDiff.TrackedArray{V,D,N,VA,DA}},N}}

    value::VA
    deriv::DA
    tape::ReverseDiff.InstructionTape
    function TrackedComplexArray{V,D,N,VA,DA}(value::AbstractArray{V,N},
                                              deriv::AbstractArray{D,N},
                                              tape::ReverseDiff.InstructionTape) where {V<:Complex, D<:Complex, N, VA, DA}

        @assert IndexStyle(value) === IndexLinear()
        @assert size(value) === size(deriv)
        return new{V,D,N,VA,DA}(value, deriv, tape)

    end

end

function TrackedComplexArray(value::AbstractArray{V,N},
                             deriv::AbstractArray{D,N},
                             tape::ReverseDiff.InstructionTape) where {V<:Complex,D,N}
    return TrackedComplexArray{V,D,N,typeof(value),typeof(deriv)}(value, deriv, tape)
end

# getters

ReverseDiff.istracked(::TrackedComplexArray) = true
@inline ReverseDiff.value(t::TrackedComplexArray) = t.value
@inline ReverseDiff.deriv(t::TrackedComplexArray) = t.deriv

@inline ReverseDiff.valtype(::TrackedComplexArray{V}) where V = V
@inline ReverseDiff.valtype(::Type{TrackedComplexArray{V,D,N,VA,DA}}) where {V,D,N,VA,DA} = V
@inline ReverseDiff.derivtype(::TrackedComplexArray{V,D}) where {V,D} = D
@inline ReverseDiff.derivtype(::Type{TrackedComplexArray{V,D,N,VA,DA}}) where {V,D,N,VA,DA} = D

@inline ReverseDiff.hastape(t::TrackedComplexArray) = ReverseDiff.tape(t) !== ReverseDiff.NULL_TAPE
@inline ReverseDiff.tape(t::TrackedComplexArray) = t.tape

# setters

@inline ReverseDiff.value!(t::TrackedComplexArray, v::AbstractArray{<:Complex}) = (copyto!(ReverseDiff.value(t), v); nothing)
@inline ReverseDiff.value!(t::TrackedComplexArray, v::AbstractArray{<:Real}) = (map!(copyto!(_t.re, _v),ReverseDiff.value(t),v); nothing) # might need to write this one out as an actual loop

@inline ReverseDiff.deriv!(t::TrackedComplexArray, v::AbstractArray{<:Complex}) = (copyto!(ReverseDiff.deriv(t), v); nothing)
@inline ReverseDiff.deriv!(t::TrackedComplexArray, v::AbstractArray{<:Real}) = (map!(copyto!(_t.re, _v),ReverseDiff.deriv(t),v); nothing) # might need to write this one out as an actual loop

# pushing/pulling functions already do nothing with TrackedArrays (and the default is also nothing), so no additional methods are needed.

# seed/unseed

ReverseDiff.seed!(t::TrackedComplexArray, i) = (t.deriv[i] = one(ReverseDiff.derivtype(t)); nothing)
ReverseDiff.unseed!(t::TrackedComplexArray) = (fill!(ReverseDiff.deriv(t), zero(ReverseDiff.derivtype(t))); nothing)
ReverseDiff.unseed!(t::TrackedComplexArray, i) = (t.deriv[i] = zero(ReverseDiff.derivtype(t)); nothing)

# capture

ReverseDiff.capture(t::TrackedComplexArray) = t

# conversion/promotion

# there are a total of 4 cases for each _convert function, since either argument can be complex or real. Only the case of converting a complex number to a real one should return an error.

#=ReverseDiff._convert(::Type{C}, t::ReverseDiff.TrackedReal) where {C <: Complex} = C(ReverseDiff.value(t))
ReverseDiff._convert(::Type{C}, t::Complex{<:ReverseDiff.TrackedReal}) where {C <: Complex} = C(ReverseDiff.value(t.re), ReverseDiff.value(t.im))
ReverseDiff._convert(::Type{R}, t::Complex{<:ReverseDiff.TrackedReal}) where {R <: Real} = error("Cannot convert t ($t), which is a $(typeof(t)), to $R!")

ReverseDiff._convert(::Type{C}, t::Real) where {C <: Complex{<:ReverseDiff.TrackedReal}} = C(t)
ReverseDiff._convert(::Type{C}, t::Complex) where {C <: Complex{<:ReverseDiff.TrackedReal}} = C(t)
ReverseDiff._convert(::Type{R}, t::Complex) where {R <: ReverseDiff.TrackedReal} = error("Cannot convert t ($t), which is a $(typeof(t)), to $R!")

ReverseDiff._convert(::Type{T1}, t::T2) where {T1 <: Number, T2 <: Number} = T1(t) # specific cases can be handled by the existing Julia conversion rules.

ReverseDiff._convert(::Type{C}, t::ReverseDiff.TrackedReal) where {C <: Complex{<:ReverseDiff.TrackedReal}} = C(t)
ReverseDiff._convert(::Type{C}, t::Complex{<:ReverseDiff.TrackedReal}) where {C <: Complex{<:ReverseDiff.TrackedReal}} = C(t)
ReverseDiff._convert(::Type{R}, t::Complex{<:ReverseDiff.TrackedReal}) where {R <: ReverseDiff.TrackedReal} = error("Cannot convert t ($t), which is a $(typeof(t)), to $R!")=#

# not sure I actually need to write the detailed convert! method, since it takes purely real inputs. I might need to write one for Complex{TrackedReal}, but we'll see.

Base.convert(::Type{T}, t::T) where {T<:TrackedComplexArray} = t

# AbstractArray interface

