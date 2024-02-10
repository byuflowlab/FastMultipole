#####
##### functions that should be overloaded for use in the FMM
#####

Base.getindex(sys, i, ::Position) = @error "getindex! not overloaded for `FLOWFMM.Position` for type $(typeof(sys)); cannot run FMM"

const WARNING_FLAG_RADIUS = Array{Bool,0}(undef)
WARNING_FLAG_RADIUS[] = true

function Base.getindex(sys, i, ::Radius)
    if WARNING_FLAG_RADIUS[] 
        @warn "getindex! not overloaded for `FLOWFMM.Radius` for type $(typeof(sys)); zero radius assumed"
        WARNING_FLAG_RADIUS[] = false
    end
    return 0.0
end

const WARNING_FLAG_SCALAR_POTENTIAL = Array{Bool,0}(undef)
WARNING_FLAG_SCALAR_POTENTIAL[] = true

function Base.getindex(sys, i, ::ScalarPotential)
    if WARNING_FLAG_SCALAR_POTENTIAL[] 
        @warn "getindex! not overloaded for `FLOWFMM.ScalarPotential` for type $(typeof(sys)); zero assumed"
        WARNING_FLAG_SCALAR_POTENTIAL[] = false
    end
    return 0.0
end

const WARNING_FLAG_VECTOR_POTENTIAL = Array{Bool,0}(undef)
WARNING_FLAG_VECTOR_POTENTIAL[] = true

function Base.getindex(sys, i, ::VectorPotential)
    if WARNING_FLAG_VECTOR_POTENTIAL[] 
        @warn "getindex! not overloaded for `FLOWFMM.VectorPotential` for type $(typeof(sys)); zero assumed"
        WARNING_FLAG_VECTOR_POTENTIAL[] = false
    end
    return SVector{3}(0.0,0.0,0.0)
end

const WARNING_FLAG_VELOCITY = Array{Bool,0}(undef)
WARNING_FLAG_VELOCITY[] = true

function Base.getindex(sys, i, ::Velocity)
    if WARNING_FLAG_VELOCITY[] 
        @warn "getindex! not overloaded for `FLOWFMM.Velocity` for type $(typeof(sys)); zero assumed"
        WARNING_FLAG_VELOCITY[] = false
    end
    return zero(SVector{3,Float64})
end

const WARNING_FLAG_VELOCITY_GRADIENT = Array{Bool,0}(undef)
WARNING_FLAG_VELOCITY_GRADIENT[] = true

function Base.getindex(sys, i, ::VelocityGradient)
    if WARNING_FLAG_VELOCITY_GRADIENT[] 
        @warn "getindex! not overloaded for `FLOWFMM.VelocityGradient` for type $(typeof(sys)); zero assumed"
        WARNING_FLAG_VELOCITY_GRADIENT[] = false
    end
    return zero(SMatrix{3,3,Float64,9})
end

const WARNING_FLAG_SCALAR_STRENGTH = Array{Bool,0}(undef)
WARNING_FLAG_SCALAR_STRENGTH[] = true

function Base.getindex(sys, i, ::ScalarStrength)
    if WARNING_FLAG_SCALAR_STRENGTH[] 
        @warn "getindex! not overloaded for `FLOWFMM.ScalarStrength` for type $(typeof(sys)); zero assumed"
        WARNING_FLAG_SCALAR_STRENGTH[] = false
    end
    return 0.0
end

const WARNING_FLAG_VECTOR_STRENGTH = Array{Bool,0}(undef)
WARNING_FLAG_VECTOR_STRENGTH[] = true

function Base.getindex(sys, i, ::VectorStrength)
    if WARNING_FLAG_VECTOR_STRENGTH[] 
        @warn "getindex! not overloaded for `FLOWFMM.VectorStrength` for type $(typeof(sys)); zero assumed"
        WARNING_FLAG_VECTOR_STRENGTH[] = false
    end
    return zero(SVector{3,Float64})
end

function Base.getindex(sys, i, ::Normal)
    @error "getindex! not overloaded for `FLOWFMM.Normal` for type $(typeof(sys))"
end

Base.getindex(sys, i, ::Vertex, i_vertex) = @error "`getindex!` not overloaded for `FLOWFMM.Vertex` for type $(typeof(sys)); cannot run FMM"

Base.setindex!(sys, val, i) = @error "setindex! not overloaded for type $(typeof(sys))"

Base.setindex!(sys, val, i, ::ScalarPotential) = nothing

Base.setindex!(sys, val, i, ::VectorPotential) = nothing

Base.setindex!(sys, val, i, ::Velocity) = nothing

Base.setindex!(sys, val, i, ::VelocityGradient) = nothing

Base.length(sys) = @error "Base.length() not overloaded for type $(typeof(sys))"

buffer_element(system) = @error "`buffer_element` not overloaded for `system::`$(typeof(system)); try using `SortWrapper(system)` or overload `buffer_element`"

function B2M!(system, branch, bodies_index, harmonics, expansion_order)
    @warn "B2M! function not overloaded for type $(typeof(system)); interaction ignored"
    return nothing
end

function direct!(target_system, target_index, source_system, source_index)
    @warn "direct! function not overloaded for type $(typeof(source_system)); interaction ignored"
    return nothing
end
