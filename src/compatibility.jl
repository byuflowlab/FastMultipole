#####
##### functions that should be overloaded for use in the FMM
#####

Base.getindex(sys, i, ::Body) = @error "getindex! not overloaded for `FastMultipole.Body` for type $(typeof(sys))"

Base.getindex(sys::AbstractArray, i, ::Body) = sys[i]

Base.getindex(sys, i, ::Position) = @error "getindex! not overloaded for `FastMultipole.Position` for type $(typeof(sys)); cannot run FMM"

const WARNING_FLAG_RADIUS = Array{Bool,0}(undef)
WARNING_FLAG_RADIUS[] = true

function Base.getindex(sys, i, ::Radius)
    if WARNING_FLAG_RADIUS[]
        @warn "getindex! not overloaded for `FastMultipole.Radius` for type $(typeof(sys)); zero radius assumed"
        WARNING_FLAG_RADIUS[] = false
    end
    return 0.0
end

const WARNING_FLAG_SCALAR_POTENTIAL = Array{Bool,0}(undef)
WARNING_FLAG_SCALAR_POTENTIAL[] = true

function Base.getindex(sys, i, ::ScalarPotential)
    if WARNING_FLAG_SCALAR_POTENTIAL[]
        @warn "getindex! not overloaded for `FastMultipole.ScalarPotential` for type $(typeof(sys)); zero assumed"
        WARNING_FLAG_SCALAR_POTENTIAL[] = false
    end
    return 0.0
end

const WARNING_FLAG_VECTOR_POTENTIAL = Array{Bool,0}(undef)
WARNING_FLAG_VECTOR_POTENTIAL[] = true

function Base.getindex(sys, i, ::VectorPotential)
    if WARNING_FLAG_VECTOR_POTENTIAL[]
        @warn "getindex! not overloaded for `FastMultipole.VectorPotential` for type $(typeof(sys)); zero assumed"
        WARNING_FLAG_VECTOR_POTENTIAL[] = false
    end
    return SVector{3}(0.0,0.0,0.0)
end

const WARNING_FLAG_VELOCITY = Array{Bool,0}(undef)
WARNING_FLAG_VELOCITY[] = true

function Base.getindex(sys, i, ::Velocity)
    if WARNING_FLAG_VELOCITY[]
        @warn "getindex! not overloaded for `FastMultipole.Velocity` for type $(typeof(sys)); zero assumed"
        WARNING_FLAG_VELOCITY[] = false
    end
    return zero(SVector{3,Float64})
end

const WARNING_FLAG_VELOCITY_GRADIENT = Array{Bool,0}(undef)
WARNING_FLAG_VELOCITY_GRADIENT[] = true

function Base.getindex(sys, i, ::VelocityGradient)
    if WARNING_FLAG_VELOCITY_GRADIENT[]
        @warn "getindex! not overloaded for `FastMultipole.VelocityGradient` for type $(typeof(sys)); zero assumed"
        WARNING_FLAG_VELOCITY_GRADIENT[] = false
    end
    return zero(SMatrix{3,3,Float64,9})
end

const WARNING_FLAG_STRENGTH = Array{Bool,0}(undef)
WARNING_FLAG_STRENGTH[] = true

function Base.getindex(sys, i, ::Strength)
    if WARNING_FLAG_STRENGTH[]
        @warn "getindex! not overloaded for `FastMultipole.Strength` for type $(typeof(sys)); zero assumed"
        WARNING_FLAG_STRENGTH[] = false
    end
    return 0.0
end

function Base.setindex(sys, i, val, ::Strength)
    if WARNING_FLAG_STRENGTH[]
        @warn "setindex! not overloaded for `FastMultipole.Strength` for type $(typeof(sys)); zero assumed"
        WARNING_FLAG_STRENGTH[] = false
    end
end

function Base.getindex(sys, i, ::Normal)
    @error "getindex! not overloaded for `FastMultipole.Normal` for type $(typeof(sys))"
end

Base.getindex(sys, i, ::Vertex, i_vertex) = @error "`getindex!` not overloaded for `FastMultipole.Vertex` for type $(typeof(sys)); cannot run FMM"

Base.setindex!(sys, val, i, ::Body) = @error "setindex! not overloaded for `FastMultipole.Body` for type $(typeof(sys))"

Base.setindex!(sys::AbstractArray, val, i, ::Body) = sys[i] = val

Base.setindex!(sys, val, i, ::ScalarPotential) = nothing

Base.setindex!(sys, val, i, ::VectorPotential) = nothing

Base.setindex!(sys, val, i, ::Velocity) = nothing

Base.setindex!(sys, val, i, ::VelocityGradient) = nothing

get_n_bodies(sys::AbstractArray) = length(sys)

get_n_bodies(sys) = @error "FastMultipole.get_n_bodies() not overloaded for type $(typeof(sys))"

buffer_element(system) = @error "`buffer_element` not overloaded for `system::`$(typeof(system)); try using `SortWrapper(system)` or overload `buffer_element`"

function B2M!(system, branch, bodies_index, harmonics, expansion_order)
    @warn "B2M! function not overloaded for type $(typeof(system)); interaction ignored"
    return nothing
end

function direct!(target_system, target_index, derivatives_switch, source_system, source_index)
    @warn "direct! function not overloaded for type $(typeof(source_system)); interaction ignored"
    return nothing
end

function direct_gpu!(target_system, target_indices, derivatives_switch, source_system, source_indices)
    @warn "direct_gpu! function not overloaded for type $(typeof(source_system)); interaction ignored"
    return nothing
end
