#####
##### functions that should be overloaded for use in the FMM
#####

Base.getindex(sys, i, ::Position) = @error "getindex! not overloaded for `FLOWFMM.Position` for type $(typeof(sys)); cannot run FMM"
function Base.getindex(sys, i, ::Radius)
    @warn "getindex! not overloaded for `FLOWFMM.Radius` for type $(typeof(sys)); zero radius assumed"
    return 0.0
end
function Base.getindex(sys, i, ::ScalarPotential)
    @warn "getindex! not overloaded for `FLOWFMM.ScalarPotential` for type $(typeof(sys)); zero assumed"
    return 0.0
end
function Base.getindex(sys, i, ::Velocity)
    @warn "getindex! not overloaded for `FLOWFMM.Velocity` for type $(typeof(sys)); zero assumed"
    return zero(SVector{3,Float64})
end
function Base.getindex(sys, i, ::VelocityGradient)
    @warn "getindex! not overloaded for `FLOWFMM.VelocityGradient` for type $(typeof(sys)); zero assumed"
    return zero(SMatrix{3,3,Float64,9})
end
function Base.getindex(sys, i, ::ScalarStrength)
    @warn "getindex! not overloaded for `FLOWFMM.ScalarStrength` for type $(typeof(sys)); zero assumed"
    return 0.0
end
function Base.getindex(sys, i, ::VectorStrength)
    @warn "getindex! not overloaded for `FLOWFMM.VectorStrength` for type $(typeof(sys)); zero assumed"
    return zero(SVector{3,Float64})
end
Base.getindex(sys, i, ::Vertex, i_vertex) = @error "`getindex!` not overloaded for `FLOWFMM.Vertex` for type $(typeof(sys)); cannot run FMM"

Base.setindex!(sys, val, i) = @warn "setindex! not overloaded for type $(typeof(sys)); octree will have errors"
Base.setindex!(sys, val, i, ::ScalarPotential) = nothing
Base.setindex!(sys, val, i, ::VectorPotential) = nothing
Base.setindex!(sys, val, i, ::Velocity) = nothing
Base.setindex!(sys, val, i, ::VelocityGradient) = nothing

Base.length(sys) = @error "Base.length() not overloaded for type $(typeof(sys))"
# Base.eltype(::SystemType) = 
buffer_element(system) = @error "`buffer_element` not overloaded for `system::`$(typeof(system)); try using `SortWrapper(system)`"

function B2M!(system, branch, bodies_index, harmonics, expansion_order)
    @warn "B2M! function not overloaded for type $(typeof(system)); interaction ignored"
    return nothing
end

function direct!(target_system, target_index, source_system, source_index)
    @warn "direct! function not overloaded for type $(typeof(source_system)); interaction ignored"
    return nothing
end
