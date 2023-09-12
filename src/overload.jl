#####
##### functions that should be overloaded for use in the FMM
#####

# Base.getindex(sys, i, ::Position) =
# Base.getindex(sys, i, ::Radius) =
# Base.getindex(sys, i, ::Potential) =
# Base.getindex(sys, i, ::ScalarPotential) =
# Base.getindex(sys, i, ::Velocity) =
# Base.getindex(sys, i, ::VelocityGradient) =
# Base.getindex(sys, i) =
function Base.setindex!(sys, val, i)
    @warn "setindex! not overloaded for type $(typeof(sys)); octree will have errors"
    return nothing
end
function Base.setindex!(sys, val, i, ::ScalarPotential)
    # @warn "setindex! (ScalarPotential) not overloaded for type $(typeof(sys))"
    return nothing
end
function Base.setindex!(sys, val, i, ::Potential)
    # @warn "setindex! (Potential) not overloaded for type $(typeof(sys))"
    return nothing
end
function Base.setindex!(sys, val, i, ::Velocity)
    # @warn "setindex! (Velocity) not overloaded for type $(typeof(sys))"
    return nothing
end
function Base.setindex!(sys, val, i, ::VelocityGradient)
    # @warn "setindex! (VelocityGradient) not overloaded for type $(typeof(sys))"
    return nothing
end
# Base.length(system) = 
# Base.eltype(::SystemType) = 

function B2M!(branch, system, bodies_index, harmonics, expansion_order)
    @warn "B2M! function not overloaded for type $(typeof(system)); interaction ignored"
    return nothing
end

function direct!(target_system, target_index, source_system, source_index)
    @warn "direct! function not overloaded for type $(typeof(source_system)); interaction ignored"
    return nothing
end