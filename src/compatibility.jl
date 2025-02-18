#------- functions that should be overloaded for use in the FMM -------#

Base.getindex(sys, i, ::Body) = @error "getindex! not overloaded for `FastMultipole.Body` for type $(typeof(sys))"

Base.getindex(sys::AbstractArray, i, ::Body) = sys[i]

Base.getindex(sys, i, ::Position) = @error "getindex! not overloaded for `FastMultipole.Position` for type $(typeof(sys)); cannot run FMM"

function Base.getindex(sys, i, ::Radius)
    if WARNING_FLAG_RADIUS[]
        @warn "getindex! not overloaded for `FastMultipole.Radius` for type $(typeof(sys)); zero radius assumed"
        WARNING_FLAG_RADIUS[] = false
    end
    return 0.0
end

function Base.getindex(sys, i, ::ScalarPotential)
    if WARNING_FLAG_SCALAR_POTENTIAL[]
        @warn "getindex! not overloaded for `FastMultipole.ScalarPotential` for type $(typeof(sys)); zero assumed"
        WARNING_FLAG_SCALAR_POTENTIAL[] = false
    end
    return 0.0
end

#function Base.getindex(sys, i, ::VectorPotential)
#    if WARNING_FLAG_VECTOR_POTENTIAL[]
#        @warn "getindex! not overloaded for `FastMultipole.VectorPotential` for type $(typeof(sys)); zero assumed"
#        WARNING_FLAG_VECTOR_POTENTIAL[] = false
#    end
#    return SVector{3}(0.0,0.0,0.0)
#end

function Base.getindex(sys, i, ::Velocity)
    if WARNING_FLAG_VELOCITY[]
        @warn "getindex! not overloaded for `FastMultipole.Velocity` for type $(typeof(sys)); zero assumed"
        WARNING_FLAG_VELOCITY[] = false
    end
    return zero(SVector{3,Float64})
end

function Base.getindex(sys, i, ::VelocityGradient)
    if WARNING_FLAG_VELOCITY_GRADIENT[]
        @warn "getindex! not overloaded for `FastMultipole.VelocityGradient` for type $(typeof(sys)); zero assumed"
        WARNING_FLAG_VELOCITY_GRADIENT[] = false
    end
    return zero(SMatrix{3,3,Float64,9})
end

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

#Base.setindex!(sys, val, i, ::VectorPotential) = nothing

Base.setindex!(sys, val, i, ::Velocity) = nothing

Base.setindex!(sys, val, i, ::VelocityGradient) = nothing

get_n_bodies(sys::AbstractArray) = length(sys)

get_n_bodies(sys) = @error "FastMultipole.get_n_bodies() not overloaded for type $(typeof(sys))"

buffer_element(system) = @error "`buffer_element` not overloaded for `system::`$(typeof(system)); try using `SortWrapper(system)` or overload `buffer_element`"

function body_to_multipole!(system, branch, bodies_index, harmonics, expansion_order)
    if WARNING_FLAG_B2M[]
        @warn "body_to_multipole! not overloaded for type $(typeof(system)); interaction ignored"
        WARNING_FLAG_B2M[] = false
    end
    return nothing
end

function direct!(target_system, target_index, derivatives_switch, source_system, source_index)
    if WARNING_FLAG_DIRECT[]
        @warn "direct! not overloaded for type $(typeof(source_system)); interaction ignored"
        WARNING_FLAG_DIRECT[] = false
    end
    return nothing
end

#------- access functions for use with a matrix of targets used as input to direct! -------#

const i_POSITION = 1:3
const i_SCALAR_POTENTIAL = 4
const i_VELOCITY = 5:7
const i_VELOCITY_GRADIENT = 8:16

#--- getters ---#

Base.getindex(sys::Matrix{TF}, i, ::Position) where TF = SVector{3,TF}(sys[i_POSITION[1], i], sys[i_POSITION[2], i], sys[i_POSITION[3], i])
Base.getindex(sys::Matrix{TF}, i, ::ScalarPotential) where TF = sys[i_SCALAR_POTENTIAL, i]
Base.getindex(sys::Matrix{TF}, i, ::Velocity) where TF =
    SVector{3,TF}(sys[i_VELOCITY[1], i], sys[i_VELOCITY[2], i], sys[i_VELOCITY[3], i])
Base.getindex(sys::Matrix{TF}, i, ::VelocityGradient) where TF =
    SMatrix{3,3,TF,9}(sys[i_VELOCITY_GRADIENT[1], i], sys[i_VELOCITY_GRADIENT[2], i], sys[i_VELOCITY_GRADIENT[3], i],
    sys[i_VELOCITY_GRADIENT[4], i], sys[i_VELOCITY_GRADIENT[5], i], sys[i_VELOCITY_GRADIENT[6], i],
    sys[i_VELOCITY_GRADIENT[7], i], sys[i_VELOCITY_GRADIENT[8], i], sys[i_VELOCITY_GRADIENT[9], i])

get_n_bodies(sys::Matrix) = size(sys, 2)

#--- setters ---#

Base.setindex!(sys::Matrix, val, i, ::Position) = @inbounds sys[i_POSITION, i] .= val
Base.setindex!(sys::Matrix, val, i, ::ScalarPotential) = sys[i_SCALAR_POTENTIAL, i] = val
Base.setindex!(sys::Matrix, val, i, ::Velocity) = @inbounds sys[i_VELOCITY, i] .= val

function Base.setindex!(sys::Matrix, val, i, ::VelocityGradient)
    for (jj,j) in enumerate(i_VELOCITY_GRADIENT)
        sys[j, i] = val[jj]
    end
end

#--- auxilliary functions ---#

function reset!(systems::Tuple)
    for system in systems
        reset!(system)
    end
end

function reset!(system)
    z_potential = zero(eltype(system))
    z_velocity = zero(SVector{3,eltype(system)})
    z_velocity_gradient = zero(SMatrix{3,3,eltype(system),9})
    for i in 1:get_n_bodies(system)
        system[i, ScalarPotential()] = z_potential
        system[i, Velocity()] = z_velocity
        system[i, VelocityGradient()] = z_velocity_gradient
    end
end
