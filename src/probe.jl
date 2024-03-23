function ProbeSystem(n_probes::Int, TF=Float64; scalar_potential=false, vector_potential=false, velocity=false, velocity_gradient=false)
    positions = Vector{SVector{3,TF}}(undef,n_probes)
    return ProbeSystem(positions; scalar_potential, vector_potential, velocity, velocity_gradient)
end

function ProbeSystem(positions::Vector{SVector{3,TF}}; scalar_potential=false, vector_potential=false, velocity=false, velocity_gradient=false) where TF
    n_probes = length(positions)
    scalar_potentials = scalar_potential ? zeros(TF,n_probes) : nothing
    vector_potentials = vector_potential ? zeros(SVector{3,TF},n_probes) : nothing
    velocities = velocity ? zeros(SVector{3,TF},n_probes) : nothing
    velocity_gradients = velocity_gradient ? zeros(SMatrix{3,3,TF,9},n_probes) : nothing
    return ProbeSystem(positions, scalar_potentials, vector_potentials, velocities, velocity_gradients)
end

function ProbeSystem(positions::Matrix{TF}; scalar_potential=false, vector_potential=false, velocity=false, velocity_gradient=false) where TF
    n_probes = size(positions,2)
    positions = [SVector{3,TF}(positions[1,i], positions[2,i], positions[3,i]) for i in 1:n_probes]
    return ProbeSystem(positions; scalar_potential, vector_potential, velocity, velocity_gradient)
end

# fast, flexible access that avoids excess storage
@inline get_scalar_potential(probe_system::ProbeSystem{TF,Nothing,<:Any,<:Any,<:Any}, i) where TF = zero(TF)
@inline get_scalar_potential(probe_system::ProbeSystem, i) = probe_system.scalar_potential[i]
@inline get_vector_potential(probe_system::ProbeSystem{TF,<:Any,Nothing,<:Any,<:Any}, i) where TF = zero(SVector{3,TF})
@inline get_vector_potential(probe_system::ProbeSystem, i) = probe_system.vector_potential[i]
@inline get_velocity(probe_system::ProbeSystem{TF,<:Any,<:Any,Nothing,<:Any}, i) where TF = zero(SVector{3,TF})
@inline get_velocity(probe_system::ProbeSystem, i) = probe_system.velocity[i]
@inline get_velocity_gradient(probe_system::ProbeSystem{TF,<:Any,<:Any,<:Any,Nothing}, i) where TF = zero(SMatrix{3,3,TF,9})
@inline get_velocity_gradient(probe_system::ProbeSystem, i) = probe_system.velocity_gradient[i]
@inline set_scalar_potential!(probe_system::ProbeSystem{<:Any,Nothing,<:Any,<:Any,<:Any}, val, i) = nothing
@inline set_scalar_potential!(probe_system::ProbeSystem, val, i) = probe_system.scalar_potential[i] = val
@inline set_vector_potential!(probe_system::ProbeSystem{<:Any,<:Any,Nothing,<:Any,<:Any}, val, i) = nothing
@inline set_vector_potential!(probe_system::ProbeSystem, val, i) = probe_system.vector_potential[i] = val
@inline set_velocity!(probe_system::ProbeSystem{<:Any,<:Any,<:Any,Nothing,<:Any}, val, i) = nothing
@inline set_velocity!(probe_system::ProbeSystem, val, i) = probe_system.velocity[i] = val
@inline set_velocity_gradient!(probe_system::ProbeSystem{<:Any,<:Any,<:Any,<:Any,Nothing}, val, i) = nothing
@inline set_velocity_gradient!(probe_system::ProbeSystem, val, i) = probe_system.velocity_gradient[i] = val

# compatibility access functions
Base.getindex(probe_system::ProbeSystem, i, ::Position) = probe_system.position[i]
Base.getindex(probe_system::ProbeSystem{TF,<:Any,<:Any,<:Any,<:Any}, i, ::Radius) where TF = zero(TF)
Base.getindex(probe_system::ProbeSystem, i, ::ScalarPotential) = get_scalar_potential(probe_system, i)
Base.getindex(probe_system::ProbeSystem, i, ::VectorPotential) = get_vector_potential(probe_system, i)
Base.getindex(probe_system::ProbeSystem, i, ::Velocity) = get_velocity(probe_system, i)
Base.getindex(probe_system::ProbeSystem, i, ::VelocityGradient) = get_velocity_gradient(probe_system, i)
Base.getindex(probe_system::ProbeSystem{TF,<:Any,<:Any,<:Any,<:Any}, i, ::ScalarStrength) where TF = zero(TF)
Base.getindex(probe_system::ProbeSystem, i, ::Body) = probe_system.position[i], get_scalar_potential(probe_system, i), get_vector_potential(probe_system, i), get_velocity(probe_system, i), get_velocity_gradient(probe_system, i)

function Base.setindex!(probe_system::ProbeSystem, val, i, ::Body)
    (position, scalar_potential, vector_potential, velocity, velocity_gradient) = val
    probe_system.position[i] = position
    set_scalar_potential!(probe_system, scalar_potential, i)
    set_vector_potential!(probe_system, vector_potential, i)
    set_velocity!(probe_system, velocity, i)
    set_velocity_gradient!(probe_system, velocity_gradient, i)

    return nothing
end
function Base.setindex!(probe_system::ProbeSystem, val, i, ::ScalarPotential)
    set_scalar_potential!(probe_system, val, i)
end
function Base.setindex!(probe_system::ProbeSystem, val, i, ::VectorPotential)
    set_vector_potential!(probe_system, val, i)
end
function Base.setindex!(probe_system::ProbeSystem, val, i, ::Velocity)
    set_velocity!(probe_system, val, i)
end
function Base.setindex!(probe_system::ProbeSystem, val, i, ::VelocityGradient)
    set_velocity_gradient!(probe_system, val, i)
end
get_n_bodies(probe_system::ProbeSystem) = length(probe_system.position)
Base.eltype(::ProbeSystem{TF,<:Any,<:Any,<:Any,<:Any}) where TF = TF

buffer_element(probe_system::ProbeSystem) = probe_system.position[1], get_scalar_potential(probe_system,1), get_vector_potential(probe_system,1), get_velocity(probe_system,1), get_velocity_gradient(probe_system,1)

B2M!(system::ProbeSystem, args...) = nothing

direct!(target_system, target_index, source_system::ProbeSystem, source_index) = nothing

function reset!(values::Vector{TV}) where TV
    for i in eachindex(values)
        values[i] = zero(TV)
    end
end

reset!(::Nothing) = nothing

function reset!(probes::ProbeSystem)
    reset!(probes.scalar_potential)
    reset!(probes.vector_potential)
    reset!(probes.velocity)
    reset!(probes.velocity_gradient)
end