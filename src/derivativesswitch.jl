"""
    DerivativesSwitch(scalar_potential, velocity, velocity_gradient)

Constructs a tuple of [`DerivativesSwitch`](@ref) objects.

# Arguments

- `scalar_potential::Vector{Bool}`: a vector of `::Bool` indicating whether the scalar potential should be computed for each target system
- `velocity::Vector{Bool}`: a vector of `::Bool` indicating whether the velocity should be computed for each target system
- `velocity_gradient::Vector{Bool}`: a vector of `::Bool` indicating whether the velocity gradient should be computed for each target system

"""
function DerivativesSwitch(scalar_potential, velocity, velocity_gradient)
    return Tuple(DerivativesSwitch(ps,vps,vs,gs) for (ps,vps,vs,gs) in zip(scalar_potential, velocity, velocity_gradient))
end

"""
    DerivativesSwitch(scalar_potential, velocity, velocity_gradient)

Constructs a single [`DerivativesSwitch`](@ref) object.

# Arguments

- `scalar_potential::Bool`: a `::Bool` indicating whether the scalar potential should be computed for the target system
- `velocity::Bool`: a `::Bool` indicating whether the velocity should be computed for the target system
- `velocity_gradient::Bool`: a `::Bool` indicating whether the velocity gradient should be computed for the target system

"""
function DerivativesSwitch(scalar_potential::Bool, velocity::Bool, velocity_gradient::Bool)
    return DerivativesSwitch{scalar_potential, velocity, velocity_gradient}()
end

"""
    DerivativesSwitch(scalar_potential, velocity, velocity_gradient, target_systems)

Constructs a `::Tuple` of indentical [`DerivativesSwitch`](@ref) objects of the same length as `target_systems` (if it is a `::Tuple`), or a single [`DerivativesSwitch`](@ref) (if `target_system` is not a `::Tuple`)

# Arguments

- `scalar_potential::Bool`: a `::Bool` indicating whether the scalar potential should be computed for each target system
- `velocity::Bool`: a `::Bool` indicating whether the velocity should be computed for each target system
- `velocity_gradient::Bool`: a `::Bool` indicating whether the velocity gradient should be computed for each target system

"""
function DerivativesSwitch(scalar_potential::Bool, velocity::Bool, velocity_gradient::Bool, target_systems::Tuple)
    return Tuple(DerivativesSwitch{scalar_potential, velocity, velocity_gradient}() for _ in target_systems)
end

function DerivativesSwitch(scalar_potential, velocity, velocity_gradient, target_systems::Tuple)
    @assert length(scalar_potential) == length(velocity) == length(velocity_gradient) == length(target_systems) "length of inputs to DerivativesSwitch inconsistent"
    return Tuple(DerivativesSwitch{scalar_potential[i], velocity[i], velocity_gradient[i]}() for i in eachindex(target_systems))
end

function DerivativesSwitch(scalar_potential::Bool, velocity::Bool, velocity_gradient::Bool, target_system)
    return DerivativesSwitch{scalar_potential, velocity, velocity_gradient}()
end

DerivativesSwitch() = DerivativesSwitch{true, true, true}()

