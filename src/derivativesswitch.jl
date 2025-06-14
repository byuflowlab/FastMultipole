"""
    DerivativesSwitch(scalar_potential, vector_field, vector_gradient)

Constructs a tuple of [`DerivativesSwitch`](@ref) objects.

# Arguments

- `scalar_potential::Vector{Bool}`: a vector of `::Bool` indicating whether the scalar potential should be computed for each target system
- `vector_field::Vector{Bool}`: a vector of `::Bool` indicating whether the vector field should be computed for each target system
- `vector_gradient::Vector{Bool}`: a vector of `::Bool` indicating whether the vector gradient should be computed for each target system

"""
function DerivativesSwitch(scalar_potential, vector_field, vector_gradient)
    return Tuple(DerivativesSwitch(ps,vs,gs) for (ps,vs,gs) in zip(scalar_potential, vector_field, vector_gradient))
end

"""
    DerivativesSwitch(scalar_potential, vector_field, vector_gradient)

Constructs a single [`DerivativesSwitch`](@ref) object.

# Arguments

- `scalar_potential::Bool`: a `::Bool` indicating whether the scalar potential should be computed for the target system
- `vector_field::Bool`: a `::Bool` indicating whether the vector field should be computed for the target system
- `vector_gradient::Bool`: a `::Bool` indicating whether the vector gradient should be computed for the target system

"""
function DerivativesSwitch(scalar_potential::Bool, vector_field::Bool, vector_gradient::Bool)
    return DerivativesSwitch{scalar_potential, vector_field, vector_gradient}()
end

"""
    DerivativesSwitch(scalar_potential, vector_field, vector_gradient, target_systems)

Constructs a `::Tuple` of indentical [`DerivativesSwitch`](@ref) objects of the same length as `target_systems` (if it is a `::Tuple`), or a single [`DerivativesSwitch`](@ref) (if `target_system` is not a `::Tuple`)

# Arguments

- `scalar_potential::Bool`: a `::Bool` indicating whether the scalar potential should be computed for each target system
- `vector_field::Bool`: a `::Bool` indicating whether the vector field should be computed for each target system
- `vector_gradient::Bool`: a `::Bool` indicating whether the vector gradient should be computed for each target system

"""
function DerivativesSwitch(scalar_potential::Bool, vector_field::Bool, vector_gradient::Bool, target_systems::Tuple)
    return Tuple(DerivativesSwitch{scalar_potential, vector_field, vector_gradient}() for _ in target_systems)
end

function DerivativesSwitch(scalar_potential, vector_field, vector_gradient, target_systems::Tuple)
    @assert length(scalar_potential) == length(vector_field) == length(vector_gradient) == length(target_systems) "length of inputs to DerivativesSwitch inconsistent"
    return Tuple(DerivativesSwitch{scalar_potential[i], vector_field[i], vector_gradient[i]}() for i in eachindex(target_systems))
end

function DerivativesSwitch(scalar_potential::Bool, vector_field::Bool, vector_gradient::Bool, target_system)
    return DerivativesSwitch{scalar_potential, vector_field, vector_gradient}()
end

DerivativesSwitch() = DerivativesSwitch{true, true, true}()

