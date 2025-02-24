#------- functions that should be overloaded for each user-defined system for use in the FMM -------#

#--- buffer functions ---#

"""
    source_system_to_buffer!(buffer::Matrix, i_buffer, system::{UserDefinedSystem}, i_body)

Compatibility function used to sort source systems. It should be overloaded for each system (where `{UserDefinedSystem}` is replaced with the type of the user-defined system) to be used as a source and should behave as follows. For the `i_body`th body contained inside of `system`,

* `buffer[1:3, i_buffer]` should be set to the x, y, and z coordinates of the body position used for sorting into the octree
* `buffer[4, i_buffer]` should be set to the radius beyond which a multipole expansion is allowed to be evaluated (e.g. for panels, or other bodies of finite area/volume)
* `buffer[5:4+strength_dims, i_buffer]` should be set to the body strength, which is a vector of length `strength_dims`

Any additional information required for either forming multipole expansions or computing direct interactions should be stored in the rest of the column.

If a body contains vertices that are required for, e.g. computing multipole coefficients of dipole panels, these must be stored immediately following the body strength, and should be listed in a counter-clockwise order. For example, if I am using vortex tri-panels with `strength_dims=3`, I would set `buffer[8:10,i] .= v1`, `buffer[11:13,i] .= v2`, and `bufer[14:16,i] .= v3`, where `v1`, `v2`, and `v3` are listed according to the right-hand-rule with thumb aligned with the panel normal vector.

Note that any system acting only as a target need not overload `source_system_to_buffer!`.

"""
function source_system_to_buffer!(buffer, i_buffer, system, i_body)
    throw("source_system_to_buffer! not overloaded for type $(typeof(system))")
end

"""
    data_per_body(system::{UserDefinedSystem})

Returns the number of values used to represent a single body in a source system. Should be overloaded for each user-defined system object (where `{UserDefinedSystem}` is replaced with the type of the user-defined system).

"""
function data_per_body(system)
    throw("data_per_body not overloaded for type $(typeof(system))")
end

#--- getters ---#

"""
    get_position(system::{UserDefinedSystem}, i)

Returns a (static) vector of length 3 containing the x, y, and z coordinates of the position of the `i`th body. Should be overloaded for each user-defined system object (where `{UserDefinedSystem}` is replaced with the type of the user-defined system).

"""
function get_position(system, i)
    throw("get_position not overloaded for type $(typeof(system))")
end

"""
    strength_dims(system::{UserDefinedSystem})

Returns the cardinality of the vector used to define the strength of each body inside `system`. E.g., a point mass would return 1, and a point dipole would return 3. Should be overloaded for each user-defined system object (where `{UserDefinedSystem}` is replaced with the type of the user-defined system).

"""
function strength_dims(system)
    throw("strength_dims() not overloaded for type $(typeof(system))")
end

"""
    get_normal(source_buffer, source_system::{UserDefinedSystem}, i)

**OPTIONAL OVERLOAD:**

Returns the unit normal vector for the `i`th body of `source_buffer`. May be (optionally) overloaded for a user-defined system object (where `{UserDefinedSystem}` is replaced with the type of the user-defined system); otherwise, the default behavior assumes counter-clockwise ordered vertices. Note that whatever method is used should match `source_system_to_buffer!` for each system.

"""
function get_normal(source_buffer, source_system, i_body)
    v1 = get_vertex(source_buffer, source_system, i_body, 1)
    v2 = get_vertex(source_buffer, source_system, i_body, 2)
    v3 = get_vertex(source_buffer, source_system, i_body, 3)
    normal = cross(v2-v1, v3-v1)

    return normal / norm(normal)
end

"""
    get_n_bodies(system::{UserDefinedSystem})

Returns the number of bodies contained inside `system`. Should be overloaded for each user-defined system object (where `{UserDefinedSystem}` is replaced with the type of the user-defined system).
"""
get_n_bodies(system) = throw("FastMultipole.get_n_bodies() not overloaded for type $(typeof(system))")

# function get_scalar_potential(system, i)
#     if WARNING_FLAG_SCALAR_POTENTIAL[]
#         @warn "get_scalar_potential not overloaded for type $(typeof(system)); zero assumed"
#         WARNING_FLAG_SCALAR_POTENTIAL[] = false
#     end
#     return zero(eltype(system))
# end

#function Base.getindex(sys, i, ::VectorPotential)
#    if WARNING_FLAG_VECTOR_POTENTIAL[]
#        @warn "getindex! not overloaded for `FastMultipole.VectorPotential` for type $(typeof(sys)); zero assumed"
#        WARNING_FLAG_VECTOR_POTENTIAL[] = false
#    end
#    return SVector{3}(0.0,0.0,0.0)
#end

# function get_velocity(system, i)
#     if WARNING_FLAG_VELOCITY[]
#         @warn "get_velocity not overloaded for type $(typeof(system)); zero assumed"
#         WARNING_FLAG_VELOCITY[] = false
#     end
#     return zero(SVector{3,eltype(system)})
# end

# function get_velocity_gradient(system, i)
#     if WARNING_FLAG_VELOCITY_GRADIENT[]
#         @warn "get_velocity_gradient not overloaded for type $(typeof(system)); zero assumed"
#         WARNING_FLAG_VELOCITY_GRADIENT[] = false
#     end
#     return zero(SMatrix{3,3,eltype(system),9})
# end

#--- setters ---#

"""
    body_to_multipole!(system::{UserDefinedSystem}, multipole_coefficients, buffer, expansion_center, bodies_index, harmonics, expansion_order)

Calculates the multipole coefficients due to the bodies contained in `buffer[:,bodies_index]` and accumulates them in `multipole_coefficients`. Should be overloaded for each user-defined system object (where `{UserDefinedSystem}` is replaced with the type of the user-defined system).

Typically, this is done using one of the convience functions contained within FastMultipole in one line, as

```julia
body_to_multipole!(system::MySystem, args...) = body_to_multipole!(Point{Vortex}, system, args...)
```

"""
function body_to_multipole!(system, multipole_coefficients, buffer, expansion_center, bodies_index, harmonics, expansion_order)
    if WARNING_FLAG_B2M[]
        @warn "body_to_multipole! not overloaded for type $(typeof(system)); multipole expansions from this system ignored"
        WARNING_FLAG_B2M[] = false
    end
    return nothing
end

"""
    direct!(target_system, target_index, derivatives_switch::DerivativesSwitch{PS,VS,GS}, ::{UserDefinedSystem}, source_buffer, source_index) where {PS,VS,GS}

Calculates direct (nearfield) interactions of `source_system` on `target_system`. Should be overloaded or each user-defined system object (where `{UserDefinedSystem}` is replaced with the type of the user-defined system), for all source bodies in `source_index`, at all target bodies in `target_index`, as follows:

```julia
# loop over source bodies
for i_source in source_index

    # extract source body information here...

    # loop over target bodies
    for i_target in target_index

        # get target position
        target_position = get_position(target_system, i_target)

        # evaluate influence here...

        # update appropriate quantities
        if PS
            set_scalar_potential!(target_system, i_target, scalar_potential)
        end
        if VS
            set_velocity!(target_system, i_target, velocity)
        end
        if GS
            set_velocity_gradient!(target_system, i_target, velocity_gradient)
        end

    end
end
```

Note that `::{UserDefinedSystem}` is used purely for overloading the method for the appropriate system, and should NOT be accessed in this function, since it will NOT be indexed according to `source_index`. Rather, `source_buffer`, which is updated using `source_system_to_buffer!`, should be accessed.

"""
function direct!(target_system, target_index, derivatives_switch, source_system, source_buffer, source_index)
    if WARNING_FLAG_DIRECT[]
        @warn "direct! not overloaded for type $(typeof(source_system)); interaction ignored"
        WARNING_FLAG_DIRECT[] = false
    end
    return nothing
end

"""
    buffer_to_target_system!(target_system::{UserDefinedSystem}, i_target, target_buffer, i_buffer)

Compatibility function used to update target systems. It should be overloaded for each system (where `{UserDefinedSystem}` is replaced with the type of the user-defined system) to be a target and should behave as follows. For the `i_body`th body contained inside of `target_system`,

* `target_buffer[4, i_buffer]` contains the scalar potential influence to be added to the `i_target` body of `target_system`
* `target_buffer[5:7, i_buffer]` contains the velocity influence to be added to the `i_target` body of `target_system`
* `target_buffer[8:16, i_buffer]` contains the velocity gradient to be added to the `i_target` body of `target_system`

Note that any system acting only as a source (and not as a target) need not overload `buffer_to_target_system!`.

"""
function buffer_to_target_system!(target_system, i_target, derivatives_switch, target_buffer, i_buffer)

    throw("buffer_to_target_system! not overloaded for type $(typeof(target_system))")
end

#------- internal functions -------#

#--- source_buffer getters ---#

# function get_position(source_buffer, source_system, i_body::Int)
#     return SVector{3}(view(source_buffer, 1:3, i_body))
# end

function get_radius(source_buffer::Matrix, i_body)
    return source_buffer[4,i_body]
end

function get_strength(source_buffer, source_system, i_body::Int)
    strength = SVector{strength_dims(source_system)}(view(source_buffer, 5:4+strength_dims(source_system), i_body))
    return strength
end

function get_vertex(source_buffer, source_system, i_body::Int, i_vertex::Int)
    i_offset = 3 * (i_vertex - 1)
    vertex = SVector{3}(view(source_buffer, 5+strength_dims(source_system)+i_offset:7+strength_dims(source_system)+i_offset, i_body))
    return vertex
end

#--- system/buffer setters ---#

function buffer_to_target!(target_systems::Tuple, target_tree::Tree, derivatives_switches=DerivativesSwitch(true, true, true, target_systems))
    buffer_to_target!(target_systems, target_tree.buffers, derivatives_switches, target_tree.sort_index_list)
end

function buffer_to_target!(target_systems::Tuple, target_buffers, derivatives_switches, sort_index_list=Tuple(1:get_n_bodies(system) for system in target_systems))
    for (target_system, target_buffer, derivatives_switch, sort_index) in zip(target_systems, target_buffers, derivatives_switches, sort_index_list)
        buffer_to_target!(target_system, target_buffer, derivatives_switch, sort_index)
    end
end

function buffer_to_target!(target_system, target_buffer, derivatives_switch, sort_index=1:get_n_bodies(target_system))
    for i_body in 1:get_n_bodies(target_system)
        buffer_to_target_system!(target_system, sort_index[i_body], derivatives_switch, target_buffer, i_body) # TODO: check this
    end
end

function target_to_buffer!(buffers, systems::Tuple, sort_index_list=SVector{length(systems)}([1:get_n_bodies(system) for system in systems]))
    for (buffer,system,sort_index) in zip(buffers, systems, sort_index_list)
        target_to_buffer!(buffer, system, sort_index)
    end
end

function target_to_buffer!(buffer::Matrix, system, sort_index=1:get_n_bodies(system))
    for i_body in 1:get_n_bodies(system)
        buffer[1:3, i_body] .= get_position(system, sort_index[i_body])
    end
end

function target_to_buffer(systems::Tuple, sort_index_list=SVector{length(systems)}([1:get_n_bodies(system) for system in systems]))
    buffers = allocate_buffers(systems, true)
    target_to_buffer!(buffers, systems, sort_index_list)
    return buffers
end

function target_to_buffer(system, sort_index=1:get_n_bodies(system))
    buffer = allocate_target_buffer(eltype(system), system)
    target_to_buffer!(buffer, system, sort_index)
    return buffer
end

function system_to_buffer!(buffers, systems::Tuple, sort_index_list=SVector{length(systems)}([1:get_n_bodies(system) for system in systems]))
    for (buffer, system, sort_index) in zip(buffers, systems, sort_index_list)
        system_to_buffer!(buffer, system, sort_index)
    end
end

function system_to_buffer!(buffer::Matrix, system, sort_index=1:get_n_bodies(system))
    for i_body in 1:get_n_bodies(system)
        source_system_to_buffer!(buffer, i_body, system, sort_index[i_body]) # TODO: check this
    end
end

function system_to_buffer(systems::Tuple, sort_index_list=SVector{length(systems)}([1:get_n_bodies(system) for system in systems]))
    buffers = allocate_buffers(systems, false)
    system_to_buffer!(buffers, systems, sort_index_list)
    return buffers
end

function system_to_buffer(system, sort_index=1:get_n_bodies(system))
    buffer = allocate_source_buffer(eltype(system), system)
    system_to_buffer!(buffer, system, sort_index)
    return buffer
end

#--- auxilliary functions ---#

@inline function get_n_bodies(systems::Tuple)
    n_bodies = 0
    for system in systems
        n_bodies += get_n_bodies(system)
    end
    return n_bodies
end

#------- access functions for use with a matrix of targets used as input to direct! -------#

#--- getters ---#

# function get_position(source_buffer, source_system, i_body::Int)
#     return SVector{3}(view(source_buffer, 1:3, i_body))
# end

get_position(system::Matrix{TF}, i) where TF = SVector{3,TF}(system[1, i], system[2, i], system[3, i])

get_scalar_potential(system::Matrix, i) = system[4, i]

get_velocity(system::Matrix{TF}, i) where TF = SVector{3,TF}(system[5,i], system[6,i], system[7,i])

get_velocity_gradient(system::Matrix{TF}, i) where TF =
    SMatrix{3,3,TF,9}(system[8, i], system[9, i], system[10, i],
    system[11, i], system[12, i], system[13, i],
    system[14, i], system[15, i], system[16, i])

get_n_bodies(sys::Matrix) = size(sys, 2)

#--- setters ---#

"""
    set_scalar_potential!(target_buffer, i_body, scalar_potential)

Accumulates `scalar_potential` to `target_buffer`.

"""
function set_scalar_potential!(system::Matrix, i, scalar_potential)
    system[4, i] += scalar_potential
end

"""
    set_velocity!(target_buffer, i_body, velocity)

Accumulates `velocity` to `target_buffer`.

"""
function set_velocity!(system::Matrix, i, velocity)
    system[5:7, i] .+= velocity
end

"""
    set_velocity_gradient!(target_buffer, i_body, velocity_gradient)

Accumulates `velocity_gradient` to `target_buffer`.

"""
function set_velocity_gradient!(system::Matrix, i, velocity_gradient)
    system[8, i] += velocity_gradient[1]
    system[9, i] += velocity_gradient[2]
    system[10, i] += velocity_gradient[3]
    system[11, i] += velocity_gradient[4]
    system[12, i] += velocity_gradient[5]
    system[13, i] += velocity_gradient[6]
    system[14, i] += velocity_gradient[7]
    system[15, i] += velocity_gradient[8]
    system[16, i] += velocity_gradient[9]
end

#--- auxilliary functions ---#

function reset!(systems::Tuple)
    for system in systems
        reset!(system)
    end
end

function reset!(system::Matrix)
    system[4:16, :] .= zero(eltype(system))
end
