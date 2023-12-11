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
Base.getindex(sys, i, ::ScalarStrength) = 0.0
Base.getindex(sys, i, ::VectorStrength) = SVector{3}(0.0,0.0,0.0)
function Base.setindex!(sys, val, i)
    @warn "setindex! not overloaded for type $(typeof(sys)); octree will have errors"
    return nothing
end
function Base.setindex!(sys, val, i, ::ScalarPotential)
    # @warn "setindex! (ScalarPotential) not overloaded for type $(typeof(sys))"
    return nothing
end
function Base.setindex!(sys, val, i, ::VectorPotential)
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
function Base.length(sys)
    @error "Base.length() not overloaded for type $(typeof(sys))"
end
# Base.eltype(::SystemType) = 

function buffer_element(system)
    @error "`buffer_element` not overloaded for system type $(typeof(system))"
end

# function buffer_element(system)
#     this_buffer = deepcopy(system[1])
#     view_free = true
#     view_free = buffer_element(this_buffer, view_free)
#     !view_free && (error("`buffer_element` overloaded to return a `view` for system type $(typeof(system)); errors will occur in the `unsort!` function\n\tTo fix this error, overload `buffer_element(system::$(typeof(system)))` to contain arrays rather than views"))
#     return this_buffer
# end

# function buffer_element(member, view_free)
#     try
#         for submember in member
#             view_free *= buffer_element(submember, view_free)
#         end
#     finally
#         return view_free
#     end
# end

# function buffer_element(member::SubArray, view_free)
#     return false
# end

function B2M!(system, branch, bodies_index, harmonics, expansion_order)
    @warn "B2M! function not overloaded for type $(typeof(system)); interaction ignored"
    return nothing
end

function direct!(target_system, target_index, source_system, source_index)
    @warn "direct! function not overloaded for type $(typeof(source_system)); interaction ignored"
    return nothing
end

#####
##### multipole generation convenience functions
#####
@inline function B2M!_vortexpoint(system, branch, bodies_index, harmonics, expansion_order)
    c_x, c_y, c_z = branch.center
    for i_body in bodies_index
        b_x, b_y, b_z = system[i_body,POSITION]
        x = b_x - c_x
        z = b_z - c_z
        y = b_y - c_y
        qx, qy, qz = system[i_body,VECTOR_STRENGTH]
        r, theta, phi = cartesian_2_spherical(x,y,z)
        regular_harmonic!(harmonics, r ,theta, -phi, expansion_order) # Ylm^* -> -dx[3]
        # update values
        for l in 0:expansion_order
            for m in 0:l
                i_solid_harmonic = l*l + l + m + 1
                i_compressed = 1 + (l * (l + 1)) >> 1 + m # only save half as Yl{-m} = conj(Ylm)
                #branch.multipole_expansion[2,i_compressed] += harmonics[i_solid_harmonic] * qx * ONE_OVER_4PI
                #branch.multipole_expansion[3,i_compressed] += harmonics[i_solid_harmonic] * qy * ONE_OVER_4PI
                #branch.multipole_expansion[4,i_compressed] += harmonics[i_solid_harmonic] * qz * ONE_OVER_4PI
                branch.multipole_expansion[1,2,i_compressed] += harmonics[1,i_solid_harmonic] * qx * ONE_OVER_4PI
                branch.multipole_expansion[2,2,i_compressed] += harmonics[2,i_solid_harmonic] * qx * ONE_OVER_4PI
                branch.multipole_expansion[1,3,i_compressed] += harmonics[1,i_solid_harmonic] * qy * ONE_OVER_4PI
                branch.multipole_expansion[2,3,i_compressed] += harmonics[2,i_solid_harmonic] * qy * ONE_OVER_4PI
                branch.multipole_expansion[1,4,i_compressed] += harmonics[1,i_solid_harmonic] * qz * ONE_OVER_4PI
                branch.multipole_expansion[2,4,i_compressed] += harmonics[2,i_solid_harmonic] * qz * ONE_OVER_4PI
            end
        end
    end
end

@inline function B2M!_sourcepoint(system, branch, bodies_index, harmonics, expansion_order)
    c_x, c_y, c_z = branch.center
    for i_body in bodies_index
        b_x, b_y, b_z = system[i_body,POSITION]
        x = b_x - c_x
        z = b_z - c_z
        y = b_y - c_y
        q = system[i_body,SCALAR_STRENGTH]
        r, theta, phi = cartesian_2_spherical(x,y,z)
        #l = length(ReverseDiff.tape(r))
        harmonics .= regular_harmonic!(harmonics, r, theta, -phi, expansion_order) # Ylm^* -> -dx[3]
        # update values
        for l in 0:expansion_order
            for m in 0:l
                i_solid_harmonic = l*l + l + m + 1
                i_compressed = 1 + (l * (l + 1)) >> 1 + m # only save half as Yl{-m} = conj(Ylm)
                #branch.multipole_expansion[1,i_compressed] += harmonics[i_solid_harmonic] * q
                branch.multipole_expansion[1,1,i_compressed] += harmonics[1,i_solid_harmonic] * q
                branch.multipole_expansion[2,1,i_compressed] += harmonics[2,i_solid_harmonic] * q
            end
        end
        #@show length(ReverseDiff.tape(r)) - l # 451 allocations (per particle)
    end
end
