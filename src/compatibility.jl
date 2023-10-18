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
# Base.length(system) = 
# Base.eltype(::SystemType) = 

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
                branch.multipole_expansion[2,i_compressed] += harmonics[i_solid_harmonic] * qx * ONE_OVER_4PI
                branch.multipole_expansion[3,i_compressed] += harmonics[i_solid_harmonic] * qy * ONE_OVER_4PI
                branch.multipole_expansion[4,i_compressed] += harmonics[i_solid_harmonic] * qz * ONE_OVER_4PI
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
        regular_harmonic!(harmonics, r, theta, -phi, expansion_order) # Ylm^* -> -dx[3]
        # update values
        for l in 0:expansion_order
            for m in 0:l
                i_solid_harmonic = l*l + l + m + 1
                i_compressed = 1 + (l * (l + 1)) >> 1 + m # only save half as Yl{-m} = conj(Ylm)
                branch.multipole_expansion[1,i_compressed] += harmonics[i_solid_harmonic] * q
            end
        end
    end
end