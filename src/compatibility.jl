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
<<<<<<< HEAD

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
=======
>>>>>>> main
