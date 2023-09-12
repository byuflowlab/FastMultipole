import FLOWFMM as fmm
using WriteVTK
import Base: getindex, setindex!
using StaticArrays
const i_POSITION = 1:3
const i_RADIUS = 4
const i_STRENGTH = 5:8
const i_POTENTIAL = 1:4
const i_VELOCITY = 5:7
const i_VELOCITY_GRADIENT = 8:16

#####
##### gravitational kernel and mass elements
#####
struct Body{TF}
    position::SVector{3,TF}
    radius::TF
    strength::SVector{4,TF}
end

struct Gravitational{TF}
    bodies::Vector{Body{TF}}
    potential::Matrix{TF}
end

function Gravitational(bodies::Matrix)
    nbodies = size(bodies)[2]
    bodies2 = [Body(SVector{3}(bodies[1:3,i]),bodies[4,i],SVector{4}(bodies[5:8,i])) for i in 1:nbodies]
    potential = zeros(52,nbodies)
    return Gravitational(bodies2,potential)
end

Base.getindex(g::Gravitational, i, ::fmm.Position) = g.bodies[i].position
# Base.getindex(g::Gravitational, i, ::fmm.Velocity) = 
Base.getindex(g::Gravitational, i, ::fmm.Radius) = g.bodies[i].radius
# Base.getindex(g::Gravitational, i, ::fmm.Potential) = view(g.potential,1:4,i)
# Base.getindex(g::Gravitational, i, ::fmm.ScalarPotential) = view(g.potential,1,i)
# Base.getindex(g::Gravitational, i, ::fmm.Velocity) = view(g.potential,i_VELOCITY,i)
# Base.getindex(g::Gravitational, i, ::fmm.VelocityGradient) = reshape(view(g.potential,i_VELOCITY_GRADIENT,i),3,3)
Base.getindex(g::Gravitational, i) = g.bodies[i], view(g.potential,:,i)
# function Base.setindex!(g::Gravitational, val, i)
#     body, potential = val
#     g.bodies[i] = body
#     g.potential[:,i] .= potential
#     return nothing
# end
function Base.setindex!(g::Gravitational, val, i, ::fmm.ScalarPotential)
    g.potential[i_POTENTIAL[1],i] .= val
end
function Base.setindex!(g::Gravitational, val, i, ::fmm.Potential)
    g.potential[i_POTENTIAL,i] .= val
end
function Base.setindex!(g::Gravitational, val, i, ::fmm.Velocity)
    g.potential[i_VELOCITY,i] .= val
end
function Base.setindex!(g::Gravitational, val, i, ::fmm.VelocityGradient)
    reshape(g.potential[i_VELOCITY_GRADIENT,i],3,3) .= val
end
Base.length(g::Gravitational) = length(g.bodies)
Base.eltype(::Gravitational{TF}) where TF = TF

function fmm.B2M!(branch, system::Gravitational, bodies_index, harmonics, expansion_order)
    c_x, c_y, c_z = branch.center
    for i_body in bodies_index
        b_x, b_y, b_z = system[i_body,fmm.POSITION]
        x = b_x - c_x
        z = b_z - c_z
        y = b_y - c_y
        q = system.bodies[i_body].strength[1]
        r, theta, phi = fmm.cartesian_2_spherical(x,y,z)
        fmm.regular_harmonic!(harmonics, r, theta, -phi, expansion_order) # Ylm^* -> -dx[3]
        # update values
        for l in 0:expansion_order
            for m in 0:l
                i_solid_harmonic = l*l + l + m + 1
                i_compressed = 1 + (l * (l + 1)) >> 1 + m # only save half as Yl{-m} = conj(Ylm)
                branch.multipole_expansion[1][i_compressed] += harmonics[i_solid_harmonic] * q
            end
        end
    end
end

function fmm.direct!(target_system, target_index, source_system::Gravitational, source_index)
    for i_source in source_index
        source_x, source_y, source_z = source_system[i_source,fmm.POSITION]
        source_strength = source_system.bodies[i_source].strength[1]
        for j_target in target_index
            target_x, target_y, target_z = target_system[j_target,fmm.POSITION]
            dx = target_x - source_x
            dy = target_y - source_y
            dz = target_z - source_z
            r = sqrt(dx*dx + dy*dy + dz*dz)
            if r > 0
                dV = source_strength / r
                target_system[j_target,fmm.SCALAR_POTENTIAL] .+= dV
            end
        end
    end
end

function save_vtk(filename, element::Gravitational, nt=0; compress=false, extra_fields=nothing)
    _, n = size(element.bodies)
    WriteVTK.vtk_grid(filename*"_point_masses."*string(nt)*".vts", reshape(view(element.bodies,1:3,:),3,n,1,1); compress) do vtk
        vtk["strength"] = reshape(view(element.bodies,4,:), 1, n, 1, 1)
        vtk["velocity"] = reshape(element.velocity, 3, n, 1, 1)
        vtk["scalar potential"] = reshape(view(element.potential,1,:), n, 1, 1)
        vtk["vector potential"] = reshape(view(element.potential,2:4,:), 3, n, 1, 1)
        if !isnothing(extra_fields)
            for i in 1:length(extra_fields)
                vtk[extra_fields[i][1]] = extra_fields[i][2]
            end
        end
    end
end
