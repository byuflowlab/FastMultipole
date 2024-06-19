import FastMultipole
using FastMultipole
using WriteVTK
import Base: getindex, setindex!
using FastMultipole.StaticArrays
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

Base.getindex(g::Gravitational, i, ::FastMultipole.Position) = g.bodies[i].position
Base.getindex(g::Gravitational, i, ::FastMultipole.Radius) = g.bodies[i].radius
Base.getindex(g::Gravitational, i, ::FastMultipole.VectorPotential) = view(g.potential,2:4,i)
Base.getindex(g::Gravitational, i, ::FastMultipole.ScalarPotential) = g.potential[1,i]
Base.getindex(g::Gravitational, i, ::FastMultipole.Velocity) = view(g.potential,i_VELOCITY,i)
Base.getindex(g::Gravitational, i, ::FastMultipole.VelocityGradient) = reshape(view(g.potential,i_VELOCITY_GRADIENT,i),3,3)
Base.getindex(g::Gravitational, i, ::FastMultipole.Strength) = g.bodies[i].strength[1]
Base.getindex(g::Gravitational, i, ::FastMultipole.Body) = g.bodies[i], view(g.potential,:,i)
function Base.setindex!(g::Gravitational, val, i, ::FastMultipole.Body)
    body, potential = val
    g.bodies[i] = body
    g.potential[:,i] .= potential
    return nothing
end
function Base.setindex!(g::Gravitational, val, i, ::FastMultipole.ScalarPotential)
    g.potential[i_POTENTIAL[1],i] = val
end
function Base.setindex!(g::Gravitational, val, i, ::FastMultipole.VectorPotential)
    g.potential[i_POTENTIAL[2:4],i] .= val
end
function Base.setindex!(g::Gravitational, val, i, ::FastMultipole.Velocity)
    g.potential[i_VELOCITY,i] .= val
end
function Base.setindex!(g::Gravitational, val, i, ::FastMultipole.VelocityGradient)
    reshape(g.potential[i_VELOCITY_GRADIENT,i],3,3) .= val
end
FastMultipole.get_n_bodies(g::Gravitational) = length(g.bodies)
Base.eltype(::Gravitational{TF}) where TF = TF

FastMultipole.buffer_element(g::Gravitational) = (deepcopy(g.bodies[1]),zeros(eltype(g),52))

FastMultipole.B2M!(system::Gravitational, args...) = FastMultipole.B2M!_sourcepoint(system, args...)

function FastMultipole.direct!(target_system, target_index, derivatives_switch, source_system::Gravitational, source_index)
    # nbad = 0
    for i_source in source_index
        source_x, source_y, source_z = source_system[i_source,FastMultipole.POSITION]
        source_strength = source_system.bodies[i_source].strength[1]
        for j_target in target_index
            target_x, target_y, target_z = target_system[j_target,FastMultipole.POSITION]
            dx = target_x - source_x
            dy = target_y - source_y
            dz = target_z - source_z
            r = sqrt(dx*dx + dy*dy + dz*dz)
            # te = @elapsed begin
            if r > 0
                dV = source_strength / r
                target_system[j_target,FastMultipole.SCALAR_POTENTIAL] += dV
            end
        # end
        # if te > 0.00001; nbad += 1; end
        end
    end
    # println("nbad = $nbad")
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
