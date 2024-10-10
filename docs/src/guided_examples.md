# Guided Examples

`FastMultipole` is designed to incorporate easily into your existing Julia code with minimal effort. We'll walk through the process of adding `FastMultipole` to an existing code in the following guided examples.

## Gravitational Example

In order to use `FastMultipole`, certain interface functions must be defined. In this example, we will walk through the process of incorporating FMM into the gravitational point mass model used in [Quick Start](@ref). This code can also be found under `test/gravitational.jl`. Note that an additional (minimalistic) vortex particle example is included in `test/vortex.jl`.

To better understand how the `FastMultipole` interface functions, let's take a look at the data structures we'll use to define our point masses:

```@example guidedex
using FastMultipole
using FastMultipole.StaticArrays

const i_POTENTIAL = 1:4             # index of the gravitational potential
const i_VELOCITY = 5:7              # index of the velocity
const i_VELOCITY_GRADIENT = 8:16    # index of the velocity gradient

# a single point mass
struct Body{TF}
    position::SVector{3,TF}
    radius::TF
    strength::SVector{4,TF}
end

# container for a system of `Body`'s
struct Gravitational{TF}
    bodies::Vector{Body{TF}}
    potential::Matrix{TF}
end

# constructor
function Gravitational(bodies::Matrix)
    nbodies = size(bodies)[2]
    bodies2 = [Body(SVector{3}(bodies[1:3,i]),bodies[4,i],SVector{4}(bodies[5:8,i])) for i in 1:nbodies]
    potential = zeros(52,nbodies)
    return Gravitational(bodies2,potential)
end
```

Note that while the particular choice of data structures is completely arbitrary, `FastMultipole` does not create any new storage containers for values of interest such as potential, velocity, velocity gradient. If these are desired, they must exist in the user-defined system, and will be updated in-place by `FastMultipole`.

### Overloading the `body_to_multipole!` Function

The function `body_to_multipole!` is used to generate multipole expansions for our particular system. This is done be overloading `body_to_multipole!` for the data structure representing our model. Convenience functions exist within `FastMultipole` to make this complicated function into a one-liner:

```@example guidedex
FastMultipole.body_to_multipole!(system::Gravitational, args...) =
    FastMultipole.body_to_multipole!(Point{Source}, system, args...)
```

Note that we are overloading `body_to_multipole!` to operate on our `::Gravitational` system, which consists of point sources. Other convenience functions exist in `FastMultipole` for any combination of `Point`, `Filament`, or `Panel` geometries using `Source`, `Dipole`, or `Vortex` kernels. These are specified when overloading `body_to_multipole!` as `<geometry>{<kernel>}`. For example, if my model used vortex filaments, I would replace `Point{Source}` in the example above to `Filament{Vortex}`.

We use the fast recursive method of generating exact coefficients for panels as developed by [GUMEROV2023112118](@cite).

### Overloading Getters

The `Gravitational` struct needs the following getters to be overloaded to support the indexing format used by `FastMultipole`. This approach allows for great flexibility without harming performance.

```@example guidedex
import Base: getindex

Base.getindex(g::Gravitational, i, ::FastMultipole.Position) = g.bodies[i].position
Base.getindex(g::Gravitational, i, ::FastMultipole.Radius) = g.bodies[i].radius
Base.getindex(g::Gravitational, i, ::FastMultipole.ScalarPotential) = g.potential[1,i]
Base.getindex(g::Gravitational, i, ::FastMultipole.Velocity) = view(g.potential,i_VELOCITY,i)
Base.getindex(g::Gravitational, i, ::FastMultipole.VelocityGradient) = reshape(view(g.potential,i_VELOCITY_GRADIENT,i),3,3)
Base.getindex(g::Gravitational, i, ::FastMultipole.Strength) = g.bodies[i].strength[1]
Base.getindex(g::Gravitational, i, ::FastMultipole.Body) = g.bodies[i], view(g.potential,:,i)
```
It is worth noting that there are ways of defining these functions that would harm performance, e.g. by allocating an array each time the velocity is requested. It is up to the user to define these functions with the efficiency they desire.

### Overloading Setters

`Gravitational` also needs the following setters to be overloaded as well. These are used to update the potential and other quantities in-place. The same performance considerations that applied to getters apply here.

```@example guidedex
import Base: setindex!

function Base.setindex!(g::Gravitational, val, i, ::FastMultipole.Body)
    body, potential = val
    g.bodies[i] = body
    g.potential[:,i] .= potential
    return nothing
end
function Base.setindex!(g::Gravitational, val, i, ::FastMultipole.ScalarPotential)
    g.potential[i_POTENTIAL[1],i] = val
end
function Base.setindex!(g::Gravitational, val, i, ::FastMultipole.Velocity)
    g.potential[i_VELOCITY,i] .= val
end
function Base.setindex!(g::Gravitational, val, i, ::FastMultipole.VelocityGradient)
    reshape(g.potential[i_VELOCITY_GRADIENT,i],3,3) .= val
end
```

### Additional Requirements

In addition to the getters and setters listed above, each system struct must be overloaded with three additional methods. In `gravitational.jl`, these are they are overloaded as follows.

```@example guidedex
# determine the number of bodies contained by the system
FastMultipole.get_n_bodies(g::Gravitational) = length(g.bodies)

# determine the float-type used by the system
Base.eltype(::Gravitational{TF}) where TF = TF

# return all data structures required to define a single body
FastMultipole.buffer_element(g::Gravitational) = (deepcopy(g.bodies[1]),zeros(eltype(g),52))
```

Note that the last of these functions, `buffer_element`, can return a tuple of values in case multiple data structures are required to define a single body.

### Vector Potential

Ordinarily, a 3-dimensional vector potential, such as the stream function in 3-dimensional fluid dynamics, requires 3 separate expansions (1 for each dimension). In `FastMultipole`, however, we employ the Lamb-Helmholtz decomposition demonstrated by [gumerov2013efficient](@cite) to reduce the number of expansions required from 3 to 2, thus reducing computational cost.

### Non-required Functionality

Though not required to run `Fast Multipole`, the `direct!` and `save_vtk` functions are useful for debugging and visualization. Here are some examples of how this could be implemented.

```@example guidedex
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

```


## FMM Tuning Parameters

We can improve the accuracy of the FMM by altering the follwing input parameters: `multipole_threshold`, `leaf_size`, and `expansion_order`. Though not linear, `expansion_order` is positively correlated with accuracy while `multipole_threshold` and `leaf_size` are negatively correlated.

- **`multipole_threshold` (between 0 and 1) the radius of the source over the minimum distance at which multipole expansion are used**
- **`leaf_size` (greater than 1) the maximum number of bodies included in each leaf level branch**
- **`expansion_order` the number of terms included in each multipole expansion minus one**

```@example guidedex
import FastMultipole as fmm
using Random

function generate_gravitational(seed, n_bodies; radius_factor=0.1, strength_factor=1.0)
    Random.seed!(seed)
    bodies = rand(8,n_bodies)
    bodies[1:3,:] = rand(3,n_bodies) # body positions

    bodies[4,:] ./= (n_bodies^(1/3)*2) # body radii
    bodies[4,:] .*= radius_factor

    bodies[5,:] .*= strength_factor # body strengths

    system = Gravitational(bodies)
    return system
end

function measure_error(;multipole_threshold=0.4, leaf_size=50, expansion_order=5)
    fmm_system = generate_gravitational(123,1000)
    direct_system = deepcopy(fmm_system)

    fmm.fmm!(fmm_system; multipole_threshold, leaf_size, expansion_order)
    fmm.direct!(direct_system)

    percent_error = abs.((fmm_system.potential[1,:] .- direct_system.potential[1,:]) ./ direct_system.potential[1,:])
    return maximum(percent_error)
end

# default parameters
println(measure_error(multipole_threshold=0.4, leaf_size=50, expansion_order=5))

# expansion_order increased to 10
println(measure_error(multipole_threshold=0.4, leaf_size=50, expansion_order=10))
```

## Evaluate Source on Target

`FastMutipole` allows for the evaluation of source systems on target systems while leaving the source systems unaltered.

```@example guidedex
target_system = generate_gravitational(123, 100)
source_system = generate_gravitational(321, 100)

fmm.fmm!(target_system, source_system)
```

## Evaluate Multiple Sources on Multiple Targets

The FMM also supports the evaluation of multiple source systems on multiple target systems. The user is also able to evaluate a single source on multiple targets or multiple sources on a single target with any combination of supported system types.

```@example guidedex
vortex_path = normpath(joinpath(splitdir(pathof(fmm))[1], "..", "test", "vortex.jl"))
include(vortex_path)

function generate_vortex(seed, n_bodies; radius_factor=0.1, strength_factor=1.0)
    Random.seed!(seed)
    bodies = rand(8,n_bodies)
    bodies[1:3,:] = rand(3,n_bodies) # body positions

    bodies[4,:] ./= (n_bodies^(1/3)*2) # body radii
    bodies[4,:] .*= radius_factor

    bodies[5,:] .*= strength_factor # body strengths

    system = VortexParticles(bodies)
    return system
end

target_one = generate_gravitational(123, 100)
target_two = generate_vortex(124, 100)

source_one = generate_gravitational(125, 100)
source_two = generate_vortex(126, 100)

fmm.fmm!((target_one, target_two), (source_one, source_two))
```

## Non-Potential Flow Applications

As a default, target systems will be evaluated and returned with `scalar_potential`, `velocity`, and `velocity_gradient` fields populated. (Note that the vector potential is not explicitly tracked, though its induced velocity and velocity gradient are.) In some situations, only some of these values may be required. By inputting a boolean vector of the same length as target systems, the user is able to speed up the calculation by not storing unecessary values.

$\overline{V} = -\nabla \phi + \nabla \times \overline{\psi}$

|  | $\phi$ | $\overline{\psi}$ | $\overline{V}$| $\nabla \overline{V}$ |
| ----------- | ----------- | ----------- | ----------- | ----------- |
| `fmm!` Keyword Arguments | `scalar_potential` | `vector_potential` | `velocity` | `velocity_gradient` |
| Fluid Dynamics | Scalar Potential | Stream Function | Fluid Velocity | Velocity Gradient |
| Electrostatics | Electric Potential | - | Electric Field | Field Gradient Tensor |
| Magnetostatics | - | Magnetic Vector Potential | Magnetic Field | Field Gradient Tensor |
| Gravity | Gravitational Potential | - | Gravitational Acceleration | Acceleration Gradient Tensor |


```@example guidedex
target_one = generate_gravitational(123, 100)
target_two = generate_vortex(124, 100)

source_one = generate_gravitational(125, 100)

fmm.fmm!((target_one, target_two), source_one, scalar_potential=[true, false],
    velocity=[true, true], velocity_gradient=[false, false])
```

## Saving Generated Trees

A given `fmm!` call will typically return a single tree (if performed on the entire system) or two seperate source/target trees (if called on a source and a target). The function call can also be modified to save these trees.

```@example guidedex

target_filepath = "target_tree"
source_filepath = "source_tree"

target_tree, source_tree = fmm.fmm!((target_one, target_two), source_one;
    save_tree_target=true, save_name_target=target_filepath,
    save_tree_source=true, save_name_source=source_filepath)
```
