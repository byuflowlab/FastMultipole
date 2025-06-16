# Quick Start

The following tutorial shows how to use `FastMultipole` to compute the gravitational potential induced by a collection of point masses. It uses data structures located in `test/gravitational.jl`.

## Create a System

First, let's create a system of 1000 randomly spaced point masses:

```@example ex
using FastMultipole
using Random # needed for `gravitational.jl`

gravitational_path = normpath(joinpath(splitdir(pathof(FastMultipole))[1], "..", "test", "gravitational.jl"))
include(gravitational_path)

rand_seed = 123
n_bodies = 1000
system = generate_gravitational(rand_seed, n_bodies)
```

## Evaluate The Potential at Each Body

The `fmm!` function evaluates the gravitational potential induced by `system` in-place. We can control the tradeoff between performance and accuracy with a handful of tuning parameters, but we'll stick with the defaults for this example:

```@example ex; continued = true
fmm!(system)
```
The resulting potential can then be accessed in the user-defined `system` object.

```@example ex
@show system.potential[1,:]
```
Let me emphasize that `system` can be _any_ user-defined object, so long as a few interface functions are defined (we'll go over those later). This allows you to use `FastMultipole` in your existing code with almost no modifications.

## Accuracy of FMM Call

By using the `direct!` function, we can check the accuracy of the `fmm!` call by evaluating the ''N''-body problem naively, without fast multipole acceleration.

```@example ex
direct_system = deepcopy(system)
direct_system.potential .= 0.0

direct!(direct_system)

percent_error = abs.((system.potential[1,:] .- direct_system.potential[1,:]) ./ direct_system.potential[1,:])

@show maximum(percent_error)
```

## Scalar-plus-Vector Potential Applications

The use of a scalar-plus-vector potential is quite general, and is useful in a variety of physics and engineering contexts, including fluid dynamics, linear elasticity, electromagnetism, astrophysics, and others. We include a table of some common field quantities that derive from the scalar-plus-vector potential, along with their corresponding `fmm!` keyword arguments and physical interpretations.

$\vec{v} = \nabla \phi + \nabla \times \vec{\psi}$

|  | $\phi$ | $\vec{\psi}$ | $\vec{v}$| $\nabla \vec{v}$ |
| ----------- | ----------- | ----------- | ----------- | ----------- |
| `fmm!` Keyword Arguments | `scalar_potential` | `vector_potential` | `gradient` | `hessian` |
| Fluid Dynamics | Scalar Potential | Stream Function | Fluid Velocity | Velocity Gradient |
| Electrostatics | Electric Potential | - | Electric Field | Field Gradient Tensor |
| Magnetostatics | - | Magnetic Vector Potential | Magnetic Field | Field Gradient Tensor |
| Astrophysics | Gravitational Potential | - | Gravitational Acceleration | Acceleration Gradient Tensor |
