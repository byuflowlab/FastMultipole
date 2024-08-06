# Quick Start

The following tutorial shows how to use `FastMultipole` to compute the gravitational potential induced by a collection of point masses. It uses data structures located in `test/gravitational.jl`.

## Create a System

First, let's create a system of 1000 randomly spaced point masses:

```@example ex
using Random
import FastMultipole as fmm
gravitational_path = normpath(joinpath(splitdir(pathof(fmm))[1], "..", "test", "gravitational.jl"))
include(gravitational_path)

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

system = generate_gravitational(123, 1000)
```

## Evaluate The Potential at Each Body

The `fmm!` function evaluates the gravitational potential induced by `system` in-place. We can control the tradeoff between performance and accuracy by tuning a handful of parameters for our particular system, but we'll stick with the defaults for this example:

```@example ex; continued = true
fmm.fmm!(system)
```
The resulting potential can then be accessed as a field of `system`.

```@example ex
@show system.potential[1,:]
```

## Accuracy of FMM Call

By using the `direct!` function, we can check the accuracy of the `fmm!` call by evaluating the ''N''-body problem naively, without fast multipole acceleration.

```@example ex
direct_system = deepcopy(system)
direct_system.potential .= 0

fmm.direct!(direct_system)

percent_error = abs.((system.potential[1,:] .- direct_system.potential[1,:]) ./ direct_system.potential[1,:])

@show maximum(percent_error)
```
