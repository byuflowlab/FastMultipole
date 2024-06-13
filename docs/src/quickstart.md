# Quick Start

The tutorial for FastMultipole is a work in progress. While we work on creating a tutorial, see the example file `scripts/simple_gravitational.jl`.



## Create a System
We can create a system of 10 randomly spaced bodies.
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

system = generate_gravitational(123, 10)
```

## Evaluate The Potential at Each Body
The `fmm!` call evaluates the potential of `system` in place.

```@example ex
    fmm.fmm!(system)
    @show system.potential[1,:]
```
This potential can then be accessed as a member variable of `system`.

## Accuracy of FMM Call
By using the `direct!` function, we can check the accuracy of a given `fmm!` call.

```@example ex
vortex_path = normpath(joinpath(splitdir(pathof(fmm))[1], "..", "test", "vortex.jl"))
include(vortex_path)

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

fmm_system = generate_gravitational(123,1000)
direct_system = deepcopy(fmm_system)

tree = fmm.fmm!(fmm_system)
fmm.direct!(direct_system)
percent_error = abs.((fmm_system.potential[1,:] .- direct_system.potential[1,:]) ./ direct_system.potential[1,:])


@show maximum(percent_error)
```