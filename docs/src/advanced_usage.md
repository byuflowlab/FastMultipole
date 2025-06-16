# Advanced Usage

The most advanced features of `FastMultipole` include:

- the ability to solve the $n$-body problem for multiple systems simultaneously in a single FMM call
- satisfy an error tolerance by dynamically adjusting the expansion order of multipole-to-local transformations
- predict the full array of FMM tuning parameters such that an error tolerance is met
In this section, we will first describe how to run the FMM for multi-system problems. Then, we will describe how to impose an error tolerance. Finally, we will show how to automatically tune the FMM parameters.

## Simultaneous Computation of Multiple and Distinct Source and Target Systems

`FastMutipole` allows FMM to be performed efficiently on distinct source and target systems. This is done by creating two octrees: one for sources, and one for targets [yokota2013fmm2](@cite). For example, say we desire to know the influence of one system of point masses on another, but not vice versa. This is done by passing both systems to the `fmm!` function with the target system first, as:

```@example guidedex
using FastMultipole
using Random # needed for `gravitational.jl`

gravitational_path = normpath(joinpath(splitdir(pathof(FastMultipole))[1], "..", "test", "gravitational.jl"))
include(gravitational_path)

target_system = generate_gravitational(123, 1000)
source_system = generate_gravitational(321, 1000)

fmm!(target_system, source_system)
```

In practice, the source system might be a collection of systems, composed of a variety of datastructures. So long as the interface functions are defined for each system, we can pass a tuple of any number of systems as the source and or target. For example, we could evaluate the influence of a system of point masses, point vortices, and vortex filaments on two different target systems:

```@example guidedex
using LinearAlgebra

# include vortex filament and particle models and interface functions
vortex_path = normpath(joinpath(splitdir(pathof(FastMultipole))[1], "..", "test", "vortex.jl"))
include(vortex_path)
filament_path = normpath(joinpath(splitdir(pathof(FastMultipole))[1], "..", "test", "vortex_filament.jl"))
include(filament_path)

# generate systems
n_bodies = 2000
target_one = generate_gravitational(123, n_bodies)
target_two = generate_vortex(124, n_bodies)
source_one = generate_gravitational(125, n_bodies)
source_two = generate_vortex(126, n_bodies)
source_three = generate_filament_field(n_bodies, n_bodies^0.333, 127; strength_scale=1/n_bodies)

# run FMM
fmm!((target_one, target_two), (source_one, source_two, source_three); lamb_helmholtz=true,
    scalar_potential=false, gradient=true, hessian=(false, true))

v1 = target_one.potential[5:7,:]
v2 = target_two.gradient_stretching[1:3,:]

# run direct
target_one.potential .= 0.0
target_two.potential .= 0.0
target_two.gradient_stretching .= 0.0
direct!((target_one, target_two), (source_one, source_two, source_three); scalar_potential=false, gradient=true, hessian=(false, true))

# test accuracy
println("Max error in target one: ", maximum(abs.(target_one.potential[5:7,:] .- v1)))
println("Max error in target two: ", maximum(abs.(target_two.gradient_stretching[1:3,:] .- v2)))
```
Note that `scalar_potential`, `gradient`, and `hessian` can be passed as a single boolean or as a tuple of booleans, one for each target system. This allows the user to specify which values are desired for each target system, and avoids unnecessary calculations for values that are not needed. In this case, we have set `scalar_potential=false` and `gradient=true` to indicate all target systems, but a tuple `hessian=(false,true)` to indicate different settings for each. It is worth remembering that these switches must be implemented by the user when overloading the `direct!` function for each system to act as a source.

## Satisfying an Error Tolerance

`FastMultipole` can be configured to satisfy an error tolerance by dynamically adjusting the expansion order of each multipole-to-local transformation to the smallest integer that satisfies the tolerance, according to an error prediction. Users can indicate their desired error tolerance via the keyword argument `ε_tol` in the `fmm!` function. A value of `nothing` indicates no error tolerance, in which case the expansion order is fixed at whatever is passed as the `expansion_order` keyword. Otherwise, `ε_tol` should inherit from the `ErrorMethod` abstract type. The chosen type will determine how the error is predicted, and hence how the expansion order is chosen. Choices include:

- `PowerAbsolutePotential{tolerance}`: Constrains the magnitude of the potential error to be less than `tolerance` using a radially invariant upper bound.
- `PowerAbsoluteGradient{tolerance}`: Constrains the magnitude of the vector field error to be less than `tolerance` using a radially invariant upper bound.
- `PowerRelativePotential{tolerance}`: Constrains the relative error of the potential to be less than `tolerance` using a radially invariant upper bound.
- `PowerRelativeGradient{tolerance}`: Constrains the relative error of the vector field to be less than `tolerance` using a radially invariant upper bound.

Say I wanted to compute the gravitational potential to a tolerance of `1e-4` using the absolute potential error method. I would call the FMM as follows:

```@example guidedex
# create system
n_bodies, rand_seed = 50_000, 123
system = generate_gravitational(rand_seed, n_bodies)

# run FMM with error tolerance
fmm!(system; lamb_helmholtz=false, scalar_potential=true, gradient=false, hessian=false, ε_tol=PowerAbsolutePotential(1e-4))
```
This keyword argument `ε_tol=PowerAbsolutePotential(1e-4)` requests that the expansion order be dynamically adjusted to ensure that the potential error of each interaction is less than `1e-4` at each body. Note that this does ensure a perfectly conservative error bound, since many multipole-to-local transformations will be performed to any given target cell. However, the maximum error should be within an order of magnitude of the specified tolerance, and the average error will be much lower.

```@example guidedex
# verify error
phi_fmm = system.potential[1,:]
system.potential .= 0.0
direct!(system; scalar_potential=true, gradient=false, hessian=false)
phi_direct = system.potential[1,:]
println("Max error in potential: ", maximum(abs.(phi_direct .- phi_fmm)))
```
It is also worth noting that the error methods used here sometimes fail for coarse tolerances. Let's see what happens when we try to use a tolerance of `1e-2`:

```@example guidedex
# run FMM with error tolerance
system.potential .= 0.0
fmm!(system; lamb_helmholtz=false, scalar_potential=true, gradient=false, hessian=false, ε_tol=PowerAbsolutePotential(1e-2))

# check error
phi_fmm = system.potential[1,:]
println("Max error in potential: ", maximum(abs.(phi_direct .- phi_fmm)))
```
As a rule of thumb, care should be taken if 2 digits of accuracy or less is requested. The methods employed are fairly robust for 3 or more digits of accuracy. `PowerAbsolutePotential` and `PowerAbsoluteGradient` are best when the maximum error is to be constrained, while `RotatedCoefficientsAbsoluteGradient` is better if the mean error is to be constrained.


## Fully Automatic Tuning of FMM Parameters

`FastMultipole` can automatically tune the FMM parameters to satisfy an error tolerance. This is done by running the [`tune_fmm`](@ref) function, which returns a named tuple of the optimal parameters and a preallocated buffer to reduce memory allocations.

Say I wanted the optimal tuning parameters for the gravitational potential of a system of point masses, with a maximum error tolerance of `1e-4`. I would call the function as follows:

```@example guidedex
opt_params, cache = tune_fmm(system; ε_tol=PowerAbsolutePotential(1e-4), lamb_helmholtz=false, scalar_potential=true, gradient=false, hessian=false)
println("Optimal parameters: ", opt_params)
```
This will return a named tuple of the optimal parameters, which can then be passed to the `fmm!` function. The `cache` is a preallocated buffer that can be used to reduce memory allocations during the FMM call.

```@example guidedex
# run FMM without default parameters
system.potential .= 0.0
@time fmm!(system; lamb_helmholtz=false, scalar_potential=true, gradient=false, hessian=false, ε_tol=PowerAbsolutePotential(1e-4))

# verify error
phi_fmm = system.potential[1,:]
println("Max error in potential without tuning: ", maximum(abs.(phi_direct .- phi_fmm)))

# run FMM with optimal parameters
@time fmm!(system; lamb_helmholtz=false, scalar_potential=true, gradient=false, hessian=false, ε_tol=PowerAbsolutePotential(1e-4), cache..., opt_params...)

# verify error
println("Max error in potential: ", maximum(abs.(phi_direct .- phi_fmm)))
```
Note that the `tune_fmm` function iterates over each requested `multipole_acceptance` criterion, so it is not optimal to call it before every FMM call. Instead, `tune_fmm` can be called once for a representative system primarily to choose `multipole_acceptance`. Since the optimal parameters depend on the size, distribution, and strengths of the sources and targets, it is recommended to set the keyword argument `tune=true` to iteratively update the `expansion_order` and `leaf_size_source` parameters each time `fmm!` is called. This will allow the FMM to adapt to the system as it evolves, and will ensure that the optimal parameters are used for each call.
