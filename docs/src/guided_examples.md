# Guided Examples

`FastMultipole` is designed to incorporate easily into your existing Julia code with minimal effort. In particular, several interface functions must be defined. We'll walk through the process in the following guided examples.

## Gravitational Example

In this example, we review the interface functions used by the gravitational point mass model used in [Quick Start](@ref). This code can also be found under `FastMultipole/test/gravitational.jl`.

To better understand how the `FastMultipole` interface functions, let's take a look at the data structures we'll use to define our point masses:

```@example guidedex
using FastMultipole
using FastMultipole.StaticArrays

const i_POTENTIAL = 1:1   # index of the gravitational potential
const i_GRADIENT = 5:7    # index of the gravitational acceleration
const i_HESSIAN = 8:16    # index of the hessian matrix

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

In short, `Body` objects represent individual point masses, and the `Gravitational` object contains a group of them, along with a matrix of the potential, velocity, and velocity gradient at each body. The observant reader will notice that the `strength` field of `Body` is a 4-vector. The first index contains the mass of the body, used to induce a gravitational potential. The final three indices are meant to contain a vector strength, in case a vector potential is desired. For the sake of simplicy, we will assume `strength[2:4]` is zero for now. 

### Overloading `body_to_multipole!`

The function `body_to_multipole!` is used to generate multipole expansions for our particular system. This is done be overloading `FastMultipole.body_to_multipole!` for our data structure. Convenience functions exist within `FastMultipole` to simplify this complicated function into a one-liner:

```@example guidedex
FastMultipole.body_to_multipole!(system::Gravitational, args...) =
    FastMultipole.body_to_multipole!(Point{Source}, system, args...)
```

The `::Gravitational` type annotation is the key that allows `FastMultipole` to dispatch on `::Gravitational` systems. `Point{Source}` is used to indicate that our bodies are points and induce a source (i.e. ``1/r``) potential. Other convenience functions exist in `FastMultipole` for any combination of `Point`, `Filament`, or `Panel` geometries using `Source`, `Dipole`, or `Vortex` kernels, and are specified when overloading `body_to_multipole!` as `<geometry>{<kernel>}`. For example, if my model used constant vortex sheet panels, I would replace `Point{Source}` in the example above to `Panel{Vortex}`.

We use the fast recursive method of generating exact coefficients for panels as developed by [GUMEROV2023112118](@cite). We have derived our own formulae for vortex filaments and panels, which are in the process of publication.

### Buffer Interface Functions

`FastMultipole` allocates a buffer matrix in which to sort and store key system information. The buffer interface is created by overloading several functions for the user-defined system type. The first of these functions is [`FastMultipole.source_system_to_buffer!`](@ref), which populates a matrix with the important information about each body, column-wise. Overloading for `::Gravitational` systems looks like this:

```@example guidedex
function FastMultipole.source_system_to_buffer!(buffer, i_buffer, system::Gravitational, i_body)
    buffer[1:3, i_buffer] .= system.bodies[i_body].position
    buffer[4, i_buffer] = system.bodies[i_body].radius
    buffer[5, i_buffer] = system.bodies[i_body].strength[1]
end
```
Here, the `i_body`th body position, radius, and strength are copied to the correct positions in the `i_buffer`th column of the buffer matrix.

It is worth noting that there are ways of defining these functions that would harm performance, e.g. by allocating an array each time the velocity is requested. If you notice `FastMultipole`'s performance is poor, this and the other interface functions are good places to check.

Because every system can have a unique buffer structure, we need to overload a few additional functions to tell `FastMultipole` how to allocate and access it. The first function, `data_per_body`, returns the number of rows in the buffer matrix that are used to define a single body. The second function, `get_position`, returns the position of the `i`th body in the system. The third function, `strength_dims`, returns the number of rows in the buffer matrix that are used to define the strength of a single body. Finally, we overload `get_n_bodies` to return the number of bodies in our system. For our `Gravitational` system, we define:

```@example guidedex
Base.eltype(::Gravitational{TF}) where TF = TF

function FastMultipole.data_per_body(system::Gravitational)
    return 5
end

function FastMultipole.get_position(system::Gravitational, i)
    return system.bodies[i].position
end

function FastMultipole.strength_dims(system::Gravitational)
    return 1
end

FastMultipole.get_n_bodies(system::Gravitational) = length(system.bodies)
```
These functions provide enough information for `FastMultipole` to allocate the buffer matrix and access the data it needs.

Finally, we need to tell `FastMultipole` how to transfer the results of the FMM call back to the user-defined system. This is done by overloading the [`FastMultipole.buffer_to_target_system!`](@ref) function. Overloading for `::Gravitational` systems looks like this:

```@example guidedex
function FastMultipole.buffer_to_target_system!(target_system::Gravitational, i_target, ::FastMultipole.DerivativesSwitch{PS,VS,GS}, target_buffer, i_buffer) where {PS,VS,GS}
    # get values
    TF = eltype(target_buffer)
    scalar_potential = PS ? FastMultipole.get_scalar_potential(target_buffer, i_buffer) : zero(TF)
    velocity = VS ? FastMultipole.get_gradient(target_buffer, i_buffer) : zero(SVector{3,TF})
    hessian = GS ? FastMultipole.get_hessian(target_buffer, i_buffer) : zero(SMatrix{3,3,TF,9})

    # update system
    target_system.potential[i_POTENTIAL[1], i_target] = scalar_potential
    target_system.potential[i_GRADIENT, i_target] .= velocity
    for (jj,j) in enumerate(i_HESSIAN)
        target_system.potential[j, i_target] = hessian[jj]
    end
end
```
We note that the convenience functions `get_velocity` and `get_hessian` return `::SVector{3}` and `::SMatrix{3,3}`, respectively, to reduce allocations.

### Overloading `direct!`

The last required interface function is [`FastMultipole.direct!`](@ref), which evaluates the potential at a target system using a source system without multipole acceleration. This is required in an FMM call for those interactions that are too close to be approximated by expansions. It is also useful for debugging and testing the accuracy of the FMM call. Overloading for `::Gravitational` systems, we have:

```@example guidedex
function FastMultipole.direct!(target_system, target_index, ::DerivativesSwitch{PS,VS,GS}, source_system::Gravitational, source_buffer, source_index) where {PS,VS,GS}
    @inbounds for i_source in source_index
        source_x, source_y, source_z = FastMultipole.get_position(source_buffer, i_source)
        source_strength = FastMultipole.get_strength(source_buffer, source_system, i_source)[1]
        @inbounds for j_target in target_index
            target_x, target_y, target_z = FastMultipole.get_position(target_system, j_target)
            dx = target_x - source_x
            dy = target_y - source_y
            dz = target_z - source_z
            r2 = dx*dx + dy*dy + dz*dz
            if r2 > 0
                r = sqrt(r2)
                if PS
                    dϕ = source_strength / r * FastMultipole.ONE_OVER_4π
                    FastMultipole.set_scalar_potential!(target_system, j_target, dϕ)
                end
                if VS
                    dF = SVector{3}(dx,dy,dz) * source_strength / (r2 * r) * FastMultipole.ONE_OVER_4π
                    FastMultipole.set_gradient!(target_system, j_target, dF)
                end
            end
        end
    end
end
```
Note that the velocity gradient is not calculated in this function. If the velocity gradient is requested, the contribution due to `::Gravitational` systems will then be zero.

### Running the FMM

Now that we have defined the interface functions, we can run the FMM on our `Gravitational` system. The [`fmm!`](@ref) function is used to perform the FMM call. For example,

```@example guidedex
using Random

# create system
function generate_gravitational(seed, n_bodies; radius_factor=0.1, strength_scale=1/n_bodies, bodies_fun=(x)->x)
    Random.seed!(seed)
    bodies = rand(8,n_bodies)
    # bodies = rand(distribution,8,n_bodies)
    bodies[4,:] ./= (n_bodies^(1/3)*2)
    bodies[4,:] .*= radius_factor
    bodies[5,:] .*= strength_scale

    bodies_fun(bodies)

    system = Gravitational(bodies)
end

n_bodies, rand_seed = 1000, 123
system = generate_gravitational(rand_seed, n_bodies)

# run FMM
fmm!(system; lamb_helmholtz=false, scalar_potential=false, gradient=true, hessian=true)
```
We have supplied `lamb_helmholtz=false` to the FMM call because the gravitational potential is a scalar potential and does not require the Lamb-Helmholtz decomposition of a vector field. This is not strictly necessary, since the predicted potential would be identical regardless. However, it does provide a slight performance benefit by omitting the Lamb-Helmholtz operators. The `scalar_potential` and `gradient` arguments are set to `false` and `true`, respectively, to indicate that we want to compute the vector field (i.e. the gravitational field lines that translate to force and their gradient) but not the scalar potential.

### Preallocating the Buffers

The `fmm!` function incurs a small overhead for allocating the buffer matrices. If our system is large and the FMM will be called more than once, we can preallocate the buffers to avoid this overhead. This is done by returning the buffer cache the first time `fmm!` is called, and then splatting it as an optional argument to `fmm!` in subsequent calls as follows:

```@example guidedex
# allocate buffer cache
println("#--- allocating benchmark ---#\n")
@time _, cache, _ = fmm!(system; lamb_helmholtz=false, scalar_potential=false, gradient=true, hessian=true)

# run FMM with preallocated buffers
println("\n#--- preallocated benchmark ---#\n")
@time fmm!(system; lamb_helmholtz=false, scalar_potential=false, gradient=true, hessian=true, cache...)
```
Note that the second call to `fmm!` allocates less memory.

### Tuning Parameters

The FMM can be tuned for accuracy and performance by adjusting the following parameters: `multipole_acceptance`, `leaf_size`, and `expansion_order`. These parameters are passed to the `fmm!` function as keyword arguments. Let's take a look at how each affects the FMM:

```@example guidedex
# create system
n_bodies, rand_seed = 1000, 123
system = generate_gravitational(rand_seed, n_bodies)

# compute potential directly
direct!(system; gradient=true)
gradient_direct = system.potential[5:7,:]

# try varying `expansion_order`
println("#--- varying expansion order ---#\n")
println("\texpansion_order = 3")
fmm!(system; expansion_order=1, multipole_acceptance=0.5, leaf_size=30, gradient=true)
system.potential .= 0.0
@time fmm!(system; expansion_order=3, multipole_acceptance=0.5, leaf_size=30, gradient=true)
println("\t\tmax error = ", maximum(abs.(system.potential[5:7,:] .- gradient_direct)))

println("\n\texpansion_order = 8")
system.potential .= 0.0
@time fmm!(system; expansion_order=5, multipole_acceptance=0.5, leaf_size=30, gradient=true)
println("\t\tmax error = ", maximum(abs.(system.potential[5:7,:] .- gradient_direct)))

println("\n\texpansion_order = 12")
system.potential .= 0.0
@time fmm!(system; expansion_order=10, multipole_acceptance=0.5, leaf_size=30, gradient=true)
println("\t\tmax error = ", maximum(abs.(system.potential[5:7,:] .- gradient_direct)))

# try varying `multipole_acceptance`
println("\n#--- varying multipole acceptance ---#\n")
println("\n\tmultipole_acceptance = 0.2")
system.potential .= 0.0
@time fmm!(system; expansion_order=4, multipole_acceptance=0.4, leaf_size=30, gradient=true)
println("\t\tmax error = ", maximum(abs.(system.potential[5:7,:] .- gradient_direct)))
println("\n\tmultipole_acceptance = 0.4")
system.potential .= 0.0
@time fmm!(system; expansion_order=4, multipole_acceptance=0.4, leaf_size=30, gradient=true)
println("\t\tmax error = ", maximum(abs.(system.potential[5:7,:] .- gradient_direct)))
println("\n\tmultipole_acceptance = 0.8")
system.potential .= 0.0
@time fmm!(system; expansion_order=4, multipole_acceptance=0.8, leaf_size=30, gradient=true)
println("\t\tmax error = ", maximum(abs.(system.potential[5:7,:] .- gradient_direct)))

# try varying `leaf_size`
println("\n#--- varying leaf size ---#\n")
println("\n\tleaf_size = 1")
system.potential .= 0.0
@time fmm!(system; expansion_order=4, multipole_acceptance=0.5, leaf_size=10, gradient=true)
println("\t\tmax error = ", maximum(abs.(system.potential[5:7,:] .- gradient_direct)))
println("\n\tleaf_size = 10")
system.potential .= 0.0
@time fmm!(system; expansion_order=4, multipole_acceptance=0.5, leaf_size=10, gradient=true)
println("\t\tmax error = ", maximum(abs.(system.potential[5:7,:] .- gradient_direct)))
println("\n\tleaf_size = 80")
system.potential .= 0.0
@time fmm!(system; expansion_order=4, multipole_acceptance=0.5, leaf_size=80, gradient=true)
println("\t\tmax error = ", maximum(abs.(system.potential[5:7,:] .- gradient_direct)))
```

This example demonstrates some trends of how each parameter affects the  The `expansion_order` parameter is always positively correlated to accuracy but negatively correlated to cost. The parameters `multipole_acceptance` and `leaf_size` are still correlated, but less predictably so. This motivates automated tuning of the parameters. We'll explore that more in the next example.

## Vortex Filament Example

In this example, we review the interface functions used by the vortex filament model found in `FastMultipole/test/vortex_filament.jl`. First, let's take a look at the data structure:

```@example guidedex
struct VortexFilaments{TF}
    x::Matrix{SVector{3,TF}}
    strength::Vector{SVector{3,TF}}
    core_size::Vector{TF}
    ε_tol::Vector{TF}
    potential::Vector{TF}
    force::Vector{SVector{3,TF}}
    gradient::Vector{SMatrix{3,3,TF,9}}
end
```
Note that `x` is a 2-row matrix, where the first row contains the start point of the filament and the second row contains the end point. The `strength` field contains the strength of the vortex filament, which is a vector quantity. The `core_size` field contains the regularization radius for each filament. The `potential`, `force`, and `gradient` fields are used to store the results of the FMM call.

### Overloading `body_to_multipole!`

Just like we did for the gravitational example, we need to overload `FastMultipole.body_to_multipole!` for our vortex filament system. The recommended way to do this is:

```@example guidedex
FastMultipole.body_to_multipole!(system::VortexFilaments, args...) = 
    FastMultipole.body_to_multipole!(Filament{Vortex}, system, args...)
```
Ordinarily, a 3-dimensional vector potential, such as the stream function in 3-dimensional fluid dynamics, requires 3 separate expansions (1 for each dimension). It is worth noting that when `Vortex` elements are used, `FastMultipole` never actually computes the vector potential ``\vec{\psi}``; rather, it computes two scalar components ``\phi`` and ``\chi`` of the Lamb-Helmholtz decomposition of the field as demonstrated by [gumerov2013efficient](@cite). This allows us to reduce the number of expansions required to express the scalar-plus-vector potential from 4 to 2, thus reducing computational cost. It is also worth noting that ``\phi`` and ``\chi`` are not coordinate system-invariant fields, and thus cannot be used to evaluate scalar or vector potentials. In other words, when a `Vortex` kernel is used as a source in a FMM call, the `scalar_potential` should not be used.

### Buffer Interface Functions

The buffer interface functions are defined for `VortexFilament` systems as follows:

```@example guidedex
function FastMultipole.source_system_to_buffer!(buffer, i_buffer, system::VortexFilaments, i_body)
    # position
    buffer[1:3,i_buffer] .= (system.x[1,i_body]  + system.x[2,i_body]) * 0.5
    
    # get regularization radius
    Γ = system.strength[i_body]
    # Γmag = norm(Γ)
    core_size = system.core_size[i_body]
    buffer[4,i_buffer] = system.core_size[i_body] + 0.5 * norm(system.x[2,i_body] - system.x[1,i_body])

    # remainding quantities
    buffer[5:7,i_buffer] .= Γ
    buffer[8:10,i_buffer] .= system.x[1,i_body]
    buffer[11:13,i_buffer] .= system.x[2,i_body]
    buffer[14,i_buffer] = core_size
end

Base.eltype(::VortexFilaments{TF}) where TF = TF

FastMultipole.data_per_body(::VortexFilaments) = 14

FastMultipole.get_position(system::VortexFilaments, i) = (system.x[1,i] + system.x[2,i]) * 0.5

FastMultipole.strength_dims(system::VortexFilaments) = 3

FastMultipole.get_n_bodies(system::VortexFilaments) = length(system.strength)

function FastMultipole.buffer_to_target_system!(target_system::VortexFilaments, i_target, ::DerivativesSwitch{PS,VS,GS}, target_buffer, i_buffer) where {PS,VS,GS}

    # extract from buffer
    PS && (potential = FastMultipole.get_scalar_potential(target_buffer, i_buffer))
    VS && (velocity = FastMultipole.get_gradient(target_buffer, i_buffer))
    GS && (hessian = FastMultipole.get_hessian(target_buffer, i_buffer))

    # load into system
    PS && (target_system.potential[i_target] += potential)
    VS && (target_system.force[i_target] += velocity)
    GS && (target_system.gradient[i_target] += hessian)

end
```
First, we point out that the `source_system_to_buffer!` function stores a core size in row 14. Although it will not be used by any internal `FastMultipole` functions, it will be required by the `direct!` function. Since the columns of the buffer are sorted into an octree structure, they will not match the order of the original system; therefore, it is important to store essential body properties in additional rows of the buffer. We also note that the `get_position` function needs not return a previously-stored value; in this case, it computes the midpoint of the filament on the fly.

### Overloading `direct!`

The `FastMultipole.direct!` function is overloaded for `VortexFilaments` systems as follows:

```@example guidedex
function get_δ(distance, core_size)
    δ = distance < core_size ? (distance-core_size) * (distance-core_size) : zero(distance)
    return δ
end

function vortex_filament_finite_core_2(x1,x2,xt,q,core_size)
    # intermediate values
    r1 = xt - x1
    r2 = xt - x2

    nr1 = norm(r1)
    nr2 = norm(r2)

    num = cross(r1, r2)
    denom = nr1 * nr2 + dot(r1, r2)

    # core size comes into play here
    distance_1 = norm((x1+x2)*0.5 - xt)
    δ1 = get_δ(distance_1, core_size)

    distance_2 = norm(x1 - xt)
    δ2 = get_δ(distance_2, core_size)

    distance_3 = norm(x2 - xt)
    δ3 = get_δ(distance_3, core_size)

    # desingularized terms
    f1 = num/(denom + δ1)
    f2 = 1/(nr1+δ2)
    f3 = 1/(nr2+δ3)

    # evaluate velocity
    V = (f1*(f2+f3))/(4*pi) * q

    return V
end

function FastMultipole.direct!(target_system, target_index, derivatives_switch::DerivativesSwitch{PS,VS,GS}, source_system::VortexFilaments, source_buffer, source_index) where {PS,VS,GS}
    for i_source in source_index
        x1 = FastMultipole.get_vertex(source_buffer, source_system, i_source, 1)
        x2 = FastMultipole.get_vertex(source_buffer, source_system, i_source, 2)
        q = FastMultipole.get_strength(source_buffer, source_system, i_source)
        core_size = source_buffer[14, i_source]

        for i_target in target_index
            xt = FastMultipole.get_position(target_system, i_target)

            if VS
                # determine sign of q
                q_mag = norm(q) * sign(dot(q, x2-x1))

                # calculate velocity
                v = vortex_filament_finite_core_2(x1,x2,xt,q_mag,core_size)
                FastMultipole.set_gradient!(target_system, i_target, v)
            end
        end
    end
end
```
Note that `core_size` is retrieved from the buffer matrix rather than from `source_system` directly, as discussed previously.

### Running the FMM

Let's try this out on a vortex filament system. The `fmm!` function is called in the same way as for the gravitational example, but with the `VortexFilaments` system:

```@example guidedex
using LinearAlgebra

# create system
function generate_filament_field(n_filaments, length_scale; strength_scale=1/n_filaments)
    centers = rand(SVector{3,Float64}, n_filaments)
    pts = zeros(SVector{3,Float64}, 2, n_filaments)
    strength_vec= zeros(SVector{3,Float64}, n_filaments)
    for (i,center) in enumerate(centers)
        dx = (rand(SVector{3,Float64}) * 2 .- 1.0) * 0.5 * length_scale
        pts[1,i] = center - dx
        pts[2,i] = center + dx
        Γ = dx / norm(dx) * rand() * strength_scale
        strength_vec[i] = Γ
    end

    # create filaments
    core_size = fill(1e-2, n_filaments)
    ε_tol = fill(1e-4, n_filaments)
    potential = zeros(length(strength_vec))
    gradient = zeros(SVector{3,Float64}, length(strength_vec))
    hessian = zeros(SMatrix{3,3,Float64,9}, length(strength_vec))

    return VortexFilaments(pts, strength_vec, core_size, ε_tol, potential, gradient, hessian)
end

n_filaments = 1000
length_scale = 1.0 / n_filaments^(1/3)
filaments = generate_filament_field(n_filaments, length_scale; strength_scale=1/n_filaments)

# run FMM
_, cache, _ = fmm!(filaments; lamb_helmholtz=true, scalar_potential=false, gradient=true)
```
We set `lamb_helmholtz=true` because the vortex filament model induces a vector potential.

### Optimal Leaf Size

We can request `FastMultipole` to predict the optimal leaf size for our system by setting `tune=true` in the `fmm!` call. This performs best when `interaction_list_method=SelfTuning()` on a single thread, but still functions well for other choices. If `tune=true`, benchmarks will be used to estimate the leaf size at which the cost of direct calculations is approximately equal to expansion calculations. It is returned as part of a named tuple as the first returned value of `fmm!`, which can be splatted as a keyword argument in subsequent `fmm!` calls:


```@example guidedex
# get optimal leaf size
println("#--- default leaf size ---#\n")
@time optargs, _ = fmm!(filaments; tune=true, lamb_helmholtz=true, scalar_potential=false, gradient=true, cache...)

# run again with optimal leaf size
println("\n#--- optimal leaf size ---#\n")
@time fmm!(filaments; lamb_helmholtz=true, scalar_potential=false, gradient=true, cache..., optargs...)
```
Note that `optargs` contains `leaf_size_source`, `expansion_order`, and `multipole_acceptance` parameters. Only `leaf_size_source` is tuned if `isnothing(ε_tol) == true`. More complete auto-tuning will be discussed later.
