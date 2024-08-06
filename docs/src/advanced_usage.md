
# Advanced Usage

Though the typical user won't need to alter these parameters, various situations may require modifications to the evaluation process of the FMM.

## Unsort Bodies

As part of the evaluation, the `fmm!` call will sort bodies contained in the systems in order to optimize performance, thus during the fmm evaluation, the body indices are different from initial indices given by the user. At the conclusion of the `fmm!` call, these bodies are then put back in their original order for the convenience of the user. This feature can be toggled off so that the bodies are left in their sorted order after the FMM evaluation.

```@example advanced
import FastMultipole as fmm
using Random
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

target_system = generate_gravitational(123, 100)
source_system = generate_gravitational(321, 100)

fmm.fmm!(target_system, source_system,
    unsort_source_bodies=true, unsort_target_bodies=true) # standard run (defaults)
fmm.fmm!(target_system, source_system,
    unsort_source_bodies=false, unsort_target_bodies=false) # sorted bodies (not original indices)
```

## Reusing the Octree

Users are able to define the tree to be used in the FMM by passing in a previously created tree into the `fmm!` call. For convenience, the `upward_pass`, `horizontal_pass`, and `downward_pass` can be toggled off seperately. These are important to toggle when generating octrees to be reused. Additionally, `unsort_target_bodies` and `unsort_source_bodies` should be set to `false` when generating the reusable trees.

- **`upward_pass` creates the multipole expansions for each branch**
- **`horizontal_pass` evaluates the expansions in the appropriate locations**
- **`downward_pass` applies parent expansions to children branches**

```@example advanced
target_system = generate_gravitational(123, 100)
source_system = generate_gravitational(321, 100)

target_tree, source_tree = fmm.fmm!(target_system, source_system;
    upward_pass=true, horizontal_pass=false, downward_pass=false,
    unsort_source_bodies=false, unsort_target_bodies=false) # generating trees to be reused

fmm.fmm!(target_tree, target_system, source_tree, source_system;
    upward_pass=false, horizontal_pass=true, downward_pass=true) # reusing trees
```

## Partial Evaluation

The `nearfield`, `farfield`, and `self_induced` parameters allow users to bypass certain portions of the FMM.

- **`nearfield` if false, `FastMultipole` omits all direct interactions**
- **`farfield` if false, `FastMultipole` omits all multipole interactions**
- **`self_induced` if false, `FastMultipole` omits the interaction of each leaf branch on itself**

```@example advanced
system = generate_gravitational(123, 100)
fmm.fmm!(system,
    nearfield=true, farfield=true, self_induced=true) # standard run (defaults)
fmm.fmm!(system,
    nearfield=false, farfield=true, self_induced=true) # nearfield not evaluated
```

## Resize and Recenter Branches

The initial creation of the octree will create uniformly sized and spaced branches in each level. These branches are not centered or optimally sized with relation to the bodies they contain. As a default, the `fmm!` call resizes and recenters these branches to utilize the smallest possible radius that contains all the bodies in each branch[DENG2021265](@cite). This feature can be toggled with the `source_shrink_recenter` and `target_shrink_recenter` parameters.

```@example advanced
target_system = generate_gravitational(123, 100)
source_system = generate_gravitational(321, 100)

fmm.fmm!(target_system, source_system,
    source_shrink_recenter=true, target_shrink_recenter=true) # standard run (defaults)
fmm.fmm!(target_system, source_system,
    source_shrink_recenter=false, target_shrink_recenter=false) # no branch recentering
```

## References

```@bibliography
```
