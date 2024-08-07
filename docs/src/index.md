# FastMultipole

*A fast, multi-system, multi-kernel, differentiable implementation of the fast multipole method for use with scalar-plus-vector potential N-body problems in pure Julia.*

**Author**: Ryan Anderson

**Features**:

* solves ``N``-body problems governed by the Laplace (``1/r``) kernel, with work planned to support the Helmholtz kernel in the future
* incorporates seamlessly into existing Julia code without modifications (just the addition of a few interface functions)
* offers convenience functions for determining the expansion coefficients for source points, vortex points, source panels, and dipole panels (this list is growing!)
* provides velocity and velocity gradient (or their equivalent for non-fluids problems) obtained using analytic expressions (no finite difference)
* uses ``\\mathcal{O}(p^4)`` multipole-to-local translation operator (where ``p`` is the expansion order), though this may improve in the near future
* automated CPU-parallelization of expansions and direct interactions
* supports GPU-parallelization of direct interactions using [CUDA](https://github.com/JuliaGPU/CUDA.jl)
* [ForwardDiff](https://github.com/JuliaDiff/ForwardDiff.jl) and [ReverseDiff](https://github.com/JuliaDiff/ReverseDiff.jl) compatible

**Installation**:

```julia
pkg> add https://github.com/byuflowlab/FastMultipole.git
```

**Documentation**:

* learn basic useage in the [Quick Start](quickstart.md) tutorial
* discover more features in the [Guided Examples](guided_examples.md)
* fine-tune performance in the [Advanced Usage](advanced_usage.md) section (for FMM experts)
* see the full [API](reference.md)
* brush up on the [Theory](theory.md)

