# FastMultipole

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://flow.byu.edu/FastMultipole/dev)
![](https://github.com/byuflowlab/FastMultipole/actions/workflows/CI.yml/badge.svg)

*A fast, multi-system, multi-kernel, differentiable implementation of the fast multipole method for use with scalar-plus-vector potential N-body problems in pure Julia.*

**Author**: Ryan Anderson

**Features**:

* solves $N$-body problems governed by the Laplace ($1/r$) Green's function, with work planned to support the Helmholtz Green's function in the future
* incorporates seamlessly into existing Julia code without modifications (just the addition of a few interface functions)
* offers convenience functions for determining the expansion coefficients for multiple kernel functions including source points, vortex points, source panels, and dipole panels (this list is growing!)
* provides velocity and velocity gradient (or their equivalent for non-fluids problems) obtained using analytic expressions (no finite difference or AD)
* uses $\mathcal{O}(p^4)$ multipole-to-local translation operator (where $p$ is the expansion order) with more efficient methods under development
* automated CPU-parallelization of expansions and direct interactions
* supports GPU-parallelization of direct interactions using [CUDA](https://github.com/JuliaGPU/CUDA.jl)
* [ForwardDiff](https://github.com/JuliaDiff/ForwardDiff.jl) and [ReverseDiff](https://github.com/JuliaDiff/ReverseDiff.jl) compatible

**Documentation**

See the [docs](https://flow.byu.edu/FastMultipole/dev).
