# FastMultipole.jl

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://flow.byu.edu/FastMultipole)

*A fast, multi-system, multi-kernel, differentiable implementation of the fast multipole method for use with scalar-plus-vector potential N-body problems in pure Julia.*

Authors: Ryan Anderson and contributors

Features:

* Laplace $1/r$ kernel
* meant to incorporate seamlessly into existing Julia code without modifications (only additions to define a few interface functions)
* convenience functions for determining the expansion coefficients for source points, vortex points, source panels, and dipole panels (this list is growing!)
* $\mathcal{O}(p^4)$ multipole-to-local translation operator (where ```math p``` is the expansion order), though this may improve in the near future
* velocity and velocity gradient (or their equivalent for non-fluids problems) obtained using analytic expressions (no finite difference)
* CPU-parallelization for expansions and direct interactions
* GPU-parallelization for direct interactions using `CUDA.jl`
* `ForwardDiff` and `ReverseDiff` compatible

