# StartUpDG
[![Docs-stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jlchan.github.io/StartUpDG.jl/stable)
[![Docs-dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jlchan.github.io/StartUpDG.jl/dev)
[![Build status](https://github.com/jlchan/StartUpDG.jl/workflows/CI/badge.svg)](https://github.com/jlchan/StartUpDG.jl/actions)
[![Codecov](https://codecov.io/gh/jlchan/StartUpDG.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/jlchan/StartUpDG.jl)

This package contains routines to initialize reference element operators, physical mesh arrays, and connectivity arrays for nodal discontinuous Galerkin (DG) methods. It is intended to quick-start the implementation of from-scratch high order time-domain nodal DG solvers for partial differential equations (PDEs), and is designed for flexibility and experimentation. 

The codes are roughly based on *Nodal Discontinuous Galerkin Methods* by Hesthaven and Warburton (2007). The original port from Matlab to Julia was by [Yimin Lin](https://github.com/yiminllin), with subsequent modifications by Jesse Chan and other contributors. 

This package is registered and can be installed via `] add StartUpDG` or `using Pkg; Pkg.add("StartUpDG")`.

# Usage overview

Variables are contained within structs `rd::RefElemData` and `md::MeshData`, which contain quantities from `Globals1D, Globals2D, Globals3D` in the Nodal DG book codes. These can be used to compute DG derivatives, and are useful for matrix-free implementations of DG methods using explicit time-stepping.

```julia
using StartUpDG

# polynomial degree and mesh size
N = 3
K1D = 8

# init ref element and mesh
rd = RefElemData(Tri(), N)
VXY, EToV = uniform_mesh(Tri(), K1D)
md = MeshData(VXY, EToV, rd)

# Define a function by interpolation
(; x, y) = md
u = @. 2 + 0.5 * exp(-100 * (x^2 + y^2))

# Compute derivatives using geometric mapping + chain rule
(; Dr, Ds) = rd
(; rxJ, sxJ, J) = md
dudx = (rxJ .* (Dr * u) + sxJ .* (Ds * u)) ./ J
```

# Contributors

* SBP nodal points were contributed by [Ethan Kubatko](https://sites.google.com/site/chilatosu/ethan-bio) and [Jason Hicken](https://doi.org/10.1007/s10915-020-01154-8). 
* [Hendrik Ranocha](https://ranocha.de) contributed to array types used in cut-cell and hybrid meshes. 
* [Mason McCallum](https://github.com/masonamccallum) contributed Gmsh reading capabilities
* [David Knapp](https://github.com/Davknapp) contributed VTK visualization capabilities and tensor product wedges.
