# Overview

This package contains routines to initialize reference element operators, physical mesh arrays, and connectivity arrays for nodal DG methods. The codes roughly based on *Nodal Discontinuous Galerkin Methods* by Hesthaven and Warburton (2007).

StartUpDG.jl is intended mainly to aid in the implementation of `rhs!` functions for DG discretizations of time-dependent partial differential equations, which can then be used with the OrdinaryDiffEq.jl library to evolve a solution in time. For example, it has been used in most publications since 2020 by the authors [Jesse Chan](https://scholar.google.com/citations?user=rqGSShYAAAAJ&hl=en) and [Yimin Lin](https://scholar.google.com/citations?hl=en&user=yCrSttgAAAAJ), as well as in the following external publications: 
* [Efficient entropy-stable discontinuous spectral-element methods using tensor-product summation-by-parts operators on triangles and tetrahedra](https://doi.org/10.1016/j.jcp.2024.113360) by Montoya and Zingg (2024).
* [Injected Dirichlet boundary conditions for general diagonal-norm SBP operators](https://www.researchgate.net/profile/Anita-Gjesteland/publication/374234334_Injected_Dirichlet_boundary_conditions_for_general_diagonal-norm_SBP_operators/links/6515500dcce2460b6c3d6eda/Injected-Dirichlet-boundary-conditions-for-general-diagonal-norm-SBP-operators.pdf) by Gjesteland, Del Rey Fernández, and Svärd (2023).

It is also used in the [Trixi.jl](https://github.com/trixi-framework/Trixi.jl/) library.

## A short example

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
(; x, y ) = md
u = @. exp(-10 * (x^2 + y^2))

# Compute derivatives using geometric mapping + chain rule
(; Dr, Ds ) = rd
(; rxJ, sxJ, J ) = md
dudx = (rxJ .* (Dr * u) + sxJ .* (Ds * u)) ./ J
```
