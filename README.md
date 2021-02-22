# StartUpDG
[![Docs-stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jlchan.github.io/StartUpDG.jl/stable)
[![Docs-dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jlchan.github.io/StartUpDG.jl/dev)
[![Build status](https://github.com/jlchan/StartUpDG.jl/workflows/CI/badge.svg)](https://github.com/jlchan/StartUpDG.jl/actions)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/jlchan/StartUpDG.jl?svg=true)](https://ci.appveyor.com/project/jlchan/StartUpDG-jl)
[![Codecov](https://codecov.io/gh/jlchan/StartUpDG.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jlchan/StartUpDG.jl)

Routines to initialize reference element operators, physical mesh arrays, and connectivity for nodal discontinuous Galerkin (DG) methods. Codes roughly based on *Nodal Discontinuous Galerkin Methods* by Hesthaven and Warburton (2007). Original port from Matlab to Julia by [Yimin Lin](https://github.com/yiminllin).

This is not yet a registered package.

# Usage overview

Variables are contained within structs `rd::RefElemData` and `md::MeshData`. These structs contain variables from `Globals1D, Globals2D, Globals3D` in the Nodal DG book codes. Variables can be unpacked using [`@unpack`](https://github.com/mauro3/UnPack.jl).

```julia
using StartUpDG

# polynomial degree and mesh size
N = 3
K1D = 8

# init ref element and mesh
rd = RefElemData(Tri(),N)
VX,VY,EToV = uniform_mesh(Tri(),K1D)
md = MeshData(VX,VY,EToV,rd)

# Define a function by interpolation
@unpack x,y = md
u = @. 2 + .5*exp(-100*(x^2+y^2))

# Compute derivatives using geometric mapping + chain rule
@unpack Dr,Ds = rd
@unpack rxJ,sxJ,J = md
dudx = (rxJ.*(Dr*u) + sxJ.*(Ds*u))./J
```
