# StartUpDG
[![Build Status](https://travis-ci.com/jlchan/StartUpDG.jl.svg?branch=master)](https://travis-ci.com/jlchan/StartUpDG.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/jlchan/StartUpDG.jl?svg=true)](https://ci.appveyor.com/project/jlchan/StartUpDG-jl)
[![Build status](https://github.com/jlchan/StartUpDG.jl/workflows/CI/badge.svg)](https://github.com/jlchan/StartUpDG.jl/actions)
[![Codecov](https://codecov.io/gh/jlchan/StartUpDG.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jlchan/StartUpDG.jl)

Routines to initialize reference element operators, physical mesh arrays, and connectivity for nodal DG methods

Codes roughly based on Nodal Discontinuous Galerkin (NDG) methods by Hesthaven and Warburton (2007). Original port from Matlab to Julia by [Yimin Lin](https://github.com/yiminllin)

# Usage

Variables are contained within (mutable) convenience structs `rd::RefElemData` and `md::MeshData`. These structs contain variables from `Globals1D, Globals2D, Globals3D` in the NDG book codes. Variables can be unpacked from (and repacked into) `rd::RefElemData` and `md::MeshData` using [`@unpack`](https://github.com/mauro3/UnPack.jl).

```
using StartUpDG

# polynomial degree and mesh size
N = 3
K1D = 8

# init ref element and mesh
rd = init_reference_tri(N)
VX,VY,EToV = uniform_tri_mesh(K1D)
md = init_DG_mesh(VX,VY,EToV,rd)

# Define a function by interpolation
@unpack x,y = md
u = @. 2 + .5*exp(-100*(x^2+y^2))

# Compute derivatives using geometric mapping + chain rule
@unpack Dr,Ds = rd
@unpack rxJ,sxJ,J = md
dudx = (rxJ.*(Dr*u) + sxJ.*(Ds*u))./J
```
