# Overview

This package contains routines to initialize reference element operators, physical mesh arrays, and connectivities for nodal DG methods. Codes roughly based on *Nodal Discontinuous Galerkin Methods* by Hesthaven and Warburton (2007).

`StartUpDG.jl` exports structs `RefElemData{Dim,ElemShape,...}` and `MeshData{Dim}` which contain data useful for evaluating DG formulations in a matrix-free fashion. These structs contain fields similar to those in `Globals1D, Globals2D, Globals3D` in the NDG book codes.

Fields can be unpacked from using [`@unpack`](https://github.com/mauro3/UnPack.jl).

## A short example

```julia
using StartUpDG
using UnPack

# polynomial degree and mesh size
N = 3
K1D = 8

# init ref element and mesh
rd = RefElemData(Tri(),N)
VX,VY,EToV = uniform_mesh(Tri(),K1D)
md = MeshData(VX,VY,EToV,rd)

# Define a function by interpolation
@unpack x,y = md
u = @. exp(-10*(x^2+y^2))

# Compute derivatives using geometric mapping + chain rule
@unpack Dr,Ds = rd
@unpack rxJ,sxJ,J = md
dudx = (rxJ.*(Dr*u) + sxJ.*(Ds*u))./J
```
