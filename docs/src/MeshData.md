# `MeshData` type

[`MeshData`](@ref) contains the following fields
* `K`: number of elements ``K`` in the mesh.
* `FToF`: indexing vector for face-to-face connectivity (length of the vector is the total number of faces, e.g., ``N_{\rm faces} K``)
* `xyz::NTuple{Dim,...}`: nodal interpolation points mapped to physical elements. All elements of `xyz` are ``N_p \times K`` matrices, where ``N_p`` are the number of nodal points on each element.
* `xyzq::NTuple{Dim,...}, wJq`: volume quadrature points/weights mapped to physical elements. All elements these tuples are ``N_q \times K`` matrices, where ``N_q`` is the number of quadrature points on each element.
* `xyzf::NTuple{Dim,...}`: face quadrature points mapped to physical elements. All elements of `xyz` are ``N_f \times K`` matrices, where ``N_f`` is the number of face points on each element.
* `mapP,mapB`: indexing arrays for inter-element node connectivity (`mapP`) and for extracting boundary nodes from the list of face nodes `xyzf` (`mapB`). `mapP` is a matrix of size ``N_f \times K``, while the length of `mapB` is the total number of nodes on the boundary.
* `rstxyzJ::SMatrix{Dim,Dim}`: volume geometric terms ``G_{ij} = \frac{\partal x_i}{\partial \hat{x}_j}``. Each element of `rstxyzJ` is a matrix of size ``N_p \times K``.
* `J,sJ`: volume and surface Jacobians evaluated at interpolation points and surface quadrature points, respectively. `J` is a matrix of size ``N_p \times K``, while `sJ` is a matrix of size ``N_f \times K``.
* `nxyzJ::NTuple{Dim,...}`: scaled outward normals evaluated at surface quadrature points. Each element of `nxyzJ` is a matrix of size ``N_f\times K``.

# Setting up `md::MeshData`

The [`MeshData`](@ref) struct contains data for high order DG methods useful for evaluating DG formulations in a matrix-free fashion.

## Generating unstructured meshes

For convenience, simple uniform meshes are included in with `StartUpDG.jl` via [`uniform_mesh`](@ref)
```julia
using StartUpDG
Kx,Ky,Kz = 4,2,8
VX,EToV = uniform_mesh(Line(),Kx)
VX,VY,EToV = uniform_mesh(Tri(),Kx,Ky)
VX,VY,EToV = uniform_mesh(Quad(),Kx,Ky)
VX,VY,VZ,EToV = uniform_mesh(Hex(),Kx,Ky,Kz)
```
The uniform triangular mesh is constructed by creating a uniform quadrilateral mesh then bisecting each quad into two triangles.

## Initializing high order DG mesh data

Given unstructured mesh information (tuple of vertex coordinates `VXYZ` and index array `EToV`) high order DG mesh data can be constructed as follows:
```julia
md = MeshData(VXYZ...,EToV,rd)
```

## Enforcing periodic boundary conditions

Periodic boundary conditions can be enforced by calling [`make_periodic`](@ref), which returns another `MeshData` struct with modified `mapP`,`mapB`, and `FToF` indexing arrays which account for periodicity.
```julia
md = MeshData(VX,VY,EToV,rd)
md_periodic = make_periodic(md,rd) # periodic in both x and y coordinates
md_periodic_x = make_periodic(md,rd,true,false) # periodic in x direction, but not y
```
One can check which dimensions are periodic via the `is_periodic` field of `MeshData`. For example, the `md_periodic_x` example above gives
```julia
julia> md_periodic_x.is_periodic
(true, false)
```

## Creating curved meshes

It's common to generate curved meshes by first generating a linear mesh, then moving high order nodes on the linear mesh. This can be done by calling [`MeshData`](@ref) again with new `x,y` coordinates:
```julia
md = MeshData(VX,VY,EToV,rd)
@unpack x,y = md
# <-- code to modify high order nodes (x,y)
md_curved = MeshData(md,rd,x,y)
```
`MeshData(md,rd,x,y)` and `MeshData(md,rd,x,y,z)` are implemented for 2D and 3D, though this is not currently implemented in 1D.

More generally, one can create a copy of a `MeshData` with certain fields modified by using `@set` or `setproperties` from `Setfield.jl`.

## Unstructured triangular meshes using Triangulate

If `Triangulate` is also loaded, then StartUpDG will includes a few additional utilities for creating and visualizing meshes. 
