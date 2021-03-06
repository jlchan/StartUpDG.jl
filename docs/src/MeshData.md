# `MeshData` type

[`MeshData`](@ref) contains the following fields
* `K` or `num_elements`: number of elements in the mesh 
* `FToF`: indexing vector for face-to-face connectivity (length of the vector is the total number of faces, e.g., ``N_{\rm faces} K``)
* `xyz::NTuple{Dim, ...}`: nodal interpolation points mapped to physical elements. All elements of `xyz` are ``N_p \times K`` matrices, where ``N_p`` are the number of nodal points on each element.
* `xyzq::NTuple{Dim, ...}, wJq`: volume quadrature points/weights mapped to physical elements. All elements these tuples are ``N_q \times K`` matrices, where ``N_q`` is the number of quadrature points on each element.
* `xyzf::NTuple{Dim, ...}`: face quadrature points mapped to physical elements. All elements of `xyz` are ``N_f \times K`` matrices, where ``N_f`` is the number of face points on each element.
* `mapP, mapB`: indexing arrays for inter-element node connectivity (`mapP`) and for extracting boundary nodes from the list of face nodes `xyzf` (`mapB`). `mapP` is a matrix of size ``N_f \times K``, while the length of `mapB` is the total number of nodes on the boundary.
* `rstxyzJ::SMatrix{Dim, Dim}`: volume geometric terms ``G_{ij} = \frac{\partial x_i}{\partial \hat{x}_j}``. Each element of `rstxyzJ` is a matrix of size ``N_p \times K``.
* `J, sJ`: volume and surface Jacobians evaluated at interpolation points and surface quadrature points, respectively. `J` is a matrix of size ``N_p \times K``, while `sJ` is a matrix of size ``N_f \times K``. 
* `nxyzJ::NTuple{Dim, ...}`: scaled outward normals evaluated at surface quadrature points. Each element of `nxyzJ` is a matrix of size ``N_f\times K``. 

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
md_periodic = make_periodic(md) # periodic in both x and y coordinates
md_periodic_x = make_periodic(md,true,false) # periodic in x direction, but not y
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
md_curved = MeshData(rd,md,x,y)
```
`MeshData(rd,md,x,y)` and `MeshData(rd,md,x,y,z)` are implemented for 2D and 3D, though this is not currently implemented in 1D.

More generally, one can create a copy of a `MeshData` with certain fields modified by using `@set` or `setproperties` from `Setfield.jl`.

## Unstructured triangular meshes using Triangulate

If `Triangulate` is also loaded, then StartUpDG will includes a few additional utilities for creating and visualizing meshes. 

## Tagging boundary faces

One can "tag" boundary faces by specifying boolean functions which evaluate to `true` if a point is on a given boundary segment. 
```julia
using Test

rd = RefElemData(Tri(),N=3)
md = MeshData(uniform_mesh(Tri(),1)...,rd)
on_bottom_boundary(x,y,tol=1e-13) = abs(y+1) < tol
on_top_boundary(x,y,tol=1e-13) = abs(y-1) < tol

boundary_dict = tag_boundary_faces(md, Dict(:bottom=>on_bottom_boundary,:top=>on_top_boundary))
@test boundary_dict == Dict(:bottom=>[1],:top=>[4])
```

You can also specify a list of boundaries using NamedTuples 
```julia
boundary_dict = tag_boundary_faces(md,(; :bottom=>on_bottom_boundary,:top=>on_top_boundary))
```