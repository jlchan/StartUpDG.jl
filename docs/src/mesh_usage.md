# MeshData

The `MeshData` struct contains data for high order DG methods useful for evaluating DG formulations in a matrix-free fashion.

## Creating unstructured meshes

For convenience, simple uniform meshes are included in with `StartUpDG.jl`.
```julia
using StartUpDG
Kx,Ky,Kz = 4,2,8
VX,EToV = uniform_mesh(Line(),Kx)
VX,VY,EToV = uniform_mesh(Tri(),Kx,Ky)
VX,VY,EToV = uniform_mesh(Quad(),Kx,Ky)
VX,VY,VZ,EToV = uniform_mesh(Hex(),Kx,Ky,Kz)
```
The triangular mesh is constructed by creating a uniform quadrilateral mesh then bisecting each quad into two triangles.

Unstructured meshes for more complex geometries can be generated using external packages. For example, `TriangleMesh.jl` can be used as follows:
```julia
using TriangleMesh

poly = polygon_Lshape()
mesh = create_mesh(poly, set_area_max=true) # asks for maximum element size
VX,VY = mesh.point[1,:],mesh.point[2,:]
EToV = permutedims(mesh.cell)
```

## Initializing high order DG mesh data

Given unstructured mesh information (tuple of vertex coordinates `VXYZ` and index array `EToV`) high order DG mesh data can be constructed as follows:
```julia
md = MeshData(VXYZ...,EToV,rd)
```

## Enforcing periodic boundary conditions

```julia
md = MeshData(VX,VY,EToV,rd)
md_periodic = make_periodic(md,rd) # periodic in both x and y coordinates
md_periodic_x = make_periodic(md,rd,true,false) # periodic in x direction, but not y
```

## Curved meshes

```julia
md = MeshData(VX,VY,EToV,rd)
@unpack x,y = md
# <-- code to modify x,y
md_curved = MeshData(md,rd,x,y)
```
