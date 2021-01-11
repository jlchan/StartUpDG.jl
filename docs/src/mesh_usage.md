# MeshData

## Creating a simple uniform mesh

For convenience, simple uniform meshes are included in the

## Generating meshes

Unstructured meshes for more complex geometries can be generated using external packages. For example, `TriangleMesh.jl`

```julia
using TriangleMesh

poly = polygon_Lshape()
mesh = create_mesh(poly, set_area_max=true) # asks for maximum element size

VX,VY = mesh.point[1,:],mesh.point[2,:]
EToV = permutedims(mesh.cell)
```

## Initializing high order DG mesh data

## Enforcing periodic boundary conditions

## Curved meshes
