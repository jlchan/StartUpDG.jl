# Changelog

StartUpDG.jl follows the interpretation of [semantic versioning (semver)](https://julialang.github.io/Pkg.jl/dev/compatibility/#Version-specifier-format-1) used in the Julia ecosystem. Recent changes will be documented in this file for human readability.

## Changes when updating to v0.16

#### Added
* Add `Polynomial{Gauss}` type

#### Changed
* add type parameter to `Polynomial`
* Switch from `LittleDict` to `MultipleRefElemData` for hybrid meshes

## Changes when updating to v0.15

#### Added 

* the `NamedArrayPartition` array type, which is similar to `ComponentArrays` but with the storage structure of `ArrayPartition`. This is used for the storage of `MeshData` fields in hybrid and cut-cell meshes, and can be used for the storage of solutions compatible with the OrdinaryDiffEq.jl framework. 
* added precomputed differentiation, face interpolation, mass, and lifting matrices for `CutCellMesh` types. These are specified using `md = MeshData(...; precompute_operators=true)`, and are stored in `md.mesh_type.cut_cell_operators`. 

#### Changed

* The `MeshData` fields `VXYZ` and `EToV` have been moved into a `VertexMappedMesh` type. However, they can still be accessed as a property of `MeshData` (e.g., `md.VX` will still work). 
* Added a `CurvedMesh` type for `MeshData` constructed from modified nodal coordinates, e.g., using the constructor `MeshData(rd, md, xyz...)`. `CurvedMesh` stores the original mesh type as a field. Previously, there was no way to distinguish a curved `MeshData` from a non-curved one.
* Changes related to element types:
  * upstream change in NodesAndModes.jl: the abstract type `AbstractElemShape` is now parametrized by the dimension, e.g., `Line <: AbstractElemShape{1}`, `Tri <: AbstractElemShape{2}`. 
  * the spatial dimension `Dim` parameter in `MeshData` is now inferred from the element type through `elem <: AbstractElemShape{Dim}`
  * the `PhysicalFrame` type now has leading type parameter `PhysicalFrame{NDIMS} <: AbstractElemShape{NDIMS}`
* `PhysicalFrame`'s fields are now restricted to be type `SVector` only (instead of `Union{SVector, Tuple}`)

#### Removed 

* the `Nplot` field has been removed from `RefElemData`
* all usages of `ComponentArrays` have been replaced by `NamedArrayPartition`
* the deprecated `MeshData(md::MeshData, rd::RefElemData, xyz...)` constructor has been removed

#### Deprecated

* the constructor for `MeshData` with all fields as well as `VXYZ`, `EToV` as arguments has been deprecated and will be removed in the next breaking release. This would only be used when extending `MeshData`. The standard `MeshData` constructor involving `VXYZ`, `EToV`, and `rd::RefElemData` is unchanged. 
* the constructor for `RefElemData` with all fields as well as `Nplot` has been deprecated and will be removed in the next breaking release. This is only used if extending `RefElemData`. 

