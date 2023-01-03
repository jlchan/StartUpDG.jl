# Changelog

StartUpDG.jl follows the interpretation of [semantic versioning (semver)](https://julialang.github.io/Pkg.jl/dev/compatibility/#Version-specifier-format-1) used in the Julia ecosystem. Recent changes will be documented in this file for human readability.

## Changes when updating to v0.15

#### Added 

* the `NamedArrayPartition` array type, which is similar to `ComponentArrays` but with the storage structure of `ArrayPartition`. This is used for the storage of `MeshData` fields in hybrid and cut-cell meshes, and can be used for the storage of solutions compatible with the OrdinaryDiffEq.jl framework. 
* Added a `CurvedMesh` type for `MeshData` constructed from modified nodal coordinates, e.g., using the constructor `MeshData(rd, md, xyz...)`. `CurvedMesh` stores the original mesh type as a field. 

#### Changed

* Changes related to element types and `MeshData`:
  * upstream change in NodesAndModes.jl: the abstract type `AbstractElemShape` is now parametrized by the dimension, e.g., `Line <: AbstractElemShape{1}`, `Tri <: AbstractElemShape{2}`. 
  * the spatial dimension `Dim` parameter in `MeshData` is now inferred from the element type through `elem <: AbstractElemShape{Dim}`
  * the `PhysicalFrame` type now has leading type parameter `PhysicalFrame{NDIMS} <: AbstractElemShape{NDIMS}`
* the `Nplot` field has been removed from `RefElemData`
* `PhysicalFrame`'s fields are now restricted to be type `SVector` only (instead of `Union{SVector, Tuple}`)
* all usages of `ComponentArrays` have been replaced by `NamedArrayPartition`

#### Deprecated

* the old constructor for `RefElemData` with field `Nplot` has been deprecated and will be removed in the next breaking release. 

