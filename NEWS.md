# Changelog

StartUpDG.jl follows the interpretation of [semantic versioning (semver)](https://julialang.github.io/Pkg.jl/dev/compatibility/#Version-specifier-format-1) used in the Julia ecosystem. Recent changes will be documented in this file for human readability.

## Changes when updating to v1.0.0

Most of the major changes are tracked in this [PR](https://github.com/jlchan/StartUpDG.jl/pull/160). Some descriptions of other changes are listed below. 

#### Added

* Generation of cut-cell meshes using `Subtriangulation` quadrature by default, which ensures positive quadrature weights. The old behavior is retained by specifying a `MomentFitting` quadrature type. 
* Added `subcell_limiting_operators`, which constructs multi-dimensional versions of subcell operators used in [Lin, Chan (2024)](https://doi.org/10.1016/j.jcp.2023.112677). These subcell operators are constructed from sparse operators returned by `sparse_low_order_SBP_operators(rd::RefElemData)`.

#### Changed

* `NamedArrayPartition` was moved to RecursiveArrayTools.jl. 
* The required Julia version was increased to v1.10. This was to make StartUpDG.jl compatibility with RecursiveArrayTools.jl v3.4+ (see above).
* Removed SimpleUnpack.jl as a dependency. Loading StartUpDG.jl will no longer reexport `@unpack`, since destructuring via `(; propertyname) = x` is supported natively in Julia 1.7 and up.
* Updated to NodesAndModes v1.0+, which changed the ordering of triangle nodes to make them consistent with tet node ordering. 
* We introduced a `MultidimensionalQuadrature` type. All `Polynomial` approximation types now utilize either `MultidimensionalQuadrature` or `TensorProductQuadrature` as a type parameter. The previous type parameter `DefaultPolynomialType` is now simply used to determine the default quadrature type parameter. Note that this is internal behavior and should not impact standard usage of StartUpDG.jl.
* Removed Requires.jl in favor of [package extensions](https://pkgdocs.julialang.org/v1/creating-packages/#Conditional-loading-of-code-in-packages-(Extensions)) for Plots.jl and SummationByPartsOperators.jl dependencies. 


## Changes when updating to v0.17

#### Added
* Added support for pyramid reference elements and meshes
* Introduced unified `read_Gmsh_2D(...)` functions 
* Specializations of `RefElemData` for `Polynomial(TensorProductQuadrature(quad_rule_1D))` approximation types
* Add docs on difference between StartUpDG.jl and Nodal DG Matlab codes

#### Changed
* removed `approximationType` as a property of `RefElemData` 
* Made `Base.show` output less verbose. 
* Deprecated CamlCase Gmsh read functions

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

