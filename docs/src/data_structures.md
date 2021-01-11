# Data structures `RefElemData` and `MeshData`

The structs `RefElemData{Dim,ElemShape,...}` and `MeshData{Dim}` contain fields which are useful for computing right hand sides for DG methods in a matrix-free fashion. `RefElemData` is immutable, while `MeshData` is mutable (in order to accomodate changes to geometric terms, e.g., curving meshes).

## `RefElemData`

`RefElemData` contains the following fields
* `elemShape::ElemShape`: element shape. `Line, Tri, Quad, Hex` currently supported.
* `Nfaces`: number of faces
* `fv`: list of vertices defining faces, e.g., `[1,2],[2,3],[3,1]` for a triangle
* `rst::NTuple{Dim,...}`: list of high order nodal interpolation points
* `rstq::NTuple{Dim,...}`,`wq`,`Vq`: volume quadrature points, weights, and interpolation matrix
* `rstp::NTuple{Dim,...}`, `Vp`: list of plotting points and plotting interpolation matrix
* `rstf::NTuple{Dim,...}`,`wf`,`Vf`: face quadrature points, weights, and interpolation matrix
* `nrstJ::NTuple{Dim,...}`: outward reference normals, scaled by face Jacobian
* `M`: mass matrix computed using quadrature
* `Pq`: quadrature-based ``L^2`` projection matrix
* `Drst::NTuple{Dim,...}`, `LIFT`: differentiation and lifting matrices

This list is incomplete; other fields are currently just used for internal computations.

By default, `RefEemData` is constructed for a nodal basis (in order to facilitate curved meshes, connectivity, etc).

Mass, differentiation, lifting, and interpolation matrices specialize on the type of matrix. For example, these matrices are dense `Matrix{T}` type for lines and triangles, but are stored as sparse matrices for quadrilaterals and hexahedra.

## `MeshData`

`MeshData` contains the following fields
* `K`: number of elements
* `FToF`: indexing array for face-to-face connectivity
* `xyz::NTuple{Dim,...}`: nodal interpolation points mapped to physical elements
* `xyzq::NTuple{Dim,...}, wJq`: volume quadrature points/weights mapped to physical elements
* `xyzf::NTuple{Dim,...}`: face quadrature points mapped to physical elements
* `mapP,mapB`: indexing arrays for inter-element node connectivity (`mapP`) and for extracting boundary nodes from the list of face nodes `xyzf` (`mapB`).
* `rstxyzJ::SMatrix{Dim,Dim}`: volume geometric terms ``G_{ij} = \frac{\partal x_i}{\partial \hat{x}_j}``
* `J,sJ`: volume and surface Jacobians evaluated at interpolation points and surface quadrature points, respectively.
* `nxyzJ::NTuple{Dim,...}`: scaled outward normals evaluated at surface quadrature points.

## Fieldname aliases

`Base.getproperty` is overloaded, so you can reference `rd.Dr` instead of `md.Drst[1]` or `md.rxJ` instead of `md.rstxyzJ[1,1]`. `@unpack` also works with aliased fieldnames.

The aliases for `RefElemData` are below:
1. `r` for `rst[1]` (similarly for `s` and `t`)
2. `rq` for `rstq[1]` (similarly for `sq` and `tq`)
3. `rf` for `rstf[1]` (similarly for `sf` and `tf`)
4. `rp` for `rstp[1]` (similarly for `sp` and `tp`)
5. `Dr` for `Drst[1]` (similarly for `Ds` and `Dt`)

The aliases for `MeshData` are below:
1. `VX` for `VXYZ[1]` (similarly for `VY`, `VZ`)
2. `x` for `xyz[1]` (similarly for `y`, `z`)
3. `xf` for `xyzf[1]` (similarly for `yf`, `zf`)
4. `xq` for `xyzq[1]` (similarly for `yq`, `zq`)
5. `nxJ` for `nxyzJ[1]` (similarly for `nyJ`, `nzJ`)
6. `rxJ` for `rstxyzJ[1,1]` (similarly for `sxJ,txJ,ryJ,syJ,tyJ,rzJ,szJ,tzJ`)
