# `RefElemData` type

[`RefElemData`](@ref) contains the following fields
* `elemShape::ElemShape`: element shape. `Line, Tri, Quad, Hex` currently supported.
* `approxType`: approximation type. Defaults to `Polynomial()`, but `SBP()` is also supported (see [RefElemData based on SBP finite differences](@ref)).
* `Nfaces`: number of faces on a given type of reference element.
* `fv`: list of vertices defining faces, e.g., `[1,2],[2,3],[3,1]` for a triangle
* `rst::NTuple{Dim,...}`: tuple of vectors of length `N_p`, each of which contains coordinates of degree ``N`` optimized polynomial interpolation points.
* `rstq::NTuple{Dim,...}`,`wq`,`Vq`: tuple of volume quadrature points, vector of weights, and quadrature interpolation matrix. Each element of `rstq` and `wq` are vectors of length ``N_q``, and `Vq` is a matrix of size ``N_q \times N_p``.
* `N_{\rm plot}`: the degree which determines the number of plotting points `N_{p,{\rm plot}}`.
* `rstp::NTuple{Dim,...}`, `Vp`: tuple of plotting points and plotting interpolation matrix. Each element of `rstp` is a vector of length ``N_{p,{\rm plot}}``, and `Vp` is a matrix of size ``N_{p,{\rm plot}} \times N_p``.
* `rstf::NTuple{Dim,...}`,`wf`,`Vf`: tuple of face quadrature points, weights, and face interpolation matrix. Each element of `rstf` and `wf` are vectors of length ``N_f``, and `Vf` is a matrix of size ``N_f \times N_p``.
* `nrstJ::NTuple{Dim,...}`: tuple of outward reference normals, scaled by face Jacobian. Each element is a vector of length ``N_f``.
* `M`: mass matrix computed using quadrature. Size ``N_p \times N_p``
* `Pq`: quadrature-based ``L^2`` projection matrix. Size ``N_p \times N_q``.
* `Drst::NTuple{Dim,...}`, `LIFT`: differentiation and lifting matrices. Differentiation matrices are size ``N_p \times N_p``, while lift matrices are size ``N_p\times N_f``.

This list is incomplete; other fields are stored or accessible but currently only used for internal computations.

Mass, differentiation, lifting, and interpolation matrices specialize on the type of matrix. For example, these matrices are dense `Matrix{T}` type for lines and triangles, but might be stored as sparse matrices for quadrilaterals and hexahedra.

# Setting up `rd::RefElemData`

The struct `rd::RefElemData` contains data for a given element type. Currently, five types of reference elements are supported: `Line`, `Tri`, `Quad`, `Tet`, and `Hex`.

To initalize a `RefElemData`, just specify the element type and polynomial degree.
```julia
N = 3
rd = RefElemData(Line(),N)
rd = RefElemData(Tri(),N)
rd = RefElemData(Quad(),N)
rd = RefElemData(Tet(),N)
rd = RefElemData(Hex(),N)
```

## Specifying different quadrature rules.

By default, [`RefElemData`](@ref) initializes volume and surface quadrature rules to be the minimum rules which exactly integrate the unweighted volume and surface mass matrices. If different quadrature rules are desired, they can be specified as follows:
```julia
N = 3

# create degree N tensor product Gauss-Lobatto rule
r1D,w1D = gauss_lobatto_quad(0,0,N)
rq,sq = vec.(StartUpDG.meshgrid(r1D))
wr,ws = vec.(StartUpDG.meshgrid(w1D))
wq = @. wr*ws

rd = RefElemData(Quad(),N; quad_rule_vol =(rq,sq,wq),  
                           quad_rule_face=(r1D,w1D))
```
This results in a DG spectral element method (DG-SEM) discretization, with a diagonal lumped mass matrix and differentiation matrices which satisfy a summation-by-parts property.

By default, `RefElemData` is constructed for a nodal basis (in order to facilitate curved meshes, connectivity, etc). There is not functionality to change interpolation nodes, since these transformations can be performed as algebraic changes of basis after setting up a `RefElemData`. 

## RefElemData based on SBP finite differences

It is also possible to construct a [`RefElemData`](@ref) based on [multi-dimensional SBP finite difference operators](https://doi.org/10.1137/15M1038360). These utilize nodes constructed by [Tianheng Chen and Chi-Wang Shu](https://doi.org/10.1016/j.jcp.2017.05.025), [Ethan Kubatko](https://sites.google.com/site/chilatosu/ethan-bio), and [Jason Hicken](https://doi.org/10.1007/s10915-020-01154-8).

Some examples:
```julia
rd = RefElemData(Tri(), SBP(), 2)
rd = RefElemData(Quad(), SBP(), 2)
rd = RefElemData(Tri(), SBP(), 2; quadrature_strength=4, quad_rule_face=:Legendre) 
```
Quadrature rules of both degree `2*N-1` (up to `N=6`) and `2*N` (up to `N=4`) are supported on triangles. For `Line`, `Quad`, and `Hex` elements, `RefElemData(...,SBP(),N)` is the same as the `RefElemData` for a DG-SEM discretization, though some fields are specialized for the SBP type. 

These SBP-based `RefElemData` objects can also be used to initialize a mesh (for example, `md = MeshData(uniform_mesh(rd.elementType,4)...,rd)`). 