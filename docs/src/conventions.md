# Background

Most high order finite element methods rely on a decomposition of a domain into a mesh of "elements" (e.g., triangles or quadrilaterals in 2D, hexahedra or tetrahedra in 3D). Each "physical" element in a mesh is assumed to be the image of single "reference" element under some geometric mapping. Using the chain rule and changes of variables, one can evaluate integrals and derivatives using only operations on the reference element and some geometric mapping data. This transformation of operations on _all_ elements to a single reference element make finite element methods efficient. 

We use the convention that coordinates on the reference element are ``r`` in 1D, ``r, s`` in 2D, or ``r, s, t`` in 3D. Physical coordinates use the standard conventions ``x``, ``x, y``, and ``x, y, z`` in 1D, 2D, and 3D. 

![Mapping](assets/mapping_diagram.png)

Derivatives of reference coordinates with respect to physical coordinates are abbreviated, e.g., ``\frac{\partial r}{\partial x} = r_x``. Additionally, ``J`` is used to denote the determinant of the Jacobian of the reference-to-physical mapping. 

## Assumptions

We make a few simplifying assumptions about the mesh:
* meshes are _conforming_ (e.g., each face of an element is shared with at most one other element). 
* the geometric mapping from reference to physical elements is the same degree polynomial as the approximation space on the reference element (e.g., the mapping is isoparametric). 

Initial experimental support for hybrid, cut-cell, and non-conforming meshes in two dimensions is also available. Please see the corresponding test sets `test/hybrid_mesh_tests.jl`, `test/cut_mesh_tests.jl`, and `noncon_mesh_tests.jl` for examples. 

## Code conventions

`StartUpDG.jl` exports structs `RefElemData{Dim, ElemShape, ...}` (which contains data associated with the reference element, such as interpolation points, quadrature rules, face nodes, normals, and differentiation/interpolation/projection matrices) and `MeshData{Dim}` (which contains geometric data associated with a mesh). These are currently used for evaluating DG formulations in a matrix-free fashion. These structs contain fields similar to those in `Globals1D, Globals2D, Globals3D` in the NDG book codes. 

We use the following code conventions:
* variables `r, s, t` and `x, y, z` correspond to values at nodal interpolation points. 
* variables ending in `q` (e.g., `rq, sq,...` and `xq, yq,...`) correspond to values at volume quadrature points. 
* variables ending in `f` (e.g., `rf, sf,...` and `xf, yf,...`) correspond to values at face quadrature points. 
* variables ending in `p` (e.g., `rp, sp,...`) correspond to equispaced plotting nodes.
* `Dr, Ds, Dt` matrices are nodal differentiation matrices with respect to the ``r, s, t`` coordinates. For example, `Dr * f.(r, s)` approximates the derivative of ``f(r, s)`` at nodal points. 
* `V` matrices correspond to interpolation matrices from nodal interpolation points. For example, `Vq` interpolates to volume quadrature points, `Vf` interpolates to face quadrature points, `Vp` interpolates to plotting nodes. 
* geometric quantities in `MeshData` are stored as matrices of dimension ``\text{number of points per element } \times \text{number of elements}``.
