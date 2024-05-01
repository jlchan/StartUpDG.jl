# Additional mesh types

In addition to more "standard" mesh types, StartUpDG.jl also has experimental support for hybrid and cut-cell meshes. Both are currently restricted to two dimensional domains.

## Hybrid meshes

!!! warning "Experimental implementation"
    This is an experimental feature and may change in future releases.

There is initial support for hybrid meshes in StartUpDG.jl. The following is a 
short example where we interpolate a polynomial and compute its derivative.
```julia
rds = RefElemData((Tri(), Quad()), N = 3)

# Simple hybrid mesh for testing
#   1  7______8______9
#      |      | 3  / |
#      |   4  |  / 5 |
#   0  4 ---- 5 ---- 6 
#      |      |      |
#      |   1  |   2  |
#   -1 1 ---- 2 ---- 3
#     -1      0      1
VX = [-1; 0; 1; -1; 0; 1; -1; 0; 1]
VY = [-1; -1; -1; 0; 0; 0; 1; 1; 1]
EToV = [[1 2 4 5], [2 3 5 6], [5 8 9], [4 5 7 8], [9 6 5]]

md = MeshData(VX, VY, EToV, rds)

# test that the local derivatives of a polynomial recover the exact derivative
(; x, y ) = md
u = @. x^3 - x^2 * y + 2 * y^3
dudx = @. 3 * x^2 - 2 * x * y

# compute local derivatives
(; rxJ, sxJ, J ) = md
dudr, duds = similar(md.x), similar(md.x)
dudr.Quad .= rds[Quad()].Dr * u.Quad
duds.Quad .= rds[Quad()].Ds * u.Quad
dudr.Tri .= rds[Tri()].Dr * u.Tri
duds.Tri .= rds[Tri()].Ds * u.Tri

@show norm(@. dudx - (rxJ * dudr + sxJ * duds) / J) # should be O(1e-14)
```

The main difference in the representation of hybrid meshes compared with standard `MeshData` objects 
is the use of `NamedArrayPartition` arrays as storage for the geometric coordinates. These arrays have 
"fields" corresponding to the element type, for example
```julia
md.x.Tri
md.x.Quad
```
but can still be indexed as linear arrays. 

The `mapP` field behaves similarly. If we interpolate the values of `u` for each element type to surface
quadrature nodes, we can use `mapP` to linearly index into the array to find neighbors. 
```julia
uf = similar(md.xf)
uf.Quad .= rds.Quad.Vf * u.Quad
uf.Tri .= rds.Tri.Vf * u.Tri
uf[md.mapP] # this returns the exterior node values
```

## Cut Meshes

!!! warning "Experimental implementation"
    This is an experimental feature and may change in future releases.

Initial support for cut-cell meshes is available via [PathIntersections.jl](https://github.com/cgt3/pathintersections.jl). By passing in a tuple of curves (defined as parametrized functions of one coordinate, see PathIntersections.jl documentation for more detail), StartUpDG.jl can compute a `MeshData` for a cut-cell mesh. 
```julia
circle = PresetGeometries.Circle(R=0.33, x0=0, y0=0)
cells_per_dimension_x, cells_per_dimension_y = 4, 4

rd = RefElemData(Quad(), N=3)
md = MeshData(rd, (circle, ), cells_per_dimension_x, cells_per_dimension_y, Subtriangulation(); precompute_operators=true)
```
Here, the final argument `quadrature_type = Subtriangulation()` determines how the quadrature on cut cells is determined. For `Subtriangulation()`, the quadrature on cut cells is constructed from a curved isoparametric subtriangulation of the cut cell. The number of quadrature points on a cut cell is then reduced (while preserving positivity) using Caratheodory pruning. If not specified, the `quadrature_type` argument defaults to `Subtriangulation()`. 

Quadrature rules can also be constructed by specifying `quadrature_type = MomentFitting()`. The quadrature points on cut cells `md.xq.cut` are determined from sampling and a pivoted QR decomposition. This is not recommended, as it can be both slower, and the cut-cell quadrature weights `md.wJq.cut` are not guaranteed to be positive. 

The interpolation points on cut cells `md.x.cut` are determined from sampled points and a pivoted QR decomposition. 

The optional keyword argument `precompute_operators` specifies whether to precompute differentiation, face interpolation, mass, and lifting matrices for each cut cell. If 
`precompute_operators=true`, these are stored in `md.mesh_type.cut_cell_operators`. 

As with hybrid meshes, the nodal coordinates `md.x`, `md.y` are `NamedArrayPartition`s with `cartesian` and `cut` fields. For example, `md.x.cartesian` and `md.x.cut` are the x-coordinates of the Cartesian and cut cells, respectively. Likewise, `md.mapP` indexes linearly into the array of face coordinates and specifies exterior node indices. For example, we can interpolate a function to face nodes and compute exterior values via the following code:
```julia
(; x, y) = md
u = @. x^rd.N - x * y^(rd.N-1) - x^(rd.N-1) * y + y^rd.N # some random function 

# interpolate the solution to face nodes
uf = similar(md.xf) 
uf.cartesian .= rd.Vf * u.cartesian
for e in 1:size(md.x.cut, 2)
    ids = md.mesh_type.cut_face_nodes[e]
    Vf = md.mesh_type.cut_cell_operators.face_interpolation_matrices[e]
    uf.cut[ids] .= Vf * u.cut[:, e]
end

uf[md.mapP] # these are "exterior" values for each entry of `uf`
```