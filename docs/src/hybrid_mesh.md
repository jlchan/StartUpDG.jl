# `Hybrid meshes` 

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
@unpack x, y = md
u = @. x^3 - x^2 * y + 2 * y^3
dudx = @. 3 * x^2 - 2 * x * y

# compute local derivatives
@unpack rxJ, sxJ, J = md
dudr, duds = similar(md.x), similar(md.x)
dudr.Quad .= rds[Quad()].Dr * u.Quad
duds.Quad .= rds[Quad()].Ds * u.Quad
dudr.Tri .= rds[Tri()].Dr * u.Tri
duds.Tri .= rds[Tri()].Ds * u.Tri

@show norm(@. dudx - (rxJ * dudr + sxJ * duds) / J) # should be O(1e-14)
```

The main difference in the representation of hybrid meshes compared with standard `MeshData` objects
is the use of [ComponentArrays.jl](https://github.com/jonniedie/ComponentArrays.jl) as storage for the
geometric coordinates. These arrays have "fields" corresponding to the element type, for example
```julia
md.x.Tri
md.x.Quad
```
but can still be indexed as linear arrays. 

The `mapP` field behaves similarly. If we interpolate the values of `u` for each element type to surface
quadrature nodes, we can use `mapP` to linearly index into the array to find neighbors. 
```julia
uf = similar(md.xf)
uf.Quad .= rds[Quad()].Vf * u.Quad
uf.Tri .= rds[Tri()].Vf * u.Tri
uf[md.mapP] # this returns the exterior node values
```

