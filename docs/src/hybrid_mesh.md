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
dudr = ArrayPartition((getproperty.(values(rds), :Dr) .* u.x)...)
duds = ArrayPartition((getproperty.(values(rds), :Ds) .* u.x)...)
@show norm(@. dudx - (rxJ * dudr + sxJ * duds) / J) # should be O(1e-14)
```

Todo: add example usage of `mapP`. 
