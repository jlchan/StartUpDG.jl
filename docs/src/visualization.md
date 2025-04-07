# Basic visualization

StartUpDG.jl has utilities for exporting fields to VTK for visualization using [Paraview](https://www.paraview.org/). Recall that StartUpDG.jl represents solution fields as `num_points` by `num_elements` arrays, where `num_points = rd.Np` and `num_elements = md.num_elements`. We can export data to a `.vtu` file using `export_to_vtk` as follows:
```julia
rd = RefElemData(Tri(), 3)
md = MeshData(uniform_mesh(Tri(), 2), rd)
(; x, y) = md
u = @. (x + y)^2
v = @. x^2 + y^2
vtu_name = export_to_vtk(rd, md, [u, v], "uv.vtu")
```
If we wish to provide labels for each field, we can do so using a `Dict`:
```julia
rd = RefElemData(Tri(), 3)
md = MeshData(uniform_mesh(Tri(), 2), rd)
(; x, y) = md
u = @. (x + y)^2
v = @. x^2 + y^2
vtu_name = export_to_vtk(rd, md, Dict("u" => u, "v" => v), "uv.vtu")
```
By default, `export_to_vtk` assumes the solution is represented using values at interpolation nodes and interpolates the values to equispaced nodes on each element. Setting the keyword argument `equi_dist_nodes` bypasses this interpolation step, for example
```julia
vtu_name = export_to_vtk(rd, md, Dict("u" => u, "v" => v), "uv.vtu"; equi_dist_nodes=false)
```


