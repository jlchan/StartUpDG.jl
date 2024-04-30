using NodesAndModes: face_basis
using Plots
using LinearAlgebra
using StaticArrays
using StartUpDG
using PathIntersections

N = 3
quad_rule_face = gauss_lobatto_quad(0, 0, N)
rd = RefElemData(Quad(), N; quad_rule_face)

cells_per_dimension = 4
circle = PresetGeometries.Circle(R=0.66, x0=0, y0=0)
objects = (circle, )

# TODO: remove eventually
cells_per_dimension_x = cells_per_dimension
cells_per_dimension_y = cells_per_dimension
vx = LinRange(-1, 1, cells_per_dimension_x + 1)
vy = LinRange(-1, 1, cells_per_dimension_y + 1)    

N = rd.N

region_flags, cutcells = StartUpDG.calculate_cutcells(vx, vy, objects)

####################################################
#          Construct Cartesian cells               # 
####################################################

xf_cartesian, yf_cartesian, nxJ_cartesian, nyJ_cartesian, wf_cartesian = 
    StartUpDG.construct_cartesian_surface_quadrature(vx, vy, region_flags, quad_rule_face)    

xq_cartesian, yq_cartesian, wJq_cartesian = 
    StartUpDG.construct_cartesian_volume_quadrature(vx, vy, region_flags, 
                                                    (rd.rq, rd.sq, rd.wq))

####################################################
#          Construct cut cells stuff               # 
####################################################

region_flags, cutcells = StartUpDG.calculate_cutcells(vx, vy, objects)

physical_frame_elements = 
    StartUpDG.construct_physical_frame_elements(region_flags, vx, vy, cutcells)

xf_cut, yf_cut, nxJ_cut, nyJ_cut, wf_cut, cut_face_node_indices = 
    StartUpDG.construct_cut_surface_quadrature(N, cutcells, quad_rule_face)

xq_pruned, yq_pruned, wJq_pruned = 
    StartUpDG.construct_cut_volume_quadrature(N, cutcells, physical_frame_elements)

#######################################
#        test weak SBP property       #
#######################################

e = 9
elem = physical_frame_elements[e]
xq, yq, wq = xq_pruned[:,e], yq_pruned[:,e], wJq_pruned[:,e]

Vq, Vxq, Vyq = basis(elem, rd.N, xq, yq)
Vf = vandermonde(elem, rd.N, xf_cut[e], yf_cut[e]) 
Qx, Qy = Vq' * diagm(wq) * Vxq, Vq' * diagm(wq) * Vyq

Bx = Diagonal(wf_cut[e] .* nxJ_cut[e])
By = Diagonal(wf_cut[e] .* nyJ_cut[e])

# represent constant vector in basis
ee = Vq \ ones(size(Vq, 1))
@show norm(ee' * Qx - ee' * Vf' * Bx * Vf)
@show norm(ee' * Qy - ee' * Vf' * By * Vf)

# plot pruned cells

# volume quadrature should be exact for degree N(N-1) + 2N-2 polynomials, 
# or N(N-1) + 2(N-1) = (N+2) * (N-1) degrees
N_phys_frame_geo = max(2 * N, (N-1) * (N + 2))     
target_degree = 2 * N
rd_tri = RefElemData(Tri(), Polynomial(MultidimensionalQuadrature()), N, 
                    quad_rule_vol=NodesAndModes.quad_nodes_tri(N_phys_frame_geo))
Np_target = StartUpDG.Np_cut(target_degree)

plot()
for e in eachindex(cutcells)
    xq, yq, _ = StartUpDG.subtriangulated_cutcell_quadrature(cutcells[e], rd_tri)
    scatter!(vec(xq), vec(yq), label="Reference quadrature"); 
    scatter!(xq_pruned[:,e], yq_pruned[:,e], markersize=8, marker=:circle, 
             z_order=:back, label="Caratheodory pruning", leg=false)
end
display(plot!(leg=false))
