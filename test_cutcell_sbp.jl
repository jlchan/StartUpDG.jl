using NodesAndModes: face_basis
using Plots
using LinearAlgebra
using StartUpDG

N = 3
quad_rule_face = gauss_quad(0, 0, 2 * N + 1)

rd = RefElemData(Quad(), N; quad_rule_face)

# objects = (PresetGeometries.Rectangle(Lx=.4, Ly=2.5, x0=1.2), )
# vx = [-1, 0, 1.1]
# vy = [-1, 0, 1]
# md = MeshData(rd, objects, vx, vy; precompute_operators=true)

cells_per_dimension = 2
circle = PresetGeometries.Circle(R=0.66, x0=0, y0=0)
md = MeshData(rd, (circle, ), cells_per_dimension; quad_rule_face)#, precompute_operators=true)

mt = md.mesh_type
(; cutcells) = mt.cut_cell_data

using Triangulate, StaticArrays


# degree N(N-1) + 2N-2 polynomial = N(N-1) + 2(N-1) = (N+2) * (N-1)
N_phys_frame_geo = max(2 * rd.N, (rd.N-1) * (rd.N + 2)) 

rd_line = RefElemData(Line(), N=rd.N, quad_rule_vol = quad_rule_face)

rd_tri = RefElemData(Tri(), N=rd.N, quad_rule_vol=NodesAndModes.quad_nodes_tri(N_phys_frame_geo))
fv = NodesAndModes.face_vertices(Tri())
tri_face_coords = (rd_tri.r[vec(rd_tri.Fmask)], rd_tri.s[vec(rd_tri.Fmask)])

# preallocate face interp node arrays
r1D = rd_line.r
tri_face_coords_x = zeros(length(r1D), 3)
tri_face_coords_y = zeros(length(r1D), 3)

# this operator performs a least squares fit, and is equivalent to isoparametric warp and blend
warp_face_points_to_interp = 
    face_basis(Tri(), rd_tri.N, rd_tri.rst...) / face_basis(Tri(), rd_tri.N, tri_face_coords...) 

using StartUpDG: map_to_interval

function compute_geometric_determinant_J(x, y, Dr, Ds)
    xr, xs = Dr * x, Ds * x
    yr, ys = Dr * y, Ds * y
    J = @. -xs * yr + xr * ys
    return J
end

x_cutcells, y_cutcells, xq_cutcells, yq_cutcells, Jq_cutcells, wJq_cutcells = 
    ntuple(_ -> Matrix{eltype(rd_tri.wq)}[], 6)


for cutcell in cutcells

    # vxy = matrix of the form [x_coordinates, y_coordinates] 
    vxy = hcat(getindex.(cutcell.(cutcell.stop_pts[1:end-1]), 1), 
               getindex.(cutcell.(cutcell.stop_pts[1:end-1]), 2))
    (VX, VY), EToV = StartUpDG.triangulate_points(permutedims(vxy))

    xq, yq, Jq, wJq = ntuple(_ -> zeros(length(rd_tri.wq), size(EToV, 1)), 6)
    
    # loop over each triangle, map 1D interpolation points to faces
    for e in axes(EToV, 1)
        ids = view(EToV, e, :)
        for (face_index, face_vertices) in enumerate(SVector{2}.(fv))
            vertices_on_face = sort(ids[face_vertices])

            # map face interpolation points to a physical element. 

            # This assumes PathIntersections.jl uses a clockwise ordering of stop curve points.
            # Since StartUpDG uses a CCW ordering, we reverse the order for 
            for i in eachindex(r1D)
                # If the vertex indices are far apart, it's the last face/boundary curve
                if (x->abs(x[2]-x[1]))(vertices_on_face) == length(VX) - 1 
                    s = map_to_interval(r1D[i], reverse(cutcell.stop_pts[end-1:end])...)
                    point = cutcell(s)
                
                # if vertex indices are consecutive, it's a boundary face    
                elseif (x->x[2]-x[1])(vertices_on_face) == 1 
                    
                    curve_id = minimum(ids[face_vertices])
                    s = map_to_interval(r1D[i], reverse(cutcell.stop_pts[curve_id:curve_id+1])...)
                    point = cutcell(s)

                else # it's a non-boundary face, it's a straight line
                    point = SVector{2}.(map_to_interval(r1D[i], VX[ids[face_vertices]]...),
                                        map_to_interval(r1D[i], VY[ids[face_vertices]]...))
                end
                tri_face_coords_x[i, face_index] = point[1]
                tri_face_coords_y[i, face_index] = point[2]
            end
        end

        # this performs a least squares fit interpolation by the face basis. It's 
        # equivalent to isoparametric warp and blend if the face node locations are 
        # representable by the face basis (e.g., polynomial). 
        tri_warped_coords_x = warp_face_points_to_interp * vec(tri_face_coords_x) 
        tri_warped_coords_y = warp_face_points_to_interp * vec(tri_face_coords_y)
        
        rxJq_e, sxJq_e, ryJq_e, syJq_e, Jq_e = 
            StartUpDG.geometric_factors(tri_warped_coords_x, tri_warped_coords_y, 
                                        rd_tri.Vq * rd_tri.Dr, rd_tri.Vq * rd_tri.Ds)

        view(xq, :, e) .= rd_tri.Vq * tri_warped_coords_x
        view(yq, :, e) .= rd_tri.Vq * tri_warped_coords_y
        view(Jq, :, e) .= Jq_e
        @. wJq[:,e] = rd_tri.wq * Jq_e
    end
    push!(xq_cutcells, xq)
    push!(yq_cutcells, yq)
    push!(Jq_cutcells, Jq)
    push!(wJq_cutcells, wJq)
end

xq_pruned, yq_pruned, wJq_pruned = ntuple(_ -> Vector{Float64}[], 3)
for e in eachindex(xq_cutcells)
    V2N = vandermonde(mt.physical_frame_elements[e], 2 * rd.N, vec(xq_cutcells[e]), vec(yq_cutcells[e]))
    w = copy(vec(wJq_cutcells[e]))
    w_pruned, inds = StartUpDG.caratheodory_pruning_qr(V2N, w)

    V = vandermonde(mt.physical_frame_elements[e], rd.N, vec(xq_cutcells[e]), vec(yq_cutcells[e]))
    # @show size(V[inds,:])
    # @show length(w), length(inds)
    @show norm(V' * diagm(w) * V - V' * diagm(w_pruned) * V)

    push!(yq_pruned, vec(yq_cutcells[e])[inds])
    push!(xq_pruned, vec(xq_cutcells[e])[inds])
    push!(wJq_pruned, vec(wJq_cutcells[e])[inds])
end

# recompute normals on cut cells
cutcell = cutcells[1]
plot()
xf, yf, nxJ, nyJ = ntuple(_ -> zeros(size(rd_line.Vq, 1), length(cutcell.stop_pts)-1), 4)
for f in 1:length(cutcell.stop_pts)-1
    points = map(s -> cutcell(map_to_interval(s, cutcell.stop_pts[f:f+1]...)), r1D)
    x = getindex.(points, 1)
    y = getindex.(points, 2)

    # compute tangent vector using polynomial mapping
    dxdr = rd_line.Vq * rd_line.Dr * x
    dydr = rd_line.Vq * rd_line.Dr * y

    tangent_vector = SVector.(dxdr, dydr)
    scaling = (cutcell.stop_pts[f+1] - cutcell.stop_pts[f]) / 2
    scaled_normal = SVector.(-dydr, dxdr)

    @. nxJ[:, f] = getindex(scaled_normal, 1)
    @. nyJ[:, f] = getindex(scaled_normal, 2)

    # interp face coordinates to face quad nodes
    xf[:, f] .= rd_line.Vq * x
    yf[:, f] .= rd_line.Vq * y

    scatter!(rd_line.Vq * x, rd_line.Vq * y)
    quiver!(rd_line.Vq * x, rd_line.Vq * y, 
            quiver=(getindex.(scaled_normal, 1), getindex.(scaled_normal, 2)))
end
plot!(leg=false, ratio=1)

plot()
for e in eachindex(xq_cutcells)
    scatter!(vec(xq_cutcells[e]), vec(yq_cutcells[e]), label="Reference quadrature"); 
    scatter!(xq_pruned[e], yq_pruned[e], markersize=8, marker=:circle, 
             z_order=:back, label="Caratheodory pruning", leg=false)
end
display(plot!(leg=false))



# # compute new normals for curved boundaries
# mapB_boundary_cartesian = findall(@. abs(abs(md.xf) - 1) < 1e2 * eps() || abs(abs(md.yf) - 1) < 1e4 * eps())
# mapB_cut = setdiff(md.mapB, mapB_boundary_cartesian)
# num_face_nodes = length(quad_rule_face[1])
# xf, yf = map(x -> reshape(x[mapB_cut], num_face_nodes, :), md.xyzf)
# D1D = rd_line.Vq * rd_line.Dr * rd_line.Pq 
# dxdr, dydr = D1D * xf, D1D * yf
# nxJ_boundary_cut, nyJ_boundary_cut = dydr, -dxdr

# test weak SBP property
(; x, y) = md
elem = md.mesh_type.physical_frame_elements[1]
VDM = vandermonde(elem, rd.N, x.cut[:, 1], y.cut[:, 1])
# xq, yq, wq = xq_pruned[1], yq_pruned[1], wJq_pruned[1]
xq, yq, wq = vec.((xq_cutcells[1], yq_cutcells[1], wJq_cutcells[1]))
Vq, Vxq, Vyq = map(A -> A / VDM, basis(elem, rd.N, xq, yq))
M = Vq' * diagm(wq) * Vq
Qx, Qy = Vq' * diagm(wq) * Vxq, Vq' * diagm(wq) * Vyq
Vf = vandermonde(elem, rd.N, xf, yf) / VDM

# Bx = Diagonal(vec(Diagonal(rd_line.wq) * reshape(md.nxJ.cut[md.mesh_type.cut_face_nodes[1]], length(rd_line.wq), :)))
# By = Diagonal(vec(Diagonal(rd_line.wq) * reshape(md.nyJ.cut[md.mesh_type.cut_face_nodes[1]], length(rd_line.wq), :)))
Bx = Diagonal(vec(Diagonal(rd_line.wq) * reshape(nxJ, length(rd_line.wq), :)))
By = Diagonal(vec(Diagonal(rd_line.wq) * reshape(nyJ, length(rd_line.wq), :)))

e = ones(size(Qx, 2))
@show norm(sum(Qx - Vf' * Bx * Vf, dims=1))
@show norm(sum(Qy - Vf' * By * Vf, dims=1))
# display(sum(Qx - Vf' * Bx * Vf, dims=1))
# display(e' * ((Qx + Qx') - Vf' * Bx * Vf))
# display(e' * ((Qy + Qy') - Vf' * By * Vf))
