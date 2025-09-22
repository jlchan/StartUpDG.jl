"""
`CutCellMesh` is used in the `MeshData` field `mesh_type` for cut cell meshes.

The field `physical_frame_elements` is a container with shifting/scaling information for 
each element. We evaluate the physical basis over each element by applying a shifting and 
scaling of the physical coordinates. The resulting shifted/scaled coordinates then fall 
into the reference element and can be used to evaluate a reference element basis. 

The field `cut_face_nodes` is a container whose elements are indices of face nodes for a 
cut element. In other words, `md.xf.cut[cut_face_nodes[1]]` returns the face nodes of the 
first element. 

We assume all cut elements have the same number of volume quadrature points (which is at 
least the dimension of a degree 2N polynomial space). 

The field `objects` contains a tuple of the objects used to define the cut region.

The field `cut_cell_operators` contains optionally precomputed operators (mass, differntiation, 
face interpolation, and lifting operators). 

The field `cut_cell_data` contains additional data from PathIntersections.
"""
struct CutCellMesh{T1, T2, T3, T4, T5, T6}
    physical_frame_elements::T1
    cut_face_nodes::T2
    objects::T3
    cut_cell_operators::T4
    cut_cell_data::T5
    quadrature_type::T6
end

function Base.show(io::IO, ::MIME"text/plain", md::MeshData{DIM, <:CutCellMesh}) where {DIM}
    @nospecialize md
    print(io,"Cut-cell MeshData of dimension $DIM with $(md.num_elements) elements " * 
             "($(num_cartesian_elements(md)) Cartesian, $(num_cut_elements(md)) cut)")
end

# maps x ∈ [-1,1] to [a,b]
map_to_interval(x, a, b) = @. a + (b-a) * 0.5 * (1 + x)

function neighbor_across_face(f, ex, ey)
    if f==1
        return (ex-1, ey)
    elseif f==2
        return (ex+1, ey)
    elseif f==3
        return (ex, ey-1)
    elseif f==4
        return (ex, ey+1)
    else
        error("Face index f = $f > 4; too large.")
    end
end

function compute_element_indices(region_flags)

    cells_per_dimension_x, cells_per_dimension_y = size(region_flags)
   
    element_indices = fill!(similar(region_flags), zero(eltype(region_flags)))

    # count Cartesian elements
    e = 1
    for ex in 1:cells_per_dimension_x, ey in 1:cells_per_dimension_y    
        if is_Cartesian(region_flags[ex, ey])
            element_indices[ex, ey] = e
            e += 1
        end
    end

    # reset counter for cut cells
    e = 1
    for ex in 1:cells_per_dimension_x, ey in 1:cells_per_dimension_y    
        if is_cut(region_flags[ex, ey])
            element_indices[ex, ey] = e
            e += 1
        end
    end

    return element_indices
end

function is_inside_domain(ex, ey, regions) 
    num_cells_x, num_cells_y = size(regions)
    inside_cartesian_domain = (0 < ex <= num_cells_x) & (0 < ey <= num_cells_y) 
    if inside_cartesian_domain 
        # check if a cell is also inside the cut geometry
        return inside_cartesian_domain & (regions[ex, ey] >= 0) 
    else
        return false
    end
end

function generate_sampling_points(objects, elem, Np_target::Int, N_sampled::Int)

    r_sampled, s_sampled = equi_nodes(Quad(), N_sampled) # oversampled nodes

    # Here, we map sampled points on [-1, 1]^2 to the smallest bounding box around the cut element. 
    # These mapped points should still lie within the background element but be clustered around the 
    # cut cell, and should be more efficient than a uniform sampling of the background element.
    x_sampled, y_sampled = map_nodes_to_cutcell_boundingbox(elem, r_sampled, s_sampled)
    is_in_domain = fill(true, length(x_sampled))
    for (index, point) in enumerate(zip(x_sampled, y_sampled))
        is_in_domain[index] = !any(map(obj -> PathIntersections.is_contained(obj, point), objects))
    end
    
    # increase number of background points until we are left with `Np_target` sampling points 
    while sum(is_in_domain) < Np_target

        N_sampled *= 2 # double degree of sampling
        r_sampled, s_sampled = equi_nodes(Quad(), N_sampled) # oversampled nodes
        x_sampled, y_sampled = map_nodes_to_cutcell_boundingbox(elem, r_sampled, s_sampled)

        # check if all the points are in all the objects
        is_in_domain = fill(true, length(x_sampled))
        for (index, point) in enumerate(zip(x_sampled, y_sampled))
            is_in_domain[index] = !any(map(obj -> PathIntersections.is_contained(obj, point), objects))
        end
    end

    ids_in_element = findall(is_in_domain)

    return x_sampled[ids_in_element], y_sampled[ids_in_element]
end

is_Cartesian(flag) = flag==0 ? true : false
is_cut(flag) = flag > 0

# returns the 1D quadrature used to build a RefElemData surface quadrature 
function get_1d_quadrature(rd::RefElemData{2, Quad})
    rf = reshape(rd.rf, :, num_faces(rd.element_type))
    wf = reshape(rd.wf, :, num_faces(rd.element_type))
    
    # face ordering on a quad is -/+ x, -/+ y. face 3 = -y
    return rf[:, 3], wf[:, 3]
end

function calculate_cutcells(vx, vy, objects, ds = 1e-3, arc_tol = 1e-10, corner_tol = 1e-10)

    stop_pts = PathIntersections.find_mesh_intersections((vx, vy), objects, ds, arc_tol, corner_tol,
                                                          closed_list=true, closure_tol=1e-12)

    # Calculate cutcells
    region_flags, cutcell_indices, cutcells = 
        PathIntersections.define_regions((vx, vy), objects, stop_pts, binary_regions=true)

    cells_per_dimension_x = length(vx) - 1
    cells_per_dimension_y = length(vy) - 1

    # sort the vector of cut cells so that they match the ordering when 
    # iterating through Cartesian mesh indices via (ex, ey).
    cutcell_ordering = zeros(Int, length(cutcells))
    e = 1
    for ex in 1:cells_per_dimension_x, ey in 1:cells_per_dimension_y 
        if is_cut(region_flags[ex, ey])
            cutcell_ordering[e] = cutcell_indices[ex, ey] 
            e += 1
        end        
    end
    permute!(cutcells, cutcell_ordering)

    return region_flags, cutcells
end

"""
    subtriangulated_cutcell_quadrature(cutcell, rd_tri::RefElemData, 
                                       r1D = gauss_lobatto_quad(0,0,rd_tri.N))

Constructs a quadrature from subtriangulations. The degree of both the quadrature 
rule and isoparametric mapping are specified by `rd_tri`. 

The optional argument `r1D` specifies the nodal points used for constructing a 
curved mapping via interpolatory warp and blend. 
"""
function subtriangulated_cutcell_quadrature(cutcell, rd_tri, 
                                            r1D = NodesAndModes.nodes(Line(), rd_tri.N))
    # vxy = matrix of the form [x_coordinates, y_coordinates] 
    vxy = hcat(getindex.(cutcell.(cutcell.stop_pts[1:end-1]), 1), 
               getindex.(cutcell.(cutcell.stop_pts[1:end-1]), 2))
    (VX, VY), EToV = StartUpDG.triangulate_points(permutedims(vxy))

    # allocate arrays for face interp coordinates
    tri_face_coords = (rd_tri.r[vec(rd_tri.Fmask)], rd_tri.s[vec(rd_tri.Fmask)])
    tri_face_coords_x = zeros(length(r1D), 3)
    tri_face_coords_y = zeros(length(r1D), 3)

    # allocate output arrays
    xq, yq, wJq = ntuple(_ -> zeros(length(rd_tri.wq), size(EToV, 1)), 3)

    # loop over each triangle, map 1D interpolation points to faces
    for e in axes(EToV, 1)
        ids = view(EToV, e, :)
        for (face_index, face_vertices) in enumerate(SVector{2}.(rd_tri.fv))
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

        # this operator performs a least squares fit, and is equivalent to isoparametric 
        # warp and blend. can be moved outside of the cutcell loop for efficiency. 
        warp_face_points_to_interp = 
            face_basis(Tri(), rd_tri.N, rd_tri.rst...) / 
            face_basis(Tri(), rd_tri.N, tri_face_coords...) 

        # this performs a least squares fit interpolation by the face basis. It's 
        # equivalent to isoparametric warp and blend if the face node locations are 
        # representable by the face basis (e.g., polynomial). 
        tri_warped_coords_x = warp_face_points_to_interp * vec(tri_face_coords_x) 
        tri_warped_coords_y = warp_face_points_to_interp * vec(tri_face_coords_y)
        
        _, _, _, _, Jq_e = 
            StartUpDG.geometric_factors(tri_warped_coords_x, tri_warped_coords_y, 
                                        rd_tri.Vq * rd_tri.Dr, rd_tri.Vq * rd_tri.Ds)

        view(xq, :, e) .= rd_tri.Vq * tri_warped_coords_x
        view(yq, :, e) .= rd_tri.Vq * tri_warped_coords_y
        @. wJq[:,e] = rd_tri.wq * Jq_e
    end
    
    return xq, yq, wJq
end

# quadrature type for cut cell meshes
struct Subtriangulation end

"""
    MeshData(rd::RefElemData, objects, 
             vx::AbstractVector, vy::AbstractVector,
             quadrature_type=Subtriangulation(); 
             quad_rule_face=get_1d_quadrature(rd), 
             precompute_operators=false)

Constructor for MeshData utilizing moment fitting. Does not guarantee positive
quadrature weights, and is slower due to the use of adaptive sampling
to construct 

!!! Warning: this may be deprecated or removed in future versions. 
"""
MeshData(rd::RefElemData, objects, cells_per_dimension, 
         quadrature_type=Subtriangulation(); kwargs...) = 
    MeshData(rd, objects, 
             cells_per_dimension, cells_per_dimension, 
             quadrature_type; kwargs...)

function MeshData(rd::RefElemData, objects, 
                  cells_per_dimension_x::Int, cells_per_dimension_y::Int,
                  quadrature_type=Subtriangulation(); 
                  coordinates_min=(-1.0, -1.0), coordinates_max=(1.0, 1.0),
                  kwargs...)

    # compute intersections of curve with a background Cartesian grid.
    vx = LinRange(coordinates_min[1], coordinates_max[1], cells_per_dimension_x + 1)
    vy = LinRange(coordinates_min[2], coordinates_max[2], cells_per_dimension_y + 1)    

    return MeshData(rd, objects, vx, vy, quadrature_type; kwargs...)
end

num_cartesian_elements(md::MeshData{DIM, <:CutCellMesh}) where {DIM} = size(md.x.cartesian, 2) 
num_cut_elements(md::MeshData{DIM, <:CutCellMesh}) where {DIM} = size(md.x.cut, 2) 

# this is used in src/MeshData.jl in `getproperty` 
function num_elements(md::MeshData{DIM, <:CutCellMesh}) where {DIM}
    # the number of elements is given by the number of columns of each component of x
    return num_cartesian_elements(md) + num_cut_elements(md)
end

"""
    construct_cut_volume_quadrature(N, cutcells; target_degree = 2 * N - 1)

Constructs volume quadrature using subtriangulations of cut cells and 
Caratheodory pruning. The resulting quadrature is exact for polynomials 
of degree `target_degree`. 

We set `target_degree` to `2 * N - 1` by default, which is sufficient to 
ensure that ∫du/dx * v is integrated exactly so that integration by parts
holds under the generated cut cell quadrature. 
"""
# TODO: time this function
function construct_cut_volume_quadrature(N, cutcells, physical_frame_elements; 
                                         target_degree = 2 * N - 1)

    # We make volume quadrature exact for degree N^2 + N(N-1) + 2(N-1) polynomials on the 
    # reference element, which ensures that ∫du/dx * v * J is integrated exactly over D̂. 
    # This is because du/dx * v ∈ P^{2N-1}(D^k) and (du/dx * v) ∘ x(r,s) ∈ P^{(2N-1) * N}). 
    # Thus, du/dx * v is a degree 2N^2-N polynomial, and J is degree 2(N-1)
    # 
    # The minimum exactness of volume quadrature is degree N(N-1) + 2N-2 polynomials,
    # which ensures Qh * 1 = 0, where Qh = hybridized SBP operator.

    # Integral     ∫    v     du/dx        J     over reference element D̂
    # TODO: why can we decrease this and still see exactness?
    N_phys_frame_geo = 2 * N^2 - N + 2 * (N-1) 

    rd_tri = RefElemData(Tri(), Polynomial(MultidimensionalQuadrature()), N, 
                         quad_rule_vol=NodesAndModes.quad_nodes_tri(N_phys_frame_geo))

    Np_target = StartUpDG.Np_cut(target_degree)
    xq_pruned, yq_pruned, wJq_pruned = ntuple(_ -> zeros(Np_target, length(cutcells)), 3)

    for (e, cutcell) in enumerate(cutcells)

        xq, yq, wJq = StartUpDG.subtriangulated_cutcell_quadrature(cutcell, rd_tri)

        # if there are not enough quadrature nodes, refine the target polynomial 
        # degree until there are. 
        N_refined = N_phys_frame_geo
        while length(wJq) < Np_target
            N_refined += 1
            rd_tri_refined = RefElemData(Tri(), Polynomial(MultidimensionalQuadrature()), N, 
                                 quad_rule_vol=NodesAndModes.quad_nodes_tri(N_refined))
            xq, yq, wJq = StartUpDG.subtriangulated_cutcell_quadrature(cutcell, rd_tri_refined)
        end

        # perform Caratheodory pruning
        Vtarget = vandermonde(physical_frame_elements[e], target_degree, vec(xq), vec(yq))
        w = vec(wJq)
        w_pruned, inds = CaratheodoryPruning.caratheodory_pruning(Vtarget, w)

        # test exactness of the pruned quadrature rule if applicable
        if target_degree >= 2 * N
            V = vandermonde(physical_frame_elements[e], N, vec(xq), vec(yq))        
            @assert norm(V' * Diagonal(w) * V - V' * Diagonal(w_pruned) * V) < 100 * eps()
        end

        # the number of pruned nodes can be smaller than Np_target; if so, we simply pad with some extra points with zero weight
        # TODO: turn off the flag for trimming extra points in CaratheodoryPruning.jl
        if length(inds) < Np_target
            ind_dest = 1:length(inds)
            @. xq_pruned[ind_dest, e]  = xq[inds]
            @. yq_pruned[ind_dest, e]  = yq[inds]
            @. wJq_pruned[ind_dest, e] = w_pruned[inds]

            ind_remainder = length(inds)+1:Np_target
            xq_pruned[ind_remainder, e]  .= xq_pruned[1,e]
            yq_pruned[ind_remainder, e]  .= yq_pruned[1,e]
            wJq_pruned[ind_remainder, e] .= 0.0
        else
            @. xq_pruned[:, e] = xq[inds]
            @. yq_pruned[:, e] = yq[inds]
            @. wJq_pruned[:, e] = w_pruned[inds]
        end
    end
    return xq_pruned, yq_pruned, wJq_pruned
end

# construct interpolation points using sampling and approximate Fekete 
# points (e.g., QR-DEIM).
function construct_cut_interpolation_nodes(N, objects, physical_frame_elements)

    N_sampled = 4 * N # this is arbitrary
    num_cut_elements = length(physical_frame_elements)
    x, y = ntuple(_ -> zeros(Np_cut(N), num_cut_elements), 2)

    # Compute interpolation points on cut elements
    for e in eachindex(physical_frame_elements)

        physical_frame_element = physical_frame_elements[e]
        
        x_sampled, y_sampled = 
            generate_sampling_points(objects, physical_frame_element, 2 * Np_cut(N), N_sampled)

        V = vandermonde(physical_frame_element, N, x_sampled, y_sampled) 

        # use pivoted QR to find good interpolation points
        QRfac = qr(V', ColumnNorm())
        ids = QRfac.p[1:Np_cut(N)]

        # if the condition number of the VDM is really bad, then increase the 
        # number of sampled points. 
        condV_original = cond(V[ids, :])
        condV = condV_original
        N_sampled = 8 * N
        iter, maxiters = 0, 100
        while condV > 1e8 && iter < maxiters
            x_sampled, y_sampled = 
                generate_sampling_points(objects, physical_frame_element, 
                                         2 * Np_cut(N), N_sampled)
            V = vandermonde(physical_frame_element, N, x_sampled, y_sampled) 

            # use pivoted QR to find good interpolation points
            QRfac = qr(V', ColumnNorm())
            ids = QRfac.p[1:Np_cut(N)]
            condV = cond(V[ids, :])
            N_sampled += 2 * N
            
            if condV < 1e8
                @info "Conditioning of old VDM for element $e is $condV_original. " * 
                "After recomputing with $(length(x_sampled)) samples: $condV."
            end

            iter += 1
        end
        if iter >= maxiters
            @warn "Adaptive sampling of cut-cell element $e: " * 
                "maximum number of iterations $maxiters exceeded."
        end

        view(x, :, e) .= x_sampled[ids]
        view(y, :, e) .= y_sampled[ids]
    end
    return x, y
end

function construct_cartesian_volume_quadrature(vx, vy, region_flags, quad_rule_volume)

    rq, sq, wq = quad_rule_volume

    num_cartesian_cells = sum(region_flags .== 0)
    xq, yq, wJq = ntuple(_ -> zeros(length(wq), num_cartesian_cells), 3)

    # compute quadrature rules for the Cartesian cells
    e = 1
    for ex in axes(region_flags, 1), ey in axes(region_flags, 2)
        if is_Cartesian(region_flags[ex, ey])
            dx = vx[ex+1] - vx[ex]
            dy = vy[ey+1] - vy[ey]
            J = dx * dy / sum(wq) 
            @. xq[:, e] = dx * 0.5 * (1 + rq) + vx[ex]
            @. yq[:, e] = dy * 0.5 * (1 + sq) + vy[ey]
            @. wJq[:, e] = wq * J
            e += 1
        end
    end
    return xq, yq, wJq    
end

function construct_cartesian_surface_quadrature(vx, vy, region_flags, 
                                                quad_rule_face)

    LX, LY = (x -> x[2] - x[1]).(extrema.((vx, vy)))
    cells_per_dimension_x, cells_per_dimension_y = size(region_flags)
    r1D, w1D = quad_rule_face

    num_cartesian_cells = sum(region_flags .== 0)
    Jf = zeros(num_faces(Quad()), num_cartesian_cells)
    xf, yf, nxJ, nyJ, wf = ntuple(_ -> zeros(num_faces(Quad()) * length(r1D), 
                                              num_cartesian_cells), 5)

    # the face Jacobian involves scaling between mapped and reference domain    
    # this is precomputed here since it's needed to compute the normals

    # face 1 and 2: ±s faces
    face_ids = [eachindex(r1D); eachindex(r1D) .+ length(r1D)]
    Jf[SVector(1, 2), :] .= (LY / (cells_per_dimension_y * sum(w1D)))

    # face 3 and 4: ±r faces
    face_ids = [eachindex(r1D) .+ 2 * length(r1D); eachindex(r1D) .+ 3 * length(r1D)]
    Jf[SVector(3, 4), :] .= (LX / (cells_per_dimension_x * sum(w1D)))

    # compute face data
    e = 1
    for ex in 1:cells_per_dimension_x, ey in 1:cells_per_dimension_y    
        if is_Cartesian(region_flags[ex, ey])
            # face 1: -r
            face_index = 1
            face_ids = eachindex(r1D) .+ (face_index - 1) * length(r1D)
            xf[face_ids, e] .= vx[ex]
            yf[face_ids, e] .= map_to_interval(r1D, vy[ey], vy[ey+1])
            nxJ[face_ids, e] .= -Jf[face_index, e]
            nyJ[face_ids, e] .= zero(eltype(nyJ))
            wf[face_ids, e] .= w1D

            # face 2: +r
            face_index = 2
            face_ids = eachindex(r1D) .+ (face_index - 1) * length(r1D)
            xf[face_ids, e] .= vx[ex+1]
            yf[face_ids, e] .= map_to_interval(r1D, vy[ey], vy[ey+1])
            nxJ[face_ids, e] .= Jf[face_index, e]
            nyJ[face_ids, e] .= zero(eltype(nyJ))
            wf[face_ids, e] .= w1D

            # face 3: -s
            face_index = 3
            face_ids = eachindex(r1D) .+ (face_index - 1) * length(r1D)
            xf[face_ids, e] .= map_to_interval(r1D, vx[ex], vx[ex+1])
            yf[face_ids, e] .= vy[ey]
            nxJ[face_ids, e] .= zero(eltype(nxJ))
            nyJ[face_ids, e] .= -Jf[face_index, e]
            wf[face_ids, e] .= w1D

            # face 3: +s
            face_index = 4
            face_ids = eachindex(r1D) .+ (face_index - 1) * length(r1D)
            xf[face_ids, e] .= map_to_interval(r1D, vx[ex], vx[ex+1])
            yf[face_ids, e] .= vy[ey+1]
            nxJ[face_ids, e] .= zero(eltype(nxJ))
            nyJ[face_ids, e] .= Jf[face_index, e]
            wf[face_ids, e] .= w1D

            e = e + 1
        end
    end   

    return xf, yf, nxJ, nyJ, wf
end

"""
    construct_physical_frame_elements(region_flags, cutcells)

Computes physical frame shifting and scaling parameters from the 
vertices of cut cells and the background cell location.
"""
function construct_physical_frame_elements(region_flags, vx, vy, cutcells)

    physical_frame_elements = typeof(PhysicalFrame())[] # populate this as we iterate through cut cells
    e = 1
    s_sampling = 0:0.001:1
    for ex in axes(region_flags, 1), ey in axes(region_flags, 2)
        if StartUpDG.is_cut(region_flags[ex, ey])

            # get extremal points (vertices) of the cut cell
            cutcell = cutcells[e]
            coordinates = cutcell.(cutcell.stop_pts[1:end-1])
            stop_pts = cutcell.(cutcell.stop_pts[1:end-1])
            boundary_sample_pts = cutcell.(s_sampling)
            coordinates = [stop_pts; boundary_sample_pts]

            # store vertex nodes and coordinates of background Cartesian cell
            physical_frame_element = 
                PhysicalFrame(getindex.(coordinates, 1), getindex.(coordinates, 2), 
                              SVector(vx[ex], vx[ex+1]), SVector(vy[ey], vy[ey+1]))

            push!(physical_frame_elements, physical_frame_element)

            e += 1
        end
    end
    return physical_frame_elements
end

"""
    construct_cut_surface_quadrature(N, cutcells, quad_rule_1D = gauss_quad(0, 0, N))

Constructs cut surface quadrature using a degree `N` geometric mapping and a reference
quadrature rule `quad_rule_1D`. Returns `xf, yf, nxJ, nyJ, wf` which are vectors, and 
`face_node_indices`, which is a `Vector{Vector{Int}}` of global face node indices (which 
index into `xf.cut`, `yf.cut`, etc) for each face of each cut element.

On boundaries of cut cells, the surface quadrature is taken to be exact for degree 
`N^2 + (N-1)` polynomials. This ensures satisfaction of a weak GSBP property.
"""
function construct_cut_surface_quadrature(N, cutcells, 
                                          quad_rule_1D = gauss_quad(0, 0, N),
                                          quad_rule_boundary = get_boundary_quadrature(N)) 

    rd_line = RefElemData(Line(), N, quad_rule_vol = quad_rule_1D)
    rq_boundary, wq_boundary = quad_rule_boundary

    interp_to_cut_boundary = vandermonde(Line(), rd_line.N, rq_boundary) / rd_line.VDM

    xf, yf, nxJ, nyJ, wf = ntuple(_ -> Vector{Float64}[], 5)

    face_node_indices = Vector{UnitRange{Int64}}[]

    # recompute normals using isoparametric mapping on cut cells
    for cutcell in cutcells
        xf_element, yf_element, nxJ_element, nyJ_element, wf_element = 
            ntuple(_ -> Vector{Float64}[], 5)

        for f in 1:length(cutcell.stop_pts)-1

            is_boundary_face = 
                !(cutcell.subcurves[f] isa PathIntersections.PresetGeometries.Line)

            # switch between lower order and higher order face quadrature 
            # if it's the curved boundary, it may be higher order
            interp_to_face_quad_pts = (is_boundary_face) ? interp_to_cut_boundary : rd_line.Vq
                
            # map 1D interpolation points to curved faces
            points = map(s -> cutcell(StartUpDG.map_to_interval(s, cutcell.stop_pts[f:f+1]...)), rd_line.r)
            xf_face = getindex.(points, 1)
            yf_face = getindex.(points, 2)

            # interp face coordinates to face quad nodes
            xfq_face = interp_to_face_quad_pts * xf_face
            yfq_face = interp_to_face_quad_pts * yf_face        

            # compute tangent vector at quadrature points using 
            # the polynomial geometric mapping.
            dxdr = interp_to_face_quad_pts * rd_line.Dr * xf_face
            dydr = interp_to_face_quad_pts * rd_line.Dr * yf_face
            scaled_normal = SVector.(-dydr, dxdr)
            nxJ_face = @. getindex(scaled_normal, 1)
            nyJ_face = @. getindex(scaled_normal, 2)
            wf_face = (is_boundary_face) ? wq_boundary : rd_line.wq
            
            push!(xf_element, xfq_face)
            push!(yf_element, yfq_face)
            push!(nxJ_element, nxJ_face)
            push!(nyJ_element, nyJ_face)
            push!(wf_element, wf_face)
        end    

        # compute indices of nodes on each face before flattening
        num_nodes_per_face = length.(xf_element)
        face_node_indices_on_this_element = 
            map((length, offset) -> (1:length) .+ offset, 
                num_nodes_per_face, [0; cumsum(num_nodes_per_face[1:end-1])])
        push!(face_node_indices, face_node_indices_on_this_element)

        # flatten the collection of quadrature points for each cutcell 
        # and push them to the arrays xf, yf, ...
        map((x,y) -> push!(x, vcat(y...)), 
            (xf, yf, nxJ, nyJ, wf), 
            (xf_element, yf_element, nxJ_element, nyJ_element, wf_element))
    end

    # compute global cut face node indices
    num_nodes_per_cut_element = map(x -> maximum(maximum.(x)), face_node_indices)
    face_node_offsets = [0; cumsum(num_nodes_per_cut_element[1:end-1])]
    for e in eachindex(face_node_indices)
        for f in eachindex(face_node_indices[e])    
            face_node_indices[e][f] = face_node_indices[e][f] .+ face_node_offsets[e]
        end
    end

    return xf, yf, nxJ, nyJ, wf, face_node_indices
end

# this assumes Cartesian cells are Quads
function compute_face_centroids(vx, vy, region_flags, 
                                cutcells::AbstractVector{<:PathIntersections.PiecewiseCurve})

    num_cartesian_cells = sum(region_flags .== 0)   
    face_centroids_cartesian = 
        zeros(SVector{2, Float64}, num_faces(Quad()), num_cartesian_cells)
    
    num_faces_per_cut_cell = map(cutcell->length(cutcell.subcurves), cutcells)        
    face_centroids_cut = [zeros(SVector{2, Float64}, num_faces_per_cut_cell[e]) 
                            for e in eachindex(cutcells)]

    # calculate Cartesian centroids
    e = 1
    for ex in axes(region_flags, 1), ey in axes(region_flags, 2)
        if is_Cartesian(region_flags[ex, ey])
            # quad faces are -r, +r, -s, +s
            face_centroids_cartesian[1, e] = SVector(vx[ex],   0.5 * (vy[ey+1] + vy[ey]))
            face_centroids_cartesian[2, e] = SVector(vx[ex+1], 0.5 * (vy[ey+1] + vy[ey]))
            face_centroids_cartesian[3, e] = SVector(0.5 * (vx[ex] + vx[ex+1]), vy[ey])
            face_centroids_cartesian[4, e] = SVector(0.5 * (vx[ex] + vx[ex+1]), vy[ey+1])
            e += 1
        end
    end

    for (e, cutcell) in enumerate(cutcells)
        for (f, face_parametrization) in enumerate(cutcell.subcurves)
            s_midpoint = 0.5 * sum(cutcell.sub_bounds[f])
            x, y = face_parametrization(s_midpoint)
            face_centroids_cut[e][f] = SVector(x, y)
        end
    end

    return face_centroids_cartesian, face_centroids_cut
end

num_faces(cutcell::PathIntersections.PiecewiseCurve) = length(cutcell.subcurves)

function connect_mesh(face_centroids, region_flags, 
                      cutcells::Vector{<:PathIntersections.PiecewiseCurve}; 
                      tol = 1e2 * eps())

    num_cartesian_cells = sum(region_flags .== 0)
    cut_faces_per_cell = num_faces.(cutcells)
    cut_face_offsets = [0; cumsum(cut_faces_per_cell)[1:end-1]] 
    num_cut_faces = sum(cut_faces_per_cell)
    num_total_faces = num_cut_faces + num_faces(Quad()) * num_cartesian_cells

    # element_indices[ex, ey] returns the global (flattened) element index into 
    # the arrays `xf.cartesian[:, e]` or `xf.cut[:, e]`
    element_indices = compute_element_indices(region_flags) 

    # compute face centroids for making face matches
    face_centroids_x, face_centroids_y = face_centroids
    
    # To determine face-to-face matches, we work with each background Cartesian element 
    # and search through the 4 neighboring background Cartesian elements for a match in 
    # the face centroids of the current cell and the face centroids of its neighbors.     
    # NOTE: this works because the cut cells can only have Cartesian neighbors across
    # flat-sided faces. It wouldn't work for meshes with curved interior interfaces.    

    FToF = collect(1:num_total_faces)
    for ex in axes(region_flags, 1), ey in axes(region_flags, 2)

        e = element_indices[ex, ey]

        # Determine face indices of current cell. The face indices are determined 
        # from the flattened (e.g., not ex, ey) element ordering. 
        # NOTE: we search for matches between all faces of `e` and all faces of 
        # `e_nbr` because the ordering of faces is different for cut elements
        # and Cartesian elements. 
        if is_Cartesian(region_flags[ex, ey])
            face_ids = (1:num_faces(Quad())) .+ (e-1) * num_faces(Quad())
        elseif is_cut(region_flags[ex, ey])
            face_ids = (1:cut_faces_per_cell[e]) .+ cut_face_offsets[e]

            # we offset by the number of Cartesian faces so we can index globally
            # into the arrays `face_centroids_x`, `face_centroids_y`.
            face_ids = face_ids .+ length(face_centroids_x.cartesian) 
        end

        if is_inside_domain(ex, ey, region_flags)

            for f in 1:4 # each Cartesian background element has 4 neighbors

                ex_nbr, ey_nbr = neighbor_across_face(f, ex, ey)
                if is_inside_domain(ex_nbr, ey_nbr, region_flags)
                    e_nbr = element_indices[ex_nbr, ey_nbr]

                    # determine face indices of neighboring cells
                    if is_Cartesian(region_flags[ex_nbr, ey_nbr])                
                        nbr_face_ids = (1:num_faces(Quad())) .+ (e_nbr-1) * num_faces(Quad())
                    elseif is_cut(region_flags[ex_nbr, ey_nbr])
                        nbr_face_ids = (1:cut_faces_per_cell[e_nbr]) .+ cut_face_offsets[e_nbr]

                        # we offset by the number of Cartesian faces so we can index globally
                        # into the arrays `face_centroids_x`, `face_centroids_y`.
                        nbr_face_ids = nbr_face_ids .+ length(face_centroids_x.cartesian)
                    end

                    # check for matches in face and neighbor face centroids.
                    # note: we index into the global `face_centroids_x/y` container 
                    # rather than the `.cut` or `.cartesian subarrays`.
                    for i in face_ids
                        xy = SVector(face_centroids_x[i], face_centroids_y[i])
                        for j in nbr_face_ids
                            xy_nbr = SVector(face_centroids_x[j], face_centroids_y[j])
                            if norm(xy - xy_nbr) < tol * max(1, norm(xy), norm(xy_nbr))
                                FToF[i] = j  
                            end
                        end
                    end

                end # if enbr is_inside_domain
            end
        end # if e is_inside_domain
    end

    return FToF
end

# If we wish to exactly integrate ∫ u * v * (nx * Jf) and ∫ u * v * (ny * Jf), 
# the reference quadrature must be exact for degree 2 * N^2 + (N-1) polynomials
# on the reference line. 
#
# Note that we don't perform Caratheodory pruning for boundary faces at the moment.
get_boundary_quadrature(rd::RefElemData) = get_boundary_quadrature(rd.N)
function get_boundary_quadrature(N)
    # N_boundary = ceil(Int, (N^2 + (N-1) - 1) / 2) # to ensure weak SBP
    N_boundary = ceil(Int, (2 * N^2 + (N - 1) - 1) / 2) # to ensure GSBP property 
    return gauss_quad(0, 0, N_boundary)
end

"""
    function MeshData(rd, geometry, vxyz...)

Creates a cut-cell mesh where the boundary is given by `geometry`, which should be a 
tuple of functions. These functions can be generated using PathIntersections.PresetGeometries, 
for example:
```julia
julia> geometry = (PresetGeometries.Circle(R=0.33, x0=0, y0=0), )
```
Here, `coordinates_min`, `coordinates_max` contain `(smallest value of x, smallest value of y)` 
and `(largest value of x, largest value of y)`, and `cells_per_dimension_x/y` is the number of 
Cartesian grid cells placed along each dimension. 
"""

function MeshData(rd::RefElemData, objects, 
                  vx::AbstractVector, vy::AbstractVector, 
                  quadrature_type::Subtriangulation; 
                  quad_rule_face=get_1d_quadrature(rd), 
                  quad_rule_boundary=get_boundary_quadrature(rd),
                  cut_quadrature_target_degree = 2 * rd.N-1,
                  precompute_operators=false)

    cells_per_dimension_x = length(vx) - 1
    cells_per_dimension_y = length(vy) - 1
    
    LX, LY = (x -> x[2] - x[1]).(extrema.((vx, vy)))

    # compute mesh intersections and cut cell elements.
    # `regions` is a matrix of dimensions `(cells_per_dimension_x, cells_per_dimension_y)` with 3 values:
    #   *  1: cut cell
    #   *  0: Cartesian cell
    #   * -1: excluded cells (in the cut-out region)
    region_flags, cutcells = StartUpDG.calculate_cutcells(vx, vy, objects)

    # pack useful cut cell information together. 
    num_cartesian_cells = sum(region_flags .== 0)
    num_cut_cells = sum(region_flags .> 0)

    N = rd.N
    
    ####################################################
    #          Construct cut cells stuff               # 
    ####################################################

    physical_frame_elements = 
        construct_physical_frame_elements(region_flags, vx, vy, cutcells)

    x_cut, y_cut = 
        construct_cut_interpolation_nodes(N, objects, physical_frame_elements)

    xq_cut, yq_cut, wJq_cut = 
        construct_cut_volume_quadrature(N, cutcells, physical_frame_elements; 
                                        target_degree = cut_quadrature_target_degree)

    # On the curved cut boundary, we want to exactly integrate degree N^2 + (N-1) 
    # polynomials since our target integrand is nxJ * u (nxJ is degree N-1 on the 
    # reference face, and u is degree N on the physical face, but mapped back to 
    # the reference face it's N^2). quad_rule_boundary defaults to this degree.
    xf_cut, yf_cut, nxJ_cut, nyJ_cut, wf_cut, cut_face_node_indices = 
        construct_cut_surface_quadrature(N, cutcells, quad_rule_face, quad_rule_boundary)

    ####################################################
    #          Construct Cartesian cells               # 
    ####################################################

    xf_cartesian, yf_cartesian, nxJ_cartesian, nyJ_cartesian, wf_cartesian = 
        construct_cartesian_surface_quadrature(vx, vy, region_flags, quad_rule_face)    

    xq_cartesian, yq_cartesian, wJq_cartesian = 
        construct_cartesian_volume_quadrature(vx, vy, region_flags, 
                                             (rd.rq, rd.sq, rd.wq))

    # reuse quadrature mapping to construct interpolation nodes                                             
    x_cartesian, y_cartesian, _ = 
        construct_cartesian_volume_quadrature(vx, vy, region_flags, 
                                             (rd.r, rd.s, ones(size(rd.r))))
                                             
    ####################################################
    #      Assemble combined cut/Cartesian arrays      # 
    ####################################################

    x, y = map((x, y) -> NamedArrayPartition((; cartesian=x, cut=y)), 
               (x_cartesian, y_cartesian), (x_cut, y_cut))    

    xq, yq, wJq = map((x, y) -> NamedArrayPartition((; cartesian=x, cut=y)), 
                      (xq_cartesian, yq_cartesian, wJq_cartesian), 
                      (xq_cut, yq_cut, wJq_cut))                                                                

    xf, yf, nxJ, nyJ, wf = 
        map((x, y) -> NamedArrayPartition((; cartesian=x, cut=y)), 
            (xf_cartesian, yf_cartesian, nxJ_cartesian, nyJ_cartesian, wf_cartesian),
            map(x -> vcat(x...), (xf_cut, yf_cut, nxJ_cut, nyJ_cut, wf_cut)))

    Jf = @. sqrt(nxJ^2 + nyJ^2)
    wJf = @. wf * Jf
            
    ####################################################
    #                Mesh connectivity                 # 
    ####################################################

    # note: face_centroids_cartesian has dims [4, num_elements]
    #       face_centroids_cut has dims [num_elements][num_faces_on_this_element]
    face_centroids_cartesian, face_centroids_cut = 
            StartUpDG.compute_face_centroids(vx, vy, region_flags, cutcells)

    face_centroids_x = NamedArrayPartition(cartesian = getindex.(face_centroids_cartesian, 1), 
                                        cut = getindex.(vcat(face_centroids_cut...), 1))
    face_centroids_y = NamedArrayPartition(cartesian = getindex.(face_centroids_cartesian, 2), 
                                        cut = getindex.(vcat(face_centroids_cut...), 2))
    face_centroids = (face_centroids_x, face_centroids_y)
    FToF = StartUpDG.connect_mesh(face_centroids, region_flags, cutcells)  

    # create arrays of node indices on each face, e.g., 
    # xf.cut[cut_face_node_indices_by_face[f]] gives nodes on face f
    cartesian_nodes_per_face = rd.Nfq ÷ rd.num_faces
    cartesian_face_node_indices_by_face = 
        [(1:cartesian_nodes_per_face) .+ (f - 1) * cartesian_nodes_per_face
            for f in 1:num_cartesian_cells * rd.num_faces]
    cut_face_node_indices_by_face = vcat(cut_face_node_indices...)

    # create global index vector (note that we offset cut cell indices 
    # by num_cartesian nodes).
    face_node_indices_by_face = vcat(cartesian_face_node_indices_by_face, 
                                    map(x -> x .+ (num_cartesian_cells * rd.Nfq), 
                                        cut_face_node_indices_by_face))

    # create node mappings
    mapM_cartesian = num_cartesian_cells > 0 ? 
        reshape(eachindex(xf.cartesian), :, num_cartesian_cells) : Int64[]
    mapM = NamedArrayPartition(cartesian = collect(mapM_cartesian), 
                            cut = collect(eachindex(xf.cut) .+ length(mapM_cartesian)))
    mapP = copy(mapM)

    # allocate a vector for points on non-boundary faces
    p = zeros(Int, length(first(quad_rule_face)))
    for f in eachindex(FToF)
        # if it's not a boundary face
        if f !== FToF[f]
            idM = view(mapM, face_node_indices_by_face[f])
            idP = view(mapM, face_node_indices_by_face[FToF[f]])
            xyzM = (view(xf, idM), view(yf, idM))
            xyzP = (view(xf, idP), view(yf, idP))
            StartUpDG.match_coordinate_vectors!(p, xyzM, xyzP)
            for (i, id) in enumerate(p)
                mapP[idM[id]] = idP[i]
            end
        end
    end
    mapB = findall(vec(mapM) .== vec(mapP)) # determine boundary nodes

    # default to non-periodic 
    is_periodic = (false, false)   
   
    # compute mapping from linear element indices to Cartesian element indices
    cartesian_to_linear_element_indices = compute_element_indices(region_flags)
    linear_to_cartesian_element_indices = (; cut=zeros(SVector{2, Int}, num_cut_cells), 
                                             cartesian=zeros(SVector{2, Int}, num_cartesian_cells))
    for ex in 1:cells_per_dimension_x, ey in 1:cells_per_dimension_y
        e = cartesian_to_linear_element_indices[ex, ey]
        if is_cut(region_flags[ex, ey])
            linear_to_cartesian_element_indices.cut[e] = SVector(ex, ey)
        elseif is_Cartesian(region_flags[ex, ey])
            linear_to_cartesian_element_indices.cartesian[e] = SVector(ex, ey)
        end
    end

    cut_cell_data = (; cutcells, region_flags, wJf,
                       cartesian_to_linear_element_indices, 
                       linear_to_cartesian_element_indices,                       
                       vxyz=(vx, vy), 
                       cells_per_dimension=(cells_per_dimension_x, 
                                            cells_per_dimension_y), 
                       cut_face_node_indices_by_elem_by_face = cut_face_node_indices,
                    ) # background Cartesian grid info

    # get flattened indices of cut face nodes. note that 
    #   `cut_face_node_indices = [[face_1_indices, face_2_indices, ...], ...]`
    cut_face_node_indices_by_cell = 
        map(x -> UnitRange(first(x[1]), last(x[end])), cut_face_node_indices)

    if precompute_operators == true

        # precompute cut-cell operators and store them in the `md.mesh_type.cut_cell_operators` field.
        differentiation_matrices = Tuple{Matrix{eltype(x)}, Matrix{eltype(x)}}[]
        volume_interpolation_matrices, face_interpolation_matrices = 
            ntuple(_ -> Matrix{eltype(x)}[], 3)

        mass_matrices, lift_matrices, projection_matrices = ntuple(_ -> Matrix{eltype(x)}[], 3)
        for (e, elem) in enumerate(physical_frame_elements)

            VDM = vandermonde(elem, rd.N, x.cut[:, e], y.cut[:, e])
            Vq, Vrq, Vsq = map(A -> A / VDM, 
                               basis(elem, rd.N, view(xq.cut, :, e), view(yq.cut, :, e)))
        
            M  = Vq' * Diagonal(wJq.cut[:, e]) * Vq
            Qr = Vq' * Diagonal(wJq.cut[:, e]) * Vrq
            Qs = Vq' * Diagonal(wJq.cut[:, e]) * Vsq    
            Dx_e, Dy_e = M \ Qr, M \ Qs
            
            xf_e = xf.cut[cut_face_node_indices_by_cell[e]]
            yf_e = yf.cut[cut_face_node_indices_by_cell[e]]
            Vf = vandermonde(elem, rd.N, xf_e, yf_e) / VDM

            # don't include jacobian scaling in LIFT matrix (for consistency 
            # with the Cartesian mesh)
            wf = wJf.cut[cut_face_node_indices_by_cell[e]] ./ 
                    Jf.cut[cut_face_node_indices_by_cell[e]]

            push!(volume_interpolation_matrices, Vq)
            push!(face_interpolation_matrices, Vf)
            push!(differentiation_matrices, (Dx_e, Dy_e))
            push!(mass_matrices, M)
            push!(lift_matrices, M \ (Vf' * Diagonal(wf)))
            push!(projection_matrices, M \ (Vq' * Diagonal(wJq.cut[:, e])))
        end
        cut_cell_operators = (; volume_interpolation_matrices, 
                                face_interpolation_matrices, 
                                differentiation_matrices, 
                                mass_matrices, lift_matrices, 
                                projection_matrices)

    else

        cut_cell_operators = nothing
    end
               

    ################################################
    #       Construct volume geometric terms       #
    ################################################

    rxJ_cartesian = LX / (2 * cells_per_dimension_x)
    syJ_cartesian = LY / (2 * cells_per_dimension_y)
    J_cartesian = (LX / cells_per_dimension_x) * (LY / cells_per_dimension_y) / 4                                               

    # Note: the volume Jacobian for cut elements is 1 since the "reference element" is the 
    # cut element itself. Similarly, geometric terms should be 1 since `basis` computes 
    # physical derivatives accounting for element scaling

    # Note: we use FillArrays.Fill instead of FillArrays.Ones and FillArrays.Zeros because 
    # we store `rxJ, sxJ, ryJ, syJ` in a single SMatrix, which assumes one homogeneous 
    # type for all the entries. Since FillArrays.Ones/Zeros are distinct types, using them
    # results in a type instability when accessing entries of md.rstxyzJ
    rxJ = NamedArrayPartition(cartesian=Fill(rxJ_cartesian, rd.Np, num_cartesian_cells), 
                              cut=Fill(1.0, Np_cut(rd.N), num_cut_cells))
    syJ = NamedArrayPartition(cartesian=Fill(syJ_cartesian, rd.Np, num_cartesian_cells), 
                              cut=Fill(1.0, Np_cut(rd.N), num_cut_cells))
    sxJ, ryJ = ntuple(_ -> NamedArrayPartition(cartesian=Fill(0.0, rd.Np, num_cartesian_cells), 
                                               cut=Fill(0.0, Np_cut(rd.N), num_cut_cells)), 2) 
    J = NamedArrayPartition(cartesian = Fill(J_cartesian, rd.Np, num_cartesian_cells), 
                            cut = Fill(1.0, Np_cut(rd.N), num_cut_cells))

    # pack geometric terms together                            
    rstxyzJ = SMatrix{2, 2, typeof(rxJ), 4}(rxJ, sxJ, ryJ, syJ) 
    
    return MeshData(CutCellMesh(physical_frame_elements, cut_face_node_indices_by_cell, objects, 
                                cut_cell_operators, cut_cell_data, quadrature_type), 
                    FToF, (x, y), (xf, yf), (xq, yq), wJq, 
                    mapM, mapP, mapB, rstxyzJ, J, (nxJ, nyJ), Jf, is_periodic)

end

"""
    hybridized_SBP_operators(md::MeshData{2, <:CutCellMesh})

This constructs hybridized SBP operators using the approach taken in Chan (2019), 
"Skew-Symmetric Entropy Stable Modal Discontinuous Galerkin Formulations". 
[https://doi.org/10.1007/s10915-019-01026-w](https://doi.org/10.1007/s10915-019-01026-w)

This function returns `hybridized_operators::Vector{Tuple{<:Matrix, <:Matrix}}` and 
`project_and_interp_operators, projection_operators, interpolation_operators`, which are 
all `Vector{<:Matrix}`, where each entry corresponds to a cut element.
"""
function hybridized_SBP_operators(md::MeshData{2, <:CutCellMesh})
    mt = md.mesh_type
    (; volume_interpolation_matrices, 
       face_interpolation_matrices, 
       differentiation_matrices, 
       mass_matrices) = mt.cut_cell_operators
    (; wJf) = mt.cut_cell_data
    wf = wJf ./ md.Jf
    (; x, nxJ, nyJ) = md

    hybridized_operators = Tuple{Matrix{eltype(x)}, Matrix{eltype(x)}}[]
    interpolation_operators, projection_operators, project_and_interp_operators = 
        ntuple(_ -> Matrix{eltype(x)}[], 3)
    for e in eachindex(differentiation_matrices)           
        Dx, Dy = differentiation_matrices[e]
        Vq = volume_interpolation_matrices[e]
        Vf = face_interpolation_matrices[e]
        M = mass_matrices[e]

        Pq = M \ (Vq' * Diagonal(md.wJq.cut[:, e]))
        E = Vf * Pq

        fids = mt.cut_face_nodes[e]
        Bx = Diagonal(wf.cut[fids] .* nxJ.cut[fids])
        By = Diagonal(wf.cut[fids] .* nyJ.cut[fids])

        Qx = Pq' * M * Dx * Pq
        Qy = Pq' * M * Dy * Pq
        Qxh = 0.5 * [Qx - Qx' E' * Bx;
                    -Bx * E   Bx     ]
        Qyh = 0.5 * [Qy - Qy' E' * By;
                    -By * E   By      ]

        # interpolation matrices
        Vh = [Vq; Vf]
        VhP = Vh * Pq
        Ph = M \ Vh'

        push!(hybridized_operators, (Qxh, Qyh))
        push!(interpolation_operators, Vh)
        push!(projection_operators, Ph)
        push!(project_and_interp_operators, VhP)
    end
    
    return hybridized_operators, project_and_interp_operators, 
           projection_operators, interpolation_operators 
end