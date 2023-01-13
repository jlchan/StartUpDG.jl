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

The field `curves` contains a tuple of the curves used to define the cut region.

The field `cut_cell_operators` contains optionally precomputed operators (mass, differntiation, 
face interpolation, and lifting operators). 

The field `cut_cell_data` contains additional data from PathIntersections.
"""
struct CutCellMesh{T1, T2, T3, T4, T5}
    physical_frame_elements::T1
    cut_face_nodes::T2
    curves::T3
    cut_cell_operators::T4
    cut_cell_data::T5
end

# TODO: add isoparametric cut cell mesh with positive quadrature points
# # This mesh type has a polynomial representation of curves, so we don't store the curve info
# struct IsoparametricCutCellMesh{T1, T2, T3, T4}
#     physical_frame_elements::T1
#     cut_face_nodes::T2
#     cut_cell_operators::T3
#     cut_cell_data::T4
# end


function Base.show(io::IO, ::MIME"text/plain", md::MeshData{DIM, <:CutCellMesh}) where {DIM}
    @nospecialize md
    print(io,"Cut-cell MeshData of dimension $DIM with $(md.num_elements) elements " * 
             "($(num_cartesian_elements(md)) Cartesian, $(num_cut_elements(md)) cut)")
end

# maps x ∈ [-1,1] to [a,b]
map_to_interval(x, a, b) = a + (b-a) * 0.5 * (1 + x)

function count_cut_faces(cutcells)
    num_cut_faces = zeros(Int, length(cutcells))
    for e in eachindex(cutcells)
        curve = cutcells[e]
        stop_points = curve.stop_pts        
        num_cut_faces[e] = length(stop_points) - 1
    end
    return num_cut_faces
end

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

function compute_face_centroids(rd, xf, yf, cutcell_data)

    @unpack region_flags, cut_faces_per_cell, cut_face_offsets = cutcell_data
    num_cut_cells = length(cut_faces_per_cell)
    num_cartesian_cells = sum(region_flags .== 0)
    num_cut_faces = sum(cut_faces_per_cell)

    num_points_per_face = length(rd.rf) ÷ num_faces(rd.element_type)
    
    face_centroids_x = NamedArrayPartition(cartesian=zeros(num_faces(rd.element_type), num_cartesian_cells), 
                                           cut=zeros(num_cut_faces))
    face_centroids_y = similar(face_centroids_x)

    for e in 1:num_cartesian_cells
        xf_element = reshape(view(xf.cartesian, :, e), num_points_per_face, num_faces(rd.element_type))
        yf_element = reshape(view(yf.cartesian, :, e), num_points_per_face, num_faces(rd.element_type))
        view(face_centroids_x.cartesian, :, e) .= vec(sum(xf_element, dims=1) / num_points_per_face)
        view(face_centroids_y.cartesian, :, e) .= vec(sum(yf_element, dims=1) / num_points_per_face)
    end

    for e in 1:num_cut_cells    
        face_node_ids = (1:(num_points_per_face * cut_faces_per_cell[e])) .+ cut_face_offsets[e] * num_points_per_face
        xf_element = reshape(view(xf.cut, face_node_ids), num_points_per_face, cut_faces_per_cell[e])
        yf_element = reshape(view(yf.cut, face_node_ids), num_points_per_face, cut_faces_per_cell[e])

        face_ids = (1:cut_faces_per_cell[e]) .+ cut_face_offsets[e]
        view(face_centroids_x.cut, face_ids) .= vec(sum(xf_element, dims=1) / num_points_per_face)
        view(face_centroids_y.cut, face_ids) .= vec(sum(yf_element, dims=1) / num_points_per_face)
    end

    return face_centroids_x, face_centroids_y
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

# generates at least Np_target sampling points within a cut cell defined by `curve`
# returns both x_sampled, y_sampled (physical points inside the cut cell), as well as 
# r_sampled, y_sampled (reference points which correspond to x_sampled, y_sampled).
function generate_sampling_points(curves, elem, rd, Np_target; N_sampled = 4 * rd.N)

    r_sampled, s_sampled = equi_nodes(rd.element_type, N_sampled) # oversampled nodes

    # map sampled points to the background Cartesian cell
    x_sampled, y_sampled = map_nodes_to_background_cell(elem, r_sampled, s_sampled)
    is_in_element = .!(is_contained.(curves, zip(x_sampled, y_sampled)))

    # increase number of background points until we are left with `Np_target` sampling points 
    while sum(is_in_element) < Np_target
        is_in_element = is_contained.(curves, zip(x_sampled, y_sampled)) .== false
        if sum(is_in_element) < Np_target
            N_sampled += rd.N
            r_sampled, s_sampled = equi_nodes(rd.element_type, N_sampled) # oversampled nodes
            x_sampled, y_sampled = map_nodes_to_background_cell(elem, r_sampled, s_sampled)
        end
    end

    ids_in_element = findall(is_in_element)

    return x_sampled[ids_in_element], y_sampled[ids_in_element]
end

# returns points (xf, yf), scaled normals (nxJ, nyJ), and face Jacobian (Jf) 
# for a curve returned from PathIntersections. 
# `out` should hold `xf, yf, nxJ, nyJ, Jf = ntuple(_ -> similar(points, (length(points), num_faces)), 5)`
function map_points_to_cut_cell_faces(points, curve)    
    num_faces = length(curve.subcurves)
    out = ntuple(_ -> similar(points, (length(points) * num_faces)), 5)
    return map_points_to_cut_cell_faces!(out, points, curve)
end

function map_points_to_cut_cell_faces!(out, points, curve)
    xf, yf, nxJ, nyJ, Jf = out
    stop_points = curve.stop_pts
    r1D = points
    for f in 1:length(stop_points)-1
        for i in eachindex(r1D)
            fid = i + (f-1) * length(r1D)

            s = map_to_interval(r1D[i], stop_points[f], stop_points[f+1])
            
            # compute tangent vector at a node, face Jacobian, and normals
            x_node, y_node = curve(s)
            xf[fid], yf[fid] = x_node, y_node                    
            tangent_vector = PathIntersections.ForwardDiff.derivative(curve, s)

            # the face Jacobian involves scaling between mapped and reference face
            # reference face = [-1, 1]
            scaling = (stop_points[f+1] - stop_points[f]) / 2
            Jf[fid] = norm(tangent_vector) * scaling

            # we have to flip the sign to get the outward normal. 
            # note: we compute the scaled normal nxJ for consistency with other meshes. 
            normal_node = SVector{2}(tangent_vector[2], -tangent_vector[1])
            nxJ[fid], nyJ[fid] = (-normal_node / norm(normal_node)) .* Jf[fid] 
        end
    end
    return vec.((xf, yf, nxJ, nyJ, Jf))
end


# Computes face geometric terms from a RefElemData, `quad_rule_face = (r1D, w1D)`, 
# the vectors of the 1D vertex nodes `vx` and `vy`, and named tuple 
# `cutcell_data is a NamedTuple containing `curves`, `region_flags`, `stop_pts``, `cutcells`. 
function compute_geometric_data(rd::RefElemData{2, Quad}, quad_rule_face, 
                                vx, vy, cutcell_data; tol=100 * eps())

    # domain size and reference face weights
    cells_per_dimension_x, cells_per_dimension_y = length(vx) - 1, length(vy) - 1
    LX, LY = (x -> x[2] - x[1]).(extrema.((vx, vy)))
    
    r1D, w1D = quad_rule_face

    @unpack curves, region_flags, cutcells, cut_faces_per_cell = cutcell_data

    # count number of cells and cut face nodes
    num_cartesian_cells = sum(region_flags .== 0)
    num_cut_cells = sum(region_flags .== 1) 
    nodes_per_face = length(r1D)
    num_cut_face_nodes = nodes_per_face * sum(cut_faces_per_cell)

    # compute face data
    xf, yf, nxJ, nyJ, Jf = ntuple(_ -> NamedArrayPartition(cartesian=zeros(rd.Nfq, num_cartesian_cells), 
                                                           cut=zeros(num_cut_face_nodes)), 5)

    # the face Jacobian involves scaling between mapped and reference domain    
    # this is precomputed here since it's needed to compute the normals
    face_ids_left_right = 1:(length(rd.rf) ÷ 2)
    face_ids_top_bottom = ((length(rd.rf) ÷ 2) + 1):length(rd.rf)
    Jf.cartesian[face_ids_top_bottom, :] .= LX / (cells_per_dimension_x * sum(w1D))
    Jf.cartesian[face_ids_left_right, :] .= LY / (cells_per_dimension_y * sum(w1D))

    # compute face data
    e = 1
    for ex in 1:cells_per_dimension_x, ey in 1:cells_per_dimension_y    
        if is_Cartesian(region_flags[ex, ey])
            vx_element = SVector(vx[ex], vx[ex + 1], vx[ex], vx[ex + 1])
            vy_element = SVector(vy[ey], vy[ey], vy[ey + 1], vy[ey + 1])
            x_element, y_element = map(x -> rd.V1 * x, (vx_element, vy_element))
            mul!(view(xf.cartesian, :, e), rd.Vf, x_element)
            mul!(view(yf.cartesian, :, e), rd.Vf, y_element)
            view(nxJ.cartesian, :, e) .= rd.nrJ .* view(Jf.cartesian, :, e)
            view(nyJ.cartesian, :, e) .= rd.nsJ .* view(Jf.cartesian, :, e)
            e = e + 1
        end
    end    

    e = 1
    offset = 0
    for (e, curve) in enumerate(cutcells)
        # map 1D quadrature points to faces
        num_cut_faces = length(curve.subcurves)
        fids = (1:length(r1D) * num_cut_faces) .+ offset
        out = map(x->view(x, fids), (xf.cut, yf.cut, nxJ.cut, nyJ.cut, Jf.cut))
        map_points_to_cut_cell_faces!(out, r1D, curve)
        offset += length(fids)             
    end                                        

    # compute face points + shifting/scaling coefficients for physical frame cut elements.
    physical_frame_elements = PhysicalFrame{2}[] # populate this as we iterate through cut cells

    # store cut-cell scaling/shifting coefficients
    @unpack cut_faces_per_cell, cut_face_offsets = cutcell_data
    num_points_per_face = length(r1D)

    e = 1
    for ex in 1:cells_per_dimension_x, ey in 1:cells_per_dimension_y    
        if is_cut(region_flags[ex, ey])

            # here, we evaluate a PhysicalFrame basis by shifting and scaling the 
            # coordinates on an element back to the reference element [-1, 1]^2.
            cut_face_node_ids = (1:num_points_per_face * cut_faces_per_cell[e]) .+ 
                                 num_points_per_face * cut_face_offsets[e]

            # store face nodes (extremal) and coordinates of background Cartesian cell
            physical_frame_element = 
                PhysicalFrame(xf.cut[cut_face_node_ids], yf.cut[cut_face_node_ids], 
                              SVector(vx[ex], vx[ex+1]), SVector(vy[ey], vy[ey+1]))

            push!(physical_frame_elements, physical_frame_element)

            e += 1
        end
    end

    # interpolation points
    x, y = ntuple(_ -> NamedArrayPartition(cartesian=zeros(rd.Np, num_cartesian_cells), 
                                           cut=zeros(Np_cut(rd.N), num_cut_cells)), 2)

    # compute interpolation points on cartesian elements
    e = 1
    for ex in 1:cells_per_dimension_x, ey in 1:cells_per_dimension_y    
        if is_Cartesian(region_flags[ex, ey])
            vx_element = SVector(vx[ex], vx[ex + 1], vx[ex], vx[ex + 1])
            vy_element = SVector(vy[ey], vy[ey], vy[ey + 1], vy[ey + 1])
            x_element, y_element = map(x -> rd.V1 * x, (vx_element, vy_element))
            view(x.cartesian, :, e) .= x_element
            view(y.cartesian, :, e) .= y_element
            e = e + 1
        end
    end
     
    # Compute interpolation points on cut elements
    for e in eachindex(physical_frame_elements)

        physical_frame_element = physical_frame_elements[e]
        
        x_sampled, y_sampled = 
            generate_sampling_points(curves, physical_frame_element, rd, 2 * Np_cut(rd.N))
        V = vandermonde(physical_frame_element, rd.N, x_sampled, y_sampled) 

        # use pivoted QR to find good interpolation points
        QRfac = qr(V', ColumnNorm())
        ids = QRfac.p[1:Np_cut(rd.N)]

        # if the condition number of the VDM is really bad, then increase the 
        # number of sampled points. 
        if cond(V[ids,:]) > 1e8
            @warn "Conditioning of VDM for element $e is $(cond(V[ids,:]));" * 
                    "recomputing with a finer set of samples."
            x_sampled, y_sampled = 
                generate_sampling_points(curves, physical_frame_element, rd, 2 * Np_cut(rd.N); 
                                            N_sampled = 100)
            V = vandermonde(physical_frame_element, rd.N, x_sampled, y_sampled) 
        end

        view(x.cut, :, e) .= x_sampled[ids]
        view(y.cut, :, e) .= y_sampled[ids]
    end

    # volume geometric terms
    rxJ_cartesian = LX / (2 * cells_per_dimension_x)
    syJ_cartesian = LY / (2 * cells_per_dimension_y)
    J_cartesian = (LX / cells_per_dimension_x) * (LY / cells_per_dimension_y) / 4                                               

    # Note: the volume Jacobian for cut elements is 1 since the "reference element" is the 
    # cut element itself. Similarly, geometric terms should be 1 since `basis` computes 
    # physical derivatives accounting for element scaling
    rxJ = NamedArrayPartition(cartesian=Fill(rxJ_cartesian, rd.Np, num_cartesian_cells), 
                              cut=Ones(Np_cut(rd.N), num_cut_cells))
    syJ = NamedArrayPartition(cartesian=Fill(syJ_cartesian, rd.Np, num_cartesian_cells), 
                              cut=Ones(Np_cut(rd.N), num_cut_cells))
    sxJ, ryJ = ntuple(_ -> NamedArrayPartition(cartesian=Zeros(rd.Np, num_cartesian_cells), 
                                               cut=Zeros(Np_cut(rd.N), num_cut_cells)), 2) 
    J = NamedArrayPartition(cartesian = Fill(J_cartesian, rd.Np, num_cartesian_cells), 
                            cut = Ones(Np_cut(rd.N), num_cut_cells))
    rstxyzJ = SMatrix{2, 2}(rxJ, sxJ, ryJ, syJ) # pack geometric terms together

    return physical_frame_elements, x, y, rstxyzJ, J, xf, yf, nxJ, nyJ, Jf
end

"""
    connect_mesh(rd, face_centroids, region_flags, cutcells; tol = 1e2 * eps())
    
Connects faces of a cut mesh to each other, returns `FToF` such that face 
`f` is connected to `FToF[f]`. 

Inputs:
- rd::RefElemData
- face_centroids = (face_centroids_x, face_centroids_y), where `face_centroids_x/y` 
                    are vectors of coordinates of face centroids
- `region_flags`, `cutcells` are return arguments from `PathIntersections.define_regions`
The keyword argument `tol` is the tolerance for matches between face centroids. 
"""    
function connect_mesh(rd, face_centroids, cutcell_data; tol = 1e2 * eps())

    @unpack region_flags, cut_faces_per_cell, cut_face_offsets = cutcell_data

    cells_per_dimension_x, cells_per_dimension_y = size(region_flags)
    num_cartesian_cells = sum(region_flags .== 0)
    num_cut_faces = sum(cut_faces_per_cell)
    num_total_faces = num_cut_faces + num_faces(rd.element_type) * num_cartesian_cells

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
    for ex in 1:cells_per_dimension_x, ey in 1:cells_per_dimension_y

        e = element_indices[ex, ey]

        # Determine face indices of current cell. The face indices are determined 
        # from the flattened (e.g., not ex, ey) element ordering. 
        # NOTE: we search for matches between all faces of `e` and all faces of 
        # `e_nbr` because the ordering of faces is different for cut elements
        # and Cartesian elements. 
        if is_Cartesian(region_flags[ex, ey])
            face_ids = (1:num_faces(rd.element_type)) .+ (e-1) * num_faces(rd.element_type)
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
                        nbr_face_ids = (1:num_faces(rd.element_type)) .+ (e_nbr-1) * num_faces(rd.element_type)
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
                                # println("match found for f = $f, e=($ex, $ey), 
                                #          enbr=($ex_nbr, $ey_nbr)")
                            end
                        end
                    end

                end # if enbr is_inside_domain
            end
        end # if e is_inside_domain
    end

    return FToF
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

function calculate_cutcells(vx, vy, curves, ds = 1e-3, arc_tol = 1e-10, corner_tol = 1e-10)

    stop_pts = find_mesh_intersections((vx, vy), curves, ds, arc_tol, corner_tol,
    closed_list=true, closure_tol=1e-12)

    # Calculate cutcells
    region_flags, cutcell_indices, cutcells = 
        define_regions((vx, vy), curves, stop_pts, binary_regions=false)

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
    function MeshData(rd, geometry, vxyz...)

Creates a cut-cell mesh where the boundary is given by `geometry`, which should be a tuple of functions. 
These functions can be generated using PathIntersections.PresetGeometries, for example:
```julia
julia> geometry = (PresetGeometries.Circle(R=0.33, x0=0, y0=0), )
```
Here, `coordinates_min`, `coordinates_max` contain `(smallest value of x, smallest value of y)` and 
`(largest value of x, largest value of y)`, and `cells_per_dimension_x/y` is the number of Cartesian grid 
cells placed along each dimension. 
"""
MeshData(rd::RefElemData, curves, cells_per_dimension;  kwargs...) = 
    MeshData(rd::RefElemData, curves, cells_per_dimension, cells_per_dimension;  kwargs...)

function MeshData(rd::RefElemData, curves, cells_per_dimension_x, cells_per_dimension_y; 
                  quad_rule_face = get_1d_quadrature(rd), 
                  coordinates_min=(-1.0, -1.0), coordinates_max=(1.0, 1.0),
                  precompute_operators=false)

    # compute intersections of curve with a background Cartesian grid.
    vx = LinRange(coordinates_min[1], coordinates_max[1], cells_per_dimension_x + 1)
    vy = LinRange(coordinates_min[2], coordinates_max[2], cells_per_dimension_y + 1)    

    # compute mesh intersections and cut cell elements.
    # `regions` is a matrix of dimensions `(cells_per_dimension_x, cells_per_dimension_y)` with 3 values:
    #   *  1: cut cell
    #   *  0: Cartesian cell
    #   * -1: excluded cells (in the cut-out region)
    region_flags, cutcells = calculate_cutcells(vx, vy, curves)

    # pack useful cut cell information together. 
    num_cartesian_cells = sum(region_flags .== 0)
    num_cut_cells = sum(region_flags .== 1)
    cut_faces_per_cell = count_cut_faces(cutcells)
    cut_face_offsets = [0; cumsum(cut_faces_per_cell)[1:end-1]] 
    cutcell_data = (; curves, region_flags, cutcells, cut_faces_per_cell, cut_face_offsets)

    # Compute volume, face points, and physical frame element scalings
    physical_frame_elements, x, y, rstxyzJ, J, xf, yf, nxJ, nyJ, Jf = 
        compute_geometric_data(rd, quad_rule_face, vx, vy, cutcell_data)
            
    # Compute face-to-face connectivity by matching face centroids
    face_centroids = compute_face_centroids(rd, xf, yf, cutcell_data)
    FToF = connect_mesh(rd, face_centroids, cutcell_data)

    # Compute node-to-node mappings
    # !!! Warning: this only works if the same quadrature rule is used for all faces! 
    num_total_faces = length(FToF)
    num_points_per_face = length(rd.rf) ÷ num_faces(rd.element_type)
    mapM = collect(reshape(1:num_points_per_face * num_total_faces, 
                           num_points_per_face, num_total_faces))
    mapP = copy(mapM)
    p = zeros(Int, num_points_per_face) # temp storage for a permutation vector
    for f in eachindex(FToF)
        idM, idP = view(mapM, :, f), view(mapM, :, FToF[f])
        xyzM = (view(xf, idM), view(yf, idM))
        xyzP = (view(xf, idP), view(yf, idP))
        StartUpDG.match_coordinate_vectors!(p, xyzM, xyzP)
        mapP[p, f] .= idP
    end
    mapB = findall(vec(mapM) .== vec(mapP)) # determine boundary nodes
        
    # compute cut-cell surface quadrature
    _, w1D = quad_rule_face
    wJf = similar(Jf)
    wJf.cartesian = reshape(Diagonal(w1D) * reshape(Jf.cartesian, length(w1D), :), size(Jf.cartesian))
    wJf.cut = reshape(Diagonal(w1D) * reshape(Jf.cut, length(w1D), :), size(Jf.cut))
    
    # The minimum number of cut cell quadrature points is `Np_cut(2 * rd.N)`. However, 
    # oversampling slightly seems to improve the conditioning of the quadrature weights.
    num_cut_quad_points = Np_cut(2 * rd.N) + 1
    xq, yq, wJq = ntuple(_ -> NamedArrayPartition(cartesian=zeros(rd.Nq, num_cartesian_cells), 
                                                  cut=zeros(num_cut_quad_points, num_cut_cells)), 3)    

    # compute quadrature rules for the Cartesian cells
    e = 1
    for ex in 1:cells_per_dimension_x, ey in 1:cells_per_dimension_y
        if is_Cartesian(region_flags[ex, ey])
            dx = vx[ex+1] - vx[ex]
            dy = vy[ey+1] - vy[ey]
            J = dx * dy / sum(rd.wq) 
            @. xq.cartesian[:, e] = dx * 0.5 * (1 + rd.rq) + vx[ex]
            @. yq.cartesian[:, e] = dy * 0.5 * (1 + rd.sq) + vy[ey]
            @. wJq.cartesian[:, e] = rd.wq * J
            e += 1
        end
    end
   
    # polynomial antidifferentiation operators for computing volume integrals
    Ix, Iy = StartUpDG.antidiff_operators(2 * rd.N)

    # refine the surface rule used to compute the volume quadrature 
    if length(quad_rule_face[1]) < 3 * rd.N + 1        
        r1D, w1D = gauss_quad(0, 0, 3 * rd.N)
    else
        r1D, w1D = quad_rule_face
    end
    # compute quadrature rules for the cut cells. integrate exactly degree 2N basis 
    for e in eachindex(cutcells)
        # compute these quantities using a higher accuracy surface quadrature rule 
        xf_element, yf_element, nxJ_element, nyJ_element, Jf_element = 
            map_points_to_cut_cell_faces(r1D, cutcells[e])                            
        nx_element = nxJ_element ./ Jf_element                
        ny_element = nyJ_element ./ Jf_element                
        wJf_element = vec(Diagonal(w1D) * reshape(Jf_element, length(w1D), :))
                    
        # compute volume integrals using numerical Green's theorem and surface integrals
        scaling = physical_frame_elements[e].scaling
        Vf = vandermonde(physical_frame_elements[e], 2 * rd.N + 1, xf_element, yf_element)
        b = 0.5 * ((Vf * Ix)' * (wJf_element .* nx_element) ./ scaling[1] + 
                   (Vf * Iy)' * (wJf_element .* ny_element) ./ scaling[2]) 
        
        # compute degree 2N basis matrix at sampled points            
        x_sampled, y_sampled = generate_sampling_points(curves, physical_frame_elements[e], 
                                                        rd, Np_cut(6 * rd.N); N_sampled = 8 * rd.N)          
        Vq = vandermonde(physical_frame_elements[e], 2 * rd.N, x_sampled, y_sampled)
        
        # naive approach to computing quadrature weights; no guarantees of positivity
        QR  = qr(Vq', ColumnNorm())
        ids = view(QR.p, 1:num_cut_quad_points)
        wq  = Vq[ids,:]' \ b
        
        quadrature_error = norm(Vq[ids,:]' * wq - b)
        quadrature_condition_number = sum(abs.(wq)) / sum(wq)
        if quadrature_condition_number > 10 || quadrature_error > 1e-13
            println("Quadrature error on element $e is $quadrature_error, " * 
                    "quadrature condition number = $quadrature_condition_number. " * 
                    "Condition number of quadrature VDM is $(cond(Vq[ids,:]')).")
        end

        view(xq.cut, :, e)  .= view(x_sampled, ids)
        view(yq.cut, :, e)  .= view(y_sampled, ids)
        view(wJq.cut, :, e) .= wq
    end

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
                       vxyz=(vx, vy), cells_per_dimension=(cells_per_dimension_x, cells_per_dimension_y), # background Cartesian grid info                       
                    )

    # get indices of cut face nodes 
    face_ids(e) = (1:(num_points_per_face * cut_faces_per_cell[e])) .+ 
                   cut_face_offsets[e] * num_points_per_face
    cut_face_node_ids = [face_ids(e) for e in 1:num_cut_cells]

    if precompute_operators == true

        # precompute cut-cell operators and store them in the `md.mesh_type.cut_cell_operators` field
        cut_face_nodes = cut_face_node_ids
        face_interpolation_matrices = Matrix{eltype(x)}[]
        lift_matrices = Matrix{eltype(x)}[]
        differentiation_matrices = Tuple{Matrix{eltype(x)}, Matrix{eltype(x)}}[]
        mass_matrices = Matrix{eltype(x)}[]
        for (e, elem) in enumerate(physical_frame_elements)

            VDM = vandermonde(elem, rd.N, x.cut[:, e], y.cut[:, e])
            Vq, Vrq, Vsq = map(A -> A / VDM, basis(elem, rd.N, xq.cut[:,e], yq.cut[:, e]))
        
            M  = Vq' * diagm(wJq.cut[:, e]) * Vq
            Qr = Vq' * diagm(wJq.cut[:, e]) * Vrq
            Qs = Vq' * diagm(wJq.cut[:, e]) * Vsq    
            Dx_e, Dy_e = M \ Qr, M \ Qs
            
            Vf = vandermonde(elem, rd.N, xf.cut[cut_face_nodes[e]], yf.cut[cut_face_nodes[e]]) / VDM

            # don't include jacobian scaling in LIFT matrix (for consistency with the Cartesian mesh)            
            _, w1D = quad_rule_face
            num_cut_faces = length(cut_face_nodes[e]) ÷ length(w1D)
            wf = vec(repeat(w1D, 1, num_cut_faces)    )

            push!(lift_matrices, M \ (Vf' * diagm(wf)))
            push!(face_interpolation_matrices, Vf)
            push!(differentiation_matrices, (Dx_e, Dy_e))
            push!(mass_matrices, M)
        end
        cut_cell_operators = (; differentiation_matrices, face_interpolation_matrices, 
                                mass_matrices, lift_matrices)

    else

        cut_cell_operators = nothing
    end
                        
    return MeshData(CutCellMesh(physical_frame_elements, cut_face_node_ids, curves, cut_cell_operators, cut_cell_data), 
                    FToF, (x, y), (xf, yf), (xq, yq), wJq, 
                    mapM, mapP, mapB, rstxyzJ, J, (nxJ, nyJ), Jf, is_periodic)

end

num_cartesian_elements(md::MeshData{DIM, <:CutCellMesh}) where {DIM} = size(md.x.cartesian, 2) 
num_cut_elements(md::MeshData{DIM, <:CutCellMesh}) where {DIM} = size(md.x.cut, 2) 

# this is used in src/MeshData.jl in `getproperty` 
function num_elements(md::MeshData{DIM, <:CutCellMesh}) where {DIM}
    # the number of elements is given by the number of columns of each component of x
    return num_cartesian_elements(md) + num_cut_elements(md)
end
