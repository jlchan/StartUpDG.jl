# get face centroids for a single coordinate array
function coordinate_face_centroids(xf, md)
    Nfaces = size(md.FToF, 1)
    Nfp = size(md.xf, 1) ÷ Nfaces
    xc = reshape(xf, Nfp, :)
    return vec(typeof(xf)(sum(xc, dims = 1) / size(xc, 1)))
end

"""
    function boundary_face_centroids(md)

Returns face centroids and `boundary_face_ids` on the boundaries of the domain given by md::MeshData.
"""
function boundary_face_centroids(md)
    compute_face_centroids(md) = map(x -> coordinate_face_centroids(x, md), md.xyzf)
    xyzc = compute_face_centroids(md)
    boundary_face_ids = findall(vec(md.FToF) .== 1:length(md.FToF))

    # compute coordinates of face centroids on the boundary
    xyzb = map(x -> x[boundary_face_ids], xyzc)
    return xyzb, boundary_face_ids
end

"""
    function tag_boundary_faces(md, boundary_name::Symbol = :entire_boundary)
    function tag_boundary_faces(md, boundary_list::Dict{Symbol, <:Function})

When called without arguments, just returns `Dict(:entire_boundary => boundary_faces)``.
    
Example usage: 
```julia
julia> rd = RefElemData(Tri(), N=1)
julia> md = MeshData(uniform_mesh(Tri(), 2)..., rd)
julia> on_bottom_boundary(x, y, tol = 1e-13) = abs(y+1) < tol
julia> on_top_boundary(x, y, tol = 1e-13) = abs(y-1) < tol
julia> tag_boundary_faces(Dict(:bottom => on_bottom_boundary,
                               :top    => on_top_boundary), md)
```
"""
tag_boundary_faces(md, ::Nothing) = tag_boundary_faces(md)

function tag_boundary_faces(md, boundary_name::Symbol = :entire_boundary)
    return Dict(boundary_name => findall(vec(md.FToF) .== 1:length(md.FToF)))
end

function _tag_boundary_faces(boundary_face_ids, boundary_list, xyzb)
    boundary_face_ids_list = Vector{Int}[]
    for boundary_face_flag in values(boundary_list)
        push!(boundary_face_ids_list, boundary_face_ids[boundary_face_flag.(zip(xyzb...))])
    end
    return boundary_face_ids_list
end

# specialization to 1D
function _tag_boundary_faces(boundary_face_ids, boundary_list,
                             xyzb::NTuple{1, Tv}) where {Tv}
    boundary_face_ids_list = Vector{Int}[]
    for boundary_face_flag in values(boundary_list)
        push!(boundary_face_ids_list, boundary_face_ids[boundary_face_flag.(xyzb[1])])
    end
    return boundary_face_ids_list
end

# todo: should I make this version with NamedTuples the default?
function tag_boundary_faces(md, boundary_list::Dict{Symbol, <:Function})
    xyzb, boundary_face_ids = boundary_face_centroids(md)
    boundary_face_ids_list = _tag_boundary_faces(boundary_face_ids, boundary_list, xyzb)
    return Dict(Pair.(keys(boundary_list), boundary_face_ids_list))
end

function tag_boundary_faces(md, boundary_list::NamedTuple)
    xyzb, boundary_face_ids = boundary_face_centroids(md)
    boundary_face_ids_list = _tag_boundary_faces(boundary_face_ids, boundary_list, xyzb)
    return NamedTuple(Pair.(keys(boundary_list), boundary_face_ids_list))
end

function tag_boundary_nodes(rd, md, boundary_list::NamedTuple)
    boundary_faces = tag_boundary_faces(md, boundary_list)
    mapM = reshape(md.mapM, size(md.mapM, 1) ÷ rd.Nfaces, rd.num_faces * md.num_elements)
    node_tags = (mapM[:, getproperty(boundary_faces, tag)] for tag in keys(boundary_faces))
    return NamedTuple(Pair.(keys(boundary_list), node_tags))
end

function tag_boundary_nodes(rd, md, boundary_list::Dict)
    boundary_faces = tag_boundary_faces(md, boundary_list)
    mapM = reshape(md.mapM, size(md.mapM, 1) ÷ rd.Nfaces, rd.num_faces * md.num_elements)
    node_tags = (mapM[:, boundary_faces[tag]] for tag in keys(boundary_faces))
    return Dict(Pair.(keys(boundary_list), node_tags))
end

"""
    function tag_boundary_faces(md::MeshData{2}, edges_dict,
                                boundary_names = :all; atol = 1e-13)
 
Map edges from Gmsh `edges_dict` to boundary edge/face indices in `md`.
This follows the approach: Compare the centroid of each DG boundary edge/face
to the precomputed (based on mesh info) midpoints of the edges/faces with the same physical tag.
 
Arguments:
- `md`: MeshData object
- `edges_dict`: Dictionary from `build_edges_dict()` with `:midpoints` field
- `boundary_names`: Symbol(s) to extract (default `:all` for all boundaries)
- `atol`: Absolute tolerance for coordinate matching (default 1e-13)
 
Returns `Dict{Symbol, Vector{Int}}` mapping boundary names to face indices in `md`
 
Example usage:
```julia
coords, EToV, edges_dict, elem_type = read_Gmsh_2D_v2("mesh.msh")
md = MeshData(coords, EToV, rd)
boundary_faces = tag_boundary_faces(md, edges_dict)
# boundary_faces[:bottom] => [1, 2, 5, 6, ...] (face indices)
```
"""
function tag_boundary_faces(md::MeshData{2}, edges_dict, 
                            boundary_names::Union{Symbol, Vector{Symbol}} = :all;
                            atol = 1e-13)
    
    # Compute boundary face centroids (using existing StartUpDG function)
    xyzb, boundary_face_ids = boundary_face_centroids(md)
    xb, yb = xyzb[1], xyzb[2] # Face centroids
    
    # Build a list of all Gmsh edge midpoints with their tags and names
    gmsh_edge_midpoints = [] # Vector of (x, y, tag, name)
    
    for (name, info) in edges_dict
        tag = info[:tag]
        midpoints = get(info, :midpoints, [])
        for (mx, my) in midpoints
            push!(gmsh_edge_midpoints, (mx, my, tag, name))
        end
    end
    
    # Determine which boundaries to include
    if boundary_names isa Symbol && boundary_names == :all
        names_to_use = keys(edges_dict) # Select all detected boundaries
    else
        names_to_use = boundary_names isa Symbol ? [boundary_names] : boundary_names
    end

    edges_per_symbol = Dict{Symbol, Vector{Int}}()
    
    for name in names_to_use
        if !haskey(edges_dict, name)
            @warn "Boundary '$name' not found in edges_dict. Skipping."
            continue
        end
        
        matched_faces = Int[] # faces per boundary name
        target_tag = edges_dict[name][:tag] # mesh info
        
        # For each Gmsh edge with this boundary tag
        for (mx, my, tag, gmsh_name) in gmsh_edge_midpoints
            if tag != target_tag # Skip edges that don't match the current boundary tag
                continue
            end

            # Find first (up to tolerance) matching DG boundary face
            for (face_idx, face_id) in enumerate(boundary_face_ids)
                face_x = xb[face_idx]
                face_y = yb[face_idx]
                dist = sqrt((face_x - mx)^2 + (face_y - my)^2)
                if dist <= atol
                    push!(matched_faces, face_id)
                    break
                end
            end
        end
        
        edges_per_symbol[name] = sort!(matched_faces) # sorting is not strictly necessary
    end
    
    return edges_per_symbol
end
