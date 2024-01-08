# get face centroids for a single coordinate array
function coordinate_face_centroids(xf, md)
    Nfaces = size(md.FToF, 1)
    Nfp = size(md.xf, 1) ÷ Nfaces
    xc = reshape(xf, Nfp, :)
    return vec(typeof(xf)(sum(xc, dims=1) / size(xc, 1)))
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
function _tag_boundary_faces(boundary_face_ids, boundary_list, xyzb::NTuple{1,Tv}) where {Tv}        
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


