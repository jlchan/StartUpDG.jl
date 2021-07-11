
"""
    function Triangulate.triangulate(triin::TriangulateIO, maxarea, minangle=20)

Convenience routine to avoid writing out `@sprintf` each time. Returns a `TriangulateIO` object.
"""
function Triangulate.triangulate(triin::TriangulateIO, maxarea, minangle=20)
    angle = @sprintf("%.15f", minangle)
    area  = @sprintf("%.15f", maxarea)
    triout, _ = triangulate("pa$(area)q$(angle)Q", triin)
    return triout
end

"""
    function triangulateIO_to_VXYEToV(triout::TriangulateIO)

Computes `VX`,`VY`,`EToV` from a `TriangulateIO` object.
"""
function triangulateIO_to_VXYEToV(triout::TriangulateIO)
    VX,VY = (triout.pointlist[i,:] for i = 1:size(triout.pointlist,1))
    EToV = permutedims(triout.trianglelist)
    Base.swapcols!(EToV,2,3) # to match MeshData ordering
    return VX,VY,Matrix{Int}(EToV)
end

"""
    function get_boundary_face_labels(triout::TriangulateIO, md::MeshData{2})

Find Triangle segment labels of boundary faces. Returns two arguments:
- `boundary_face_tags`: tags of faces on the boundary
- `boundary_faces`: list of faces on the boundary of the domain
"""
function get_boundary_face_labels(triout::TriangulateIO, rd::RefElemData{2, Tri}, md::MeshData{2})
    segmentlist = sort(triout.segmentlist,dims=1)
    boundary_faces = findall(vec(md.FToF) .== 1:length(md.FToF))
    boundary_face_tags = zeros(Int,length(boundary_faces))
    for (f,boundary_face) in enumerate(boundary_faces)
        element = (boundary_face - 1) ÷ rd.Nfaces + 1
        face    = (boundary_face - 1) % rd.Nfaces + 1
        vertices_on_face = sort(md.EToV[element,rd.fv[face]])
        tag_id = findfirst(c->view(segmentlist,:,c)==vertices_on_face,axes(segmentlist,2))
        boundary_face_tags[f] = triout.segmentmarkerlist[tag_id]
    end
    return boundary_face_tags, boundary_faces
end

"""
    function get_node_boundary_tags(triout::TriangulateIO,md::MeshData{2},rd::RefElemData{2,Tri})

Computes `node_tags` = `Nfp` x `Nfaces * num_elements` array where each entry is a Triangulate.jl tag number.
"""
function get_node_boundary_tags(triout::TriangulateIO, rd::RefElemData{2, Tri}, md::MeshData{2})
    boundary_face_tags, boundary_faces = get_boundary_face_labels(triout, rd, md)
    node_tags = zeros(Int, size(md.xf, 1) ÷ rd.Nfaces, md.K * rd.Nfaces) # make Nfp x Nfaces*num_elements
    for (i, boundary_face) in enumerate(boundary_faces)
        node_tags[:, boundary_face] .= boundary_face_tags[i]
    end
    node_tags = reshape(node_tags, size(md.xf)...)
end


function tag_boundary_faces(triout::TriangulateIO,
                            rd::RefElemData{2, Tri}, md::MeshData{2},
                            boundary_list::Dict{Symbol, Int})

    boundary_face_ids_list = _tag_boundary_faces(triout,rd,md,boundary_list)
    return Dict(Pair.(keys(boundary_list), boundary_face_ids_list))
end

"""
    function tag_boundary_faces(triout::TriangulateIO,
                                rd::RefElemData{2,Tri}, md::MeshData{2},
                                boundary_list::Union{NamedTuple,Dict{Symbol,Int}})

Here, `boundary_list` is a `Dict` (or `NamedTuple`) whose values are the boundary tags for a
`TriangulateIO` mesh format. The output is a `Dict` or `NamedTuple` with keys given by
`boundary_list` and `values` equal to vectors of faces on that given boundary.

Example usage:
```julia
julia> using Triangulate, StartUpDG
julia> triout = scramjet()
julia> rd = RefElemData(Tri(),N=1)
julia> md = MeshData(triangulateIO_to_VXYEToV(triout)...,rd)
julia> tag_boundary_faces(triout,rd,md, Dict(:wall=>1, :inflow=>2, :outflow=>3))
```
"""
function tag_boundary_faces(triout::TriangulateIO,
                            rd::RefElemData{2, Tri}, md::MeshData{2},
                            boundary_list::NamedTuple)

    boundary_face_ids_list = _tag_boundary_faces(triout,rd,md,boundary_list)
    return NamedTuple(Pair.(keys(boundary_list), boundary_face_ids_list))
end

# this version works for both boundary_list::Union{NamedTuple, Dict{Symbol,Int}}
function _tag_boundary_faces(triout::TriangulateIO,
                            rd::RefElemData{2, Tri}, md::MeshData{2},
                            boundary_list)
    boundary_face_tags, boundary_faces = get_boundary_face_labels(triout, rd, md)
    boundary_face_ids_list = Vector{Int}[]
    for boundary_face_flag in values(boundary_list)
        push!(boundary_face_ids_list, boundary_faces[findall(@. boundary_face_tags == boundary_face_flag)])
    end
    return boundary_face_ids_list
end

"""
    function refine(triout, h, href = h/2)

Refinement of a previous mesh given the current mesh size `h`. Preserves boundary/volume tags.
"""
function refine(triout, h, href = h/2)
    angle = @sprintf("%.15f",20)
    area  = @sprintf("%.15f",href^2)
    triout2,_ = triangulate("rpa$(area)q$(angle)Q", triout)
    return triout2
end

VertexMeshPlotter(triout::TriangulateIO) = VertexMeshPlotter(triangulateIO_to_VXYEToV(triout)..., face_vertices(Tri()))

"""
    BoundaryTagPlotter(triout::TriangulateIO)

Plot recipe to visualize boundary tags by color. Usage: `plot(BoundaryTagPlotter(triout))`
"""
struct BoundaryTagPlotter
    triout::TriangulateIO
end

RecipesBase.@recipe function f(m::BoundaryTagPlotter)
    triout = m.triout
    tags = unique(triout.segmentmarkerlist)
    num_colors = length(tags)
    colors = distinguishable_colors(num_colors)
    xseg = zeros(2, size(triout.segmentlist,2))
    yseg = zeros(2, size(triout.segmentlist,2))
    segcolor = HSV{Float32}[]
    for (col,segment) in enumerate(eachcol(triout.segmentlist))
        xseg[:,col] .= triout.pointlist[1,segment]
        yseg[:,col] .= triout.pointlist[2,segment]
        push!(segcolor, colors[triout.segmentmarkerlist[col]])
    end
    for i = 1:num_colors
        color_ids = findall(triout.segmentmarkerlist .== tags[i])

        # NaN separators for distinct lines
        x_i = vec([xseg[:, color_ids]; fill(NaN, length(color_ids))'])
        y_i = vec([yseg[:, color_ids]; fill(NaN, length(color_ids))'])

        @series begin
            marker --> :circle
            seriescolor --> permutedims(segcolor[color_ids]),
            ratio --> 1
            label --> string(tags[i])
            x_i,y_i
        end
    end
    # we already added all shapes in @series so we don't want to return a series
    # here (we are returning an empty series which is not added to the legend.)
    primary := false
    ()
end
