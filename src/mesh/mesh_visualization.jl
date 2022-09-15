"""
    MeshPlotter(rd::RefElemData, md::RefElemData)

Plot recipe to plot a (possibly curved) quadrilateral or triangular mesh. Usage: `plot(MeshPlotter(...))`
"""
struct MeshPlotter{Dim, RD<:RefElemData{Dim}, MD<:MeshData{Dim}}
    rd::RD
    md::MD
end

RecipesBase.@recipe function f(m::MeshPlotter{2})

    linecolor --> :black
    legend --> false
    aspect_ratio --> 1

    @unpack rd, md = m
    @unpack x, y = md
    @unpack Fmask = rd
    Fmask = reshape(Fmask, length(Fmask) ÷ rd.Nfaces, rd.Nfaces)
    get_face_nodes(u, e, f) = view(u, view(Fmask, :, f), e)    

    VDM = vandermonde(Line(), rd.N, nodes(Line(), rd.N))
    Vp1D = vandermonde(Line(), rd.N, LinRange(-1, 1, 15)) / VDM
    
    xmesh, ymesh = eltype(x)[], eltype(y)[]    
    for e in 1:md.num_elements
        for f in 1:rd.Nfaces
            x_f,y_f = (x->append!(Vp1D*x, NaN)).(get_face_nodes.((x, y), e, f))
            append!(xmesh, x_f)
            append!(ymesh, y_f)
        end
    end
    return xmesh,ymesh
end

"""
    VertexMeshPlotter((VX, VY), EToV, fv)
    VertexMeshPlotter(triout::TriangulateIO)    

Plot recipe to plot a quadrilateral or triangular mesh. Usage: `plot(VertexMeshPlotter(...))`
"""
struct VertexMeshPlotter{NDIMS, Tv, Ti, Nfaces}
    VXY::NTuple{NDIMS, Vector{Tv}}
    EToV::Matrix{Ti}
    fv::NTuple{Nfaces, Vector{Int}}
end

RecipesBase.@recipe function f(m::VertexMeshPlotter)

    @unpack VXY, EToV, fv = m
    VX, VY = VXY

    linecolor --> :black
    legend --> false
    aspect_ratio --> 1
    # title --> "$(size(EToV,1)) elements"

    xmesh = eltype(VX)[]
    ymesh = eltype(VY)[]
    for vertex_ids in eachrow(EToV)
        ids = vcat(vertex_ids, vertex_ids[1])
        for f in fv            
            append!(xmesh, [VX[ids[f]]; NaN])
            append!(ymesh, [VY[ids[f]]; NaN])
        end
    end
    return xmesh, ymesh
end

function wedge_plotter(VXYZ, EToV)
    num_elem = first(size(EToV))
    num_points = first(size(VXYZ[1]))
    points = zeros(Float64, (3, num_points))
    point_data = zeros(Int32, (num_points))
    points[1, :] = VXYZ[1]
    points[2, :] = VXYZ[2]
    points[3, :] = VXYZ[3]
    for i in 1:num_points
        point_data[i] = i
    end
    cells = [MeshCell(VTKCellTypes.VTK_WEDGE, EToV[i, :]) for i in 1:num_elem]
    vtkfile = vtk_grid("mesh", points, cells)
    vtkfile["point_id", VTKPointData()] = point_data
    vtk_save(vtkfile)
end
