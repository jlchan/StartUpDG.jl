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

RecipesBase.@recipe function f(m::VertexMeshPlotter{2})

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

"""
    Meshdata_to_vtk(md, rd, dim, filename)

Translate the given mesh into a 
"""
function Meshdata_to_vtk(md::MeshData, rd::RefElemData, dim, filename)
    # Compute the permutation between the StartUpDG order of points and vtk
    perm = SUD_to_vtk_order(rd, dim)
    # The number of points per element
    num_lagrange_points = size(perm)[1]
    # Number of elements
    num_elements = size(md.EToV)[1]
    vtk_cell_type = type_to_vtk(rd.element_type)

    # Construction of the vtkfile
    cells = [MeshCell(vtk_cell_type, perm.+((i-1)*num_lagrange_points)) for i in 1:num_elements]
    # Todo: Interpolate to equidstant points 
    coords = [vec(md.xyz[i]) for i in 1:dim]
    vtkfile = []
    if dim == 1
        vtkfile = vtk_grid(filename, coords[1], cells)
    elseif dim == 2
        vtkfile = vtk_grid(filename, coords[1], coords[2], cells)
    else
        vtkfile = vtk_grid(filename, coords[1], coords[2], coords[3], cells)
    end
    vtk_save(vtkfile)
end

