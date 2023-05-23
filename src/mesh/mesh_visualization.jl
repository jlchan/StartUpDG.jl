"""
    MeshPlotter(rd::RefElemData, md::RefElemData)

Plot recipe to plot a (possibly curved) quadrilateral or triangular mesh. Usage: `plot(MeshPlotter(...))`
"""
struct MeshPlotter{Dim, RD<:RefElemData{Dim}, MD<:MeshData{Dim}}
    rd::RD
    md::MD
end

RecipesBase.@recipe function f(rd::RefElemData, md::MeshData)
    return MeshPlotter(rd, md)
end

RecipesBase.@recipe function f(m::MeshPlotter{2})

    linecolor --> :black
    legend --> false
    aspect_ratio --> 1

    (; rd, md) = m
    (; x, y) = md
    (; Fmask) = rd
    Fmask = reshape(Fmask, length(Fmask) ÷ rd.Nfaces, rd.Nfaces)
    get_face_nodes(u, e, f) = view(u, view(Fmask, :, f), e)    

    VDM = vandermonde(Line(), rd.N, nodes(Line(), rd.N))
    Vp1D = vandermonde(Line(), rd.N, LinRange(-1, 1, 15)) / VDM
    
    xmesh, ymesh = eltype(x)[], eltype(y)[]    
    for e in 1:md.num_elements
        for f in 1:rd.Nfaces
            x_f,y_f = (x->append!(Vp1D * x, NaN)).(get_face_nodes.((x, y), e, f))
            append!(xmesh, x_f)
            append!(ymesh, y_f)
        end
    end
    return xmesh, ymesh
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

    (; VXY, EToV, fv ) = m
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
    MeshData_to_vtk(md, rd, dim, data, dataname, datatype, filename, write_data = false, equi_dist_nodes = true)

Translate the given mesh into a vtk-file.
`md` holds a `MeshData` object
`rd` holds a reference element data/`RefElemData` object. 
`data` holds an array of arrays (of size `num_nodes` by `num_elements`) with plotting data
`dataname` is an array of strings with name of the associated data
`write_data`, flag if data should be written or not (e.g., if data is not written, only the mesh will be saved as output)
`equi_dist_nodes` flag if points should be interpolated to equidstant nodes
"""
function MeshData_to_vtk(md::MeshData, rd::RefElemData{DIM}, data, dataname, filename, 
                         write_data = false, equi_dist_nodes = true) where {DIM}
                         
    # Compute the permutation between the StartUpDG order of points and vtk
    perm = SUD_to_vtk_order(rd)
    # The number of points per element
    num_lagrange_points = length(perm)
    vtk_cell_type = type_to_vtk(rd.element_type)

    # Construction of the vtkfile
    cells = [MeshCell(vtk_cell_type, perm .+ ((i-1) * num_lagrange_points)) for i in 1:md.num_elements]

    if equi_dist_nodes == true
        interpolate = vandermonde(rd.element_type, rd.N, equi_nodes(rd.element_type, rd.N)...) / rd.VDM
        coords = map(x -> vec(interpolate * x), md.xyz)
    else # don't interpolate
        coords = vec.(md.xyz) 
    end   
    
    vtkfile = vtk_grid(filename, coords..., cells)

    if write_data
        for i in 1:length(dataname)
            vtkfile[dataname[i]] = data[i]
        end
    end

    return vtk_save(vtkfile)
end


function MeshData_to_vtk(md::MeshData, rd::RefElemData, filename)
    num_tri_nodes = length(rd.approximation_type.tri.r)
    num_wedge_nodes = length(rd.r)
    node_connection = []
    wedge_add = 0
    for elems in 1:md.num_elements
        for k in 1:rd.N[1]
            add = 0
            bottom = (k-1)*num_tri_nodes
            top = k*num_tri_nodes
            for i in rd.N[2]+1:-1:1
                for j in 1:i-1
                    push!(node_connection, [j+add+bottom + wedge_add, j+1+add+bottom + wedge_add, j+i+add+bottom + wedge_add, j+add+top + wedge_add, j+1+add+top + wedge_add, j+i+add+top + wedge_add])
                    if j!=i && j!=1
                        push!(node_connection,  [j+add+bottom + wedge_add, j+i+add+bottom + wedge_add, j+i+add-1+bottom + wedge_add, j+add+top + wedge_add, j+i+add+top + wedge_add, j+i+add-1+top + wedge_add])
                    end
                end
                add = add + i
            end
        end
        wedge_add = wedge_add + num_wedge_nodes
    end
    
    total_num_elems = length(node_connection)
    vtk_cell_type = VTKCellTypes.VTK_WEDGE
    cells = [MeshCell(vtk_cell_type, node_connection[i]) for i in 1:total_num_elems]
    coords = vec.(md.xyz)
    vtkfile = vtk_grid(filename, coords..., cells)
    return vtk_save(vtkfile)
end


