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
    function export_to_vtk(rd, md, data::{AbstractArray{T}}, filename; 
                           write_data = false, equi_dist_nodes = true) where {T <: Real}
    function export_to_vtk(rd, md, data::AbstractDict{String, AbstractArray{T}}, filename; 
                           write_data = false, equi_dist_nodes = true) where {T <: Real}
    function export_to_vtk(rd, md, data::AbstractArray, dataname::AbstractArray, filename; 
                           equi_dist_nodes = true)                           

Exports `data` into a vtk-file for visualization.
- `rd` is a reference element data/`RefElemData` object. 
- `md` is a `MeshData` object
- `dataname` is an array of strings with names of the associated data
- `equi_dist_nodes` flag if points should be interpolated to equidstant nodes

The argument `data` can be any of the following:
- an array of matrices of plotting data, where each matrix is size `num_nodes` by `num_elements`.
- a `Dict{String, AbstractArray{T}} where {T <: Real}`, where the keys correspond to names of each field
"""
export_to_vtk(rd, md, data::AbstractDict, filename; kwargs...) =
    export_to_vtk(rd, md, values(data), collect(keys(data)), filename; kwargs...)

# this assumes `data` is a container (e.g., vector or tuple) or matrices
export_to_vtk(rd, md, data, filename; kwargs...) = 
    export_to_vtk(rd, md, data, "Field " .* string.(eachindex(data)), filename; kwargs...)

# this is the same interface as `MeshData_to_vtk` but with `rd, md` arguments ordered differently
# for consistency and without the `write_data` kwarg
function export_to_vtk(rd, md, data, dataname, filename; equi_dist_nodes = true)
    write_data = true 
    return MeshData_to_vtk(md, rd, data, dataname, filename,
                           write_data, equi_dist_nodes)
end

"""
    MeshData_to_vtk(md, rd, data, dataname, filename, write_data = false, equi_dist_nodes = true)

Translate the given mesh into a vtk-file.
`md` holds a `MeshData` object
`rd` holds a reference element data/`RefElemData` object. 
`data` holds an array of matrices (of size `num_nodes` by `num_elements`) with plotting data
`dataname` is an array of strings with name of the associated data
`write_data`, flag if data should be written or not (e.g., if data is not written, only the mesh will be saved as output)
`equi_dist_nodes` flag if points should be interpolated to equidstant nodes
"""
function MeshData_to_vtk(md::MeshData, rd::RefElemData, data, dataname, filename, 
                         write_data = false, equi_dist_nodes = true) 
                         
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
        for (i, data_i) in enumerate(data)
            vtkfile[dataname[i]] = data_i
        end
    end

    return vtk_save(vtkfile)
end

"""
MeshData_to_vtk(md, rd, data, dataname, filename, write_data = false, equi_dist_nodes = true)

Translate the given mesh into a vtk-file.
`md` holds a `MeshData` object
`rd` holds a reference element data/`RefElemData` of a TensorProductWedge
`data` holds an array of matrices (of size `num_nodes` by `num_elements`) with plotting data
`dataname` is an array of strings with name of the associated data
`write_data`, flag if data should be written or not (e.g., if data is not written, only the mesh will be saved as output)
`equi_dist_nodes` flag if points should be interpolated to equidstant nodes
"""
function MeshData_to_vtk(md::MeshData, rd::RefElemData{3, <:Wedge, <:TensorProductWedge}, data, dataname, filename, 
                    write_data = false, equi_dist_nodes = true)
    # Number of nodes in the triangular base
    num_tri_nodes = length(rd.approximation_type.tri.r)
    # Number of all wedges in the mesh
    num_wedge_nodes = length(rd.r)
    # Array to fill with the connectivity data for vtk
    node_connection = Vector{Int}[]
    # Shift for the number of nodes of the wedges we already computed
    wedge_add = 0
    for elems in 1:md.num_elements
        # Iterate of the number of `levels` in the TensorProductWedge. Each level produces a new set of linear 'VTK_WEDGE'
        # The node-ids increase columnwise. 
        for k in 1:rd.N[1]
            # A shift for the colums already passed. 
            add = 0
            # A shift to the buttom of the current wedge. 
            bottom = (k-1)*num_tri_nodes
            # A shift to the top of the current wedge. 
            top = k*num_tri_nodes
            # i and j iterate over the triangular base of the wedge. 
            for i in rd.N[2]+1:-1:1
                for j in 1:i-1
                    # 3 buttom and 3 top nodes form a wedge. 
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

    # The total number of written vtk-elements is higher than the number of md.num_elements
    total_num_elems = length(node_connection)
    vtk_cell_type = VTKCellTypes.VTK_WEDGE
    # Fill the cells-Array for VTK
    cells = [MeshCell(vtk_cell_type, node_connection[i]) for i in 1:total_num_elems]
    # Coordinates for VTK
    

    if equi_dist_nodes
        # Construct an interpolation matrix for the triangular basis. 
        tri_interpolate = vandermonde(rd.approximation_type.tri.element_type, rd.approximation_type.tri.N, equi_nodes(rd.approximation_type.tri.element_type, rd.approximation_type.tri.N)...)/rd.approximation_type.tri.VDM
        # Construct an interpolation matrix for the linear basis
        line_interpolate = vandermonde(rd.approximation_type.line.element_type, rd.approximation_type.line.N, collect(LinRange(-1, 1, rd.approximation_type.line.N+1))) / rd.approximation_type.line.VDM
        
        # storage for the equi-distant nodes
        coords = (similar(md.x), similar(md.y), similar(md.z))

        # Get the number of points per element 
        tri_num_points = length(rd.approximation_type.tri.r)
        line_num_points = length(rd.approximation_type.line.r)

        # equi-distant nodes for triangular basis (hence only consider x-y-coords. )
        for dim in 1:2
            # iterate over all elements
            for elem in 1:md.num_elements
                for i in 1:line_num_points
                    # interpolate each slice of the wedge
                    range = ((i-1)*tri_num_points + 1):(i*tri_num_points)
                    coords[dim][range, elem] = tri_interpolate * md.xyz[dim][range, elem]
                end
            end
        end  

        # equi-distant nodes for linear basis (hence only consider z-coords)
        for elem in 1:md.num_elements
            # Get the z-coord of each wedge-slice
            z_coords = [md.z[(i-1)*tri_num_points + 1, elem] for i in 1:line_num_points]
            # interpolate
            z_tmp = line_interpolate * z_coords
            # each slice has 'tri_num_points' nodes. repeat the z-value 'tri_num_points' times and add to z-coords. 
            for i in 1:line_num_points
                range = ((i-1)*tri_num_points + 1):(i*tri_num_points)
                coords[3][range, elem] = fill(z_tmp[i], tri_num_points)
            end
        end
        coords = vec.(coords)  
    else
        coords = vec.(md.xyz)
    end
    vtkfile = vtk_grid(filename, coords..., cells)
    
    if write_data
        for (i, data_i) in enumerate(data)
            vtkfile[dataname[i]] = data_i
        end
    end
    return vtk_save(vtkfile)
end


