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
        coords = map(x -> vec(rd.Vp * x), md.xyz)
        data = [rd.Vp * data_i for data_i in data]
    else # don't interpolate
        if (rd.approximation_type isa SBP) && (rd.element_type isa Union{Tri, Tet})
            error("Support for non-interpolated SBP approximations is not supported on simplices.")
        else # if polynomial
            coords = vec.(md.xyz) 
        end
    end   
    
    vtkfile = vtk_grid(filename, coords..., cells)

    if write_data
        for (i, data_i) in enumerate(data)
            vtkfile[dataname[i]] = data_i
        end
    end

    return vtk_save(vtkfile)
end
