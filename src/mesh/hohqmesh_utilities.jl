struct ISM end
struct ISM_V2 end

struct HOHQMeshData{NDIMS, TV, TE, TP, TC, TB}
    VXYZ::NTuple{NDIMS, TV}
    EToV::TE
    polydeg::TP
    curved_elements::TC
    boundary_tags::TB
end

using NodesAndModes: face_basis

# why does this one seem to work?
get_HOHQMesh_to_StartUpDG_face_ordering(::Quad) = SVector(2, 4, 1, 3) 
get_HOHQMesh_to_StartUpDG_face_ordering(::Tri) = SVector(1, 2, 3) 

vertex_reordering(::Quad) = SVector(1, 2, 4, 3)
vertex_reordering(::Tri) = SVector(1, 3, 2)

function get_HOHQMesh_ids(::Quad, curved_edges, f, polydeg) 
    active_face_offset = max.(0, cumsum(curved_edges) .- 1) * (polydeg + 1)
    return (1:polydeg+1) .+ active_face_offset[f]
end

function get_HOHQMesh_ids(::Tri, curved_edges, f, polydeg)
    active_face_offset = max.(0, cumsum(curved_edges) .- 1) * (polydeg + 1)
    if f == 1 || f == 2
        return (1:polydeg+1) .+ active_face_offset[f]
    else # reverse order for face 3
        return (polydeg+1:-1:1) .+ active_face_offset[f]
    end    
end

function MeshData(hmd::HOHQMeshData{2}, rd::RefElemData)
    (; VXYZ, EToV) = hmd    
    md = MeshData(VXYZ, EToV[:, vertex_reordering(rd.element_type)], rd)

    # interpolation matrix from chebyshev_to_lobatto nodes
    r_chebyshev = [-cos(j * pi / hmd.polydeg) for j in 0:hmd.polydeg]
    chebyshev_to_lobatto = 
        vandermonde(Line(), hmd.polydeg, nodes(Line(), rd.N)) / vandermonde(Line(), hmd.polydeg, r_chebyshev)

    warp_face_nodes_to_volume_nodes = 
        face_basis(rd.element_type, rd.N, rd.rst...) / face_basis(rd.element_type, rd.N, rd.r[rd.Fmask], rd.s[rd.Fmask])           

    HOHQMesh_to_StartUpDG_face_ordering = get_HOHQMesh_to_StartUpDG_face_ordering(rd.element_type)
    curved_face_coordinates = ntuple(_ -> similar(md.x[rd.Fmask, 1]), 2)
    fids = reshape(1:length(rd.Fmask), :, rd.num_faces)
    x, y = copy.(md.xyz)
    for curved_elem in hmd.curved_elements
        (; element, curved_edges, curved_edge_coordinates) = curved_elem

        # initialize face coordinates as linear coordinates
        curved_face_coordinates[1] .= md.x[rd.Fmask, element]
        curved_face_coordinates[2] .= md.y[rd.Fmask, element]

        ids = 1:rd.N+1
        for f in 1:rd.num_faces
            f_HOHQMesh = HOHQMesh_to_StartUpDG_face_ordering[f]
            if curved_edges[f_HOHQMesh] == 1     
                # if isdefined(Main, :Infiltrator)
                #     Main.Infiltrator.infiltrate(@__MODULE__, Base.@locals, @__FILE__, @__LINE__)
                # end
                ids_HOHQMesh = get_HOHQMesh_ids(rd.element_type, curved_edges, f, hmd.polydeg)
                curved_lobatto_coordinates = chebyshev_to_lobatto * curved_edge_coordinates[ids_HOHQMesh, :]
                curved_face_coordinates[1][fids[:, f]] .= curved_lobatto_coordinates[:, 1]
                curved_face_coordinates[2][fids[:, f]] .= curved_lobatto_coordinates[:, 2]

                # move onto next set of nodes if there is one
                ids = ids .+ (rd.N+1) 
            end
        end
        
        x[:, element] .= warp_face_nodes_to_volume_nodes * curved_face_coordinates[1]
        y[:, element] .= warp_face_nodes_to_volume_nodes * curved_face_coordinates[2]
    end

    md_curved = MeshData(rd, md, x, y) 

    return md_curved
end

function read_HOHQMesh(filename::String)
    f = open(filename)
    lines = readlines(f)

    if contains(lines[1], "ISM-V2")
        mesh_format = ISM_V2()
        deleteat!(lines, 1)
    else
        mesh_format = ISM()
    end

    return HOHQMeshData(_read_HOHQMesh(lines, mesh_format)...)
end

# container for curved data for a HOHQMesh element
struct CurvedHOHQMeshElement
    element::Int
    curved_edges::Vector{Int}
    curved_edge_coordinates::Matrix{Float64}
end

function _read_HOHQMesh(lines, mesh_format)
    if mesh_format isa ISM_V2
        num_nodes, num_edges, nelements, polydeg = parse.(Int, split(popat!(lines, 1)))
    else
        num_nodes, nelements, polydeg = parse.(Int, split(popat!(lines, 1)))
    end
    VX, VY, VZ = ntuple(_ -> zeros(num_nodes), 3)
    for i in 1:num_nodes
        VX[i], VY[i], VZ[i] = parse.(Float64, split(lines[i]))
    end
    deleteat!(lines, 1:num_nodes)

    if mesh_format isa ISM_V2
        # ignore edges for now 
        deleteat!(lines, 1:num_edges)
    end

    curved_elements = CurvedHOHQMeshElement[]
    nvertices = length(split(lines[1]))    
    num_faces = nvertices == 4 ? 4 : 6 # in 2D, 4 faces. In 3D, 6 faces
    EToV = zeros(Int, nelements, nvertices)
    boundary_tags = Matrix{String}(undef, nelements, num_faces)
    for e in 1:nelements
        EToV[e, :] .= parse.(Int, split(lines[1]))
        curved_edges = parse.(Bool, split(lines[2]))

        deleteat!(lines, 1:2) # move onto next lines

        if all(curved_edges .== false)
            # do nothing
        else
            num_curved_faces = count(curved_edges)
            num_face_nodes = nvertices == 4 ? (polydeg + 1) : (polydeg + 1)^2
            curved_nodes = 1:num_face_nodes * num_curved_faces
            curved_edge_coordinates = mapreduce(vcat, lines[curved_nodes]) do s
                (parse.(Float64, split(s)))'
            end
            push!(curved_elements, 
                  CurvedHOHQMeshElement(e, curved_edges, curved_edge_coordinates))
            deleteat!(lines, curved_nodes) # move on to next lines
        end

        tags = split(popat!(lines, 1))
        for i in eachindex(tags)
            boundary_tags[e, i] = tags[i]
        end        
    end

    if nvertices == 4
        return (VX, VY), EToV, polydeg, curved_elements, boundary_tags
    else
        return (VX, VY, VZ), EToV, polydeg, curved_elements, boundary_tags
    end
end

# reads a curved triangular .mesh file (not a public HOHQMesh utility yet)
function read_HOHQMesh(filename::String, element_type::Union{Tri, Tet})

    f = open(filename)
    lines = readlines(f)

    num_nodes, nelements, polydeg = parse.(Int, split(popat!(lines, 1)))
    VX, VY, VZ = ntuple(_ -> zeros(num_nodes), 3)
    for i in 1:num_nodes
        VX[i], VY[i], VZ[i] = parse.(Float64, split(lines[i]))
    end
    deleteat!(lines, 1:num_nodes)

    curved_elements = CurvedHOHQMeshElement[]
    nvertices = length(split(lines[1]))
    
    # for Tri/Tet meshes, we still have 4 faces for Tri and 6 faces for Tet. 
    # These correspond to collapsed quad/hex faces, respectively.
    num_faces = nvertices == 3 ? 4 : 6 

    EToV = zeros(Int, nelements, nvertices)
    boundary_tags = Matrix{String}(undef, nelements, num_faces)
    for e in 1:nelements
        EToV[e, :] .= parse.(Int, split(lines[1]))
        curved_edges = parse.(Bool, split(lines[2]))

        deleteat!(lines, 1:2) # move onto next lines

        if all(curved_edges .== false)
            # do nothing
        else
            num_curved_faces = count(curved_edges)            

            # note: HOHQMesh.jl uses (p+1)^2 nodes per face of the triangle, which 
            # correspond to collapsed coordinates on the quad face. 
            num_face_nodes = element_type isa Tri ? (polydeg + 1) : (polydeg + 1)^2 

            curved_nodes = 1:num_face_nodes * num_curved_faces
            curved_edge_coordinates = mapreduce(vcat, lines[curved_nodes]) do s
                (parse.(Float64, split(s)))'
            end
            push!(curved_elements, 
                  CurvedHOHQMeshElement(e, curved_edges, curved_edge_coordinates))
            deleteat!(lines, curved_nodes) # move on to next lines
        end

        tags = split(popat!(lines, 1))
        for i in eachindex(tags)
            boundary_tags[e, i] = tags[i]
        end
    end
    
    if element_type isa Tri
        return HOHQMeshData((VX, VY), EToV, polydeg, curved_elements, boundary_tags)
    elseif element_type isa Tet
        return HOHQMeshData((VX, VY, VZ), EToV, polydeg, curved_elements, boundary_tags)
    end
end        