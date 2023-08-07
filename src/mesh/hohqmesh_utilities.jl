struct HOHQMeshData{NDIMS, TV, TE, TP, TC, TB}
    VXYZ::NTuple{NDIMS, TV}
    EToV::TE
    polydeg::TP
    curved_elements::TC
    boundary_tags::TB
end

# HOHQMesh mesh formats
struct ISM end
struct ISM_V2 end

# for transfinite interpolation
using NodesAndModes: face_basis

vertex_reordering(::Quad) = SVector(1, 2, 4, 3)
vertex_reordering(::Tri) = SVector(1, 3, 2)

# see https://trixi-framework.github.io/HOHQMesh/TheISMMeshFileFormats 
# for more details on HOHQMesh's face ordering. 
get_HOHQMesh_to_StartUpDG_face_ordering(::Quad) = SVector(3, 2, 4, 1) 
get_HOHQMesh_to_StartUpDG_face_ordering(::Tri) = SVector(1, 2, 3) 
get_HOHQMesh_to_StartUpDG_face_ordering(::Tet) = 1:4

function get_HOHQMesh_ids(::Quad, curved_faces, f_HOHQMesh, f, polydeg) 
    active_face_offset = max.(0, cumsum(curved_faces) .- 1) * (polydeg + 1)
    return (1:polydeg+1) .+ active_face_offset[f_HOHQMesh]
end

function get_HOHQMesh_ids(::Tri, curved_faces, f_HOHQMesh, f, polydeg)
    active_face_offset = max.(0, cumsum(curved_faces) .- 1) * (polydeg + 1)
    if f == 1 || f == 2
        return (1:polydeg+1) .+ active_face_offset[f_HOHQMesh]
    else # reverse node ordering for face 3
        return (polydeg+1:-1:1) .+ active_face_offset[f_HOHQMesh]
    end    
end

# type for `md.mesh_type`
struct HOHQMeshType{HMD <: HOHQMeshData, TB}
    hmd::HMD 
    boundary_faces::TB
end

# for quads and tris
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
    curved_face_coordinates = ntuple(_ -> similar(reshape(md.x[rd.Fmask, 1], :, rd.num_faces)), 2)
    x, y = copy.(md.xyz)
    for curved_elem in hmd.curved_elements
        (; element, curved_faces, curved_edge_coordinates) = curved_elem

        # initialize face coordinates as linear coordinates
        curved_face_coordinates[1] .= reshape(md.x[rd.Fmask, element], :, rd.num_faces)
        curved_face_coordinates[2] .= reshape(md.y[rd.Fmask, element], :, rd.num_faces)

        for (f_HOHQMesh, f) in enumerate(HOHQMesh_to_StartUpDG_face_ordering)
            if curved_faces[f_HOHQMesh] == 1                     
                ids_HOHQMesh = get_HOHQMesh_ids(rd.element_type, curved_faces, f_HOHQMesh, f, hmd.polydeg)

                curved_lobatto_coordinates = chebyshev_to_lobatto * curved_edge_coordinates[ids_HOHQMesh, :]
                curved_face_coordinates[1][:, f] .= curved_lobatto_coordinates[:, 1]
                curved_face_coordinates[2][:, f] .= curved_lobatto_coordinates[:, 2]
            end
        end
       
        x[:, element] .= warp_face_nodes_to_volume_nodes * vec(curved_face_coordinates[1])
        y[:, element] .= warp_face_nodes_to_volume_nodes * vec(curved_face_coordinates[2])
    end

    md_curved = MeshData(rd, md, x, y)

    # find boundary tags
    boundary_strings = unique(hmd.boundary_tags)
    deleteat!(boundary_strings, findfirst(boundary_strings .=== "---")) # remove non-boundary faces
    face_indices = [map(x -> global_face_index(x, rd.num_faces, HOHQMesh_to_StartUpDG_face_ordering), 
                        findall(hmd.boundary_tags .== name)) for name in boundary_strings]
    boundary_faces = NamedTuple(Pair.(Symbol.(boundary_strings), face_indices))

    return @set md_curved.mesh_type = HOHQMeshType(hmd, boundary_faces)
end

function global_face_index(local_elem_and_face::CartesianIndex, num_faces_per_element, face_ordering) 
    element, face = local_elem_and_face.I
    return face_ordering[face] + (element - 1) * num_faces_per_element
end

# see https://trixi-framework.github.io/HOHQMesh/TheISMMeshFileFormats 
# for more details on HOHQMesh's face ordering. 
get_HOHQMesh_to_StartUpDG_face_ordering(::Hex) = SVector(3, 4, 5, 2, 6, 1)
vertex_reordering(::Hex) = SVector(1, 4, 2, 3, 5, 8, 6, 7)

function get_HOHQMesh_ids(::Hex, curved_faces, f_HOHQMesh, f, polydeg) 
    num_face_nodes = (polydeg + 1)^2
    active_face_offset = max.(0, cumsum(curved_faces) .- 1) * num_face_nodes
    if f == 5|| f == 6
        # permute node ordering for these faces
        ids = vec(permutedims(reshape(1:num_face_nodes, polydeg+1, polydeg+1)))
        return ids .+ active_face_offset[f_HOHQMesh]
    else
        return (1:num_face_nodes) .+ active_face_offset[f_HOHQMesh]
    end
end

# for hexahedra
function MeshData(hmd::HOHQMeshData{3}, rd::RefElemData{3, <:Hex})
    (; VXYZ, EToV) = hmd    
    md = MeshData(VXYZ, EToV[:, vertex_reordering(rd.element_type)], rd)
    
    # interpolation matrix from chebyshev_to_lobatto nodes
    r_chebyshev = vec([-cos(j * pi / hmd.polydeg) for j in 0:hmd.polydeg, i in 0:hmd.polydeg])
    s_chebyshev = vec([-cos(j * pi / hmd.polydeg) for i in 0:hmd.polydeg, j in 0:hmd.polydeg])

    chebyshev_to_lobatto = 
        vandermonde(face_type(rd.element_type), hmd.polydeg, nodes(face_type(rd.element_type), rd.N)...) / vandermonde(face_type(rd.element_type), hmd.polydeg, r_chebyshev, s_chebyshev)

    warp_face_nodes_to_volume_nodes = 
        face_basis(rd.element_type, rd.N, rd.rst...) / face_basis(rd.element_type, rd.N, rd.r[rd.Fmask], rd.s[rd.Fmask], rd.t[rd.Fmask])

    HOHQMesh_to_StartUpDG_face_ordering = get_HOHQMesh_to_StartUpDG_face_ordering(rd.element_type)
    curved_face_coordinates = ntuple(_ -> similar(reshape(md.x[rd.Fmask, 1], :, rd.num_faces)), 3)
    x, y, z = copy.(md.xyz)
    for curved_elem in hmd.curved_elements
        (; element, curved_faces, curved_edge_coordinates) = curved_elem

        # initialize face coordinates as linear coordinates
        curved_face_coordinates[1] .= reshape(md.x[rd.Fmask, element], :, rd.num_faces)
        curved_face_coordinates[2] .= reshape(md.y[rd.Fmask, element], :, rd.num_faces)
        curved_face_coordinates[3] .= reshape(md.z[rd.Fmask, element], :, rd.num_faces)

        for (f_HOHQMesh, f) in enumerate(HOHQMesh_to_StartUpDG_face_ordering)            
            if curved_faces[f_HOHQMesh] == 1
                # incorporate any permutations or face reorderings into indices
                ids_HOHQMesh = get_HOHQMesh_ids(rd.element_type, curved_faces, f_HOHQMesh, f, hmd.polydeg)

                curved_lobatto_coordinates = chebyshev_to_lobatto * curved_edge_coordinates[ids_HOHQMesh, :]
                curved_face_coordinates[1][:, f] .= curved_lobatto_coordinates[:, 1]
                curved_face_coordinates[2][:, f] .= curved_lobatto_coordinates[:, 2]
                curved_face_coordinates[3][:, f] .= curved_lobatto_coordinates[:, 3]
            end
        end        

        x[:, element] .= warp_face_nodes_to_volume_nodes * vec(curved_face_coordinates[1])
        y[:, element] .= warp_face_nodes_to_volume_nodes * vec(curved_face_coordinates[2])
        z[:, element] .= warp_face_nodes_to_volume_nodes * vec(curved_face_coordinates[3])
    end

    md_curved = MeshData(rd, md, x, y, z) 

    # find boundary tags
    boundary_strings = unique(hmd.boundary_tags)
    deleteat!(boundary_strings, findfirst(boundary_strings .=== "---")) # remove non-boundary faces
    face_indices = [map(x -> global_face_index(x, rd.num_faces, HOHQMesh_to_StartUpDG_face_ordering), 
                        findall(hmd.boundary_tags .== name)) for name in boundary_strings]
    boundary_faces = NamedTuple(Pair.(Symbol.(boundary_strings), face_indices))

    return @set md_curved.mesh_type = HOHQMeshType(hmd, boundary_faces)
end

vertex_reordering(::Tet) = SVector(1, 2, 3, 4)
function get_HOHQMesh_ids(::Tet, curved_faces, f_HOHQMesh, f, polydeg) 
    num_face_nodes = (polydeg + 1)^2
    return (1:num_face_nodes) .+ (f_HOHQMesh - 1) * num_face_nodes #active_face_offset[f_HOHQMesh]
end

function abtors(a, b)
    s = b
    r = @. (1 - s) * 0.5 * (a + 1) - 1
    return r, s
end

# for tetrahedral meshes
function MeshData(hmd::HOHQMeshData{3}, rd::RefElemData{3, <:Tet})
    (; VXYZ, EToV) = hmd    
    md = MeshData(VXYZ, EToV[:, vertex_reordering(rd.element_type)], rd)

    # interpolation matrix from chebyshev_to_lobatto nodes
    a_chebyshev = vec([-cos(j * pi / hmd.polydeg) for j in 0:hmd.polydeg, i in 0:hmd.polydeg])
    b_chebyshev = vec([-cos(j * pi / hmd.polydeg) for i in 0:hmd.polydeg, j in 0:hmd.polydeg])
    r_chebyshev, s_chebyshev = abtors(a_chebyshev, b_chebyshev)

    chebyshev_to_lobatto = 
        vandermonde(face_type(rd.element_type), hmd.polydeg, nodes(face_type(rd.element_type), rd.N)...) / vandermonde(face_type(rd.element_type), hmd.polydeg, r_chebyshev, s_chebyshev)

    # permute rows of chebyshev_to_lobatto according to the permutation between StartUpDG rd.Fmask 
    # ordering  and nodes(face_type(Tet())) ordering
    Fmask = reshape(rd.Fmask, :, 4) # 4 = number of tet faces
    rf, sf = nodes(face_type(rd.element_type), rd.N)
    rFmask, sFmask = rd.r[Fmask[:, 4]], rd.s[Fmask[:, 4]] # use face 4 nodes for which t = -1
    permuted_face_ids = match_coordinate_vectors((rf, sf), (rFmask, sFmask))    
    chebyshev_to_lobatto = chebyshev_to_lobatto[permuted_face_ids, :]                

    warp_face_nodes_to_volume_nodes = 
        face_basis(rd.element_type, rd.N, rd.rst...) / face_basis(rd.element_type, rd.N, rd.r[rd.Fmask], rd.s[rd.Fmask], rd.t[rd.Fmask])            

    HOHQMesh_to_StartUpDG_face_ordering = get_HOHQMesh_to_StartUpDG_face_ordering(rd.element_type)
    curved_face_coordinates = ntuple(_ -> similar(reshape(md.x[rd.Fmask, 1], :, rd.num_faces)), 3)
    x, y, z = copy.(md.xyz)
    for curved_elem in hmd.curved_elements
        (; element, curved_faces, curved_edge_coordinates) = curved_elem

        # initialize face coordinates as linear coordinates
        curved_face_coordinates[1] .= reshape(md.x[rd.Fmask, element], :, rd.num_faces)
        curved_face_coordinates[2] .= reshape(md.y[rd.Fmask, element], :, rd.num_faces)
        curved_face_coordinates[3] .= reshape(md.z[rd.Fmask, element], :, rd.num_faces)

        for (f_HOHQMesh, f) in enumerate(HOHQMesh_to_StartUpDG_face_ordering)            
            if curved_faces[f_HOHQMesh] == 1
                # incorporate any permutations or face reorderings into indices
                ids_HOHQMesh = get_HOHQMesh_ids(rd.element_type, curved_faces, f_HOHQMesh, f, hmd.polydeg)

                curved_lobatto_coordinates = chebyshev_to_lobatto * curved_edge_coordinates[ids_HOHQMesh, :]
                curved_face_coordinates[1][:, f] .= curved_lobatto_coordinates[:, 1]
                curved_face_coordinates[2][:, f] .= curved_lobatto_coordinates[:, 2]
                curved_face_coordinates[3][:, f] .= curved_lobatto_coordinates[:, 3]
            end
        end        

        x[:, element] .= warp_face_nodes_to_volume_nodes * vec(curved_face_coordinates[1])
        y[:, element] .= warp_face_nodes_to_volume_nodes * vec(curved_face_coordinates[2])
        z[:, element] .= warp_face_nodes_to_volume_nodes * vec(curved_face_coordinates[3])
    end

    md_curved = MeshData(rd, md, x, y, z) 
    
    # find boundary tags
    boundary_tags = hmd.boundary_tags[:, 1:num_faces(rd.element_type)] # ignore `#undef` entries in last columns 5 and 6
    boundary_strings = unique(boundary_tags)
    deleteat!(boundary_strings, findfirst(boundary_strings .=== "---")) # remove non-boundary faces
    face_indices = [map(x -> global_face_index(x, rd.num_faces, HOHQMesh_to_StartUpDG_face_ordering), 
                        findall(boundary_tags .== name)) for name in boundary_strings]
    boundary_faces = NamedTuple(Pair.(Symbol.(boundary_strings), face_indices))

    return @set md_curved.mesh_type = HOHQMeshType(hmd, boundary_faces)
end

function read_HOHQMesh(filename::String, element_type=Quad())
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
    curved_faces::Vector{Int}
    curved_edge_coordinates::Matrix{Float64} # TODO: change this to curved_face_coordinates for 2D/3D consistency 
end

# for hexes and quads
function _read_HOHQMesh(lines, mesh_format)
    if mesh_format isa ISM_V2
        num_nodes, num_edges, nelements, polydeg = parse.(Int, split(popat!(lines, 1)))
    else
        num_nodes, nelements, polydeg = parse.(Int, split(popat!(lines, 1)))
    end
    VX, VY, VZ = ntuple(_ -> zeros(num_nodes), 3)    

    has_z_coordinate = length(split(lines[1])) == 3 
    for i in 1:num_nodes
        if has_z_coordinate
            VX[i], VY[i], VZ[i] = parse.(Float64, split(lines[i]))
        else
            VX[i], VY[i] = parse.(Float64, split(lines[i]))
        end
    end
    deleteat!(lines, 1:num_nodes)

    if mesh_format isa ISM_V2
        # ignore edge data for now
        deleteat!(lines, 1:num_edges)
    end

    curved_elements = CurvedHOHQMeshElement[]
    nvertices = length(split(lines[1]))    
    num_faces = nvertices == 4 ? 4 : 6 # in 2D, 4 faces. In 3D, 6 faces
    EToV = zeros(Int, nelements, nvertices)
    boundary_tags = Matrix{String}(undef, nelements, num_faces)
    for e in 1:nelements
        EToV[e, :] .= parse.(Int, split(lines[1]))
        curved_faces = parse.(Bool, split(lines[2]))

        deleteat!(lines, 1:2) # move onto next lines

        if all(curved_faces .== false)
            # do nothing
        else
            num_curved_faces = count(curved_faces)
            num_face_nodes = nvertices == 4 ? (polydeg + 1) : (polydeg + 1)^2
            curved_nodes = 1:num_face_nodes * num_curved_faces
            curved_edge_coordinates = mapreduce(vcat, lines[curved_nodes]) do s
                (parse.(Float64, split(s)))'
            end

            push!(curved_elements, 
                  CurvedHOHQMeshElement(e, curved_faces, curved_edge_coordinates))
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
    
    # for Tri/Tet meshes. Tets still have 6 faces each because they're 
    # constructed from a collapsed mapping of a hex.
    num_faces = nvertices == 3 ? 3 : 6 # TODO: fix when David fixes the number of tet faces

    EToV = zeros(Int, nelements, nvertices)
    boundary_tags = Matrix{String}(undef, nelements, num_faces)
    for e in 1:nelements
        EToV[e, :] .= parse.(Int, split(lines[1]))
        curved_faces = parse.(Bool, split(lines[2]))

        deleteat!(lines, 1:2) # move onto next lines

        if all(curved_faces .== false)
            # do nothing
        else
            num_curved_faces = count(curved_faces)            

            # note: HOHQMesh.jl uses (p+1)^2 nodes per face of the triangle, which 
            # correspond to collapsed coordinates on the quad face. 
            num_face_nodes = element_type isa Tri ? (polydeg + 1) : (polydeg + 1)^2 

            curved_nodes = 1:num_face_nodes * num_curved_faces
            curved_edge_coordinates = mapreduce(vcat, lines[curved_nodes]) do s
                (parse.(Float64, split(s)))'
            end

            push!(curved_elements, 
                  CurvedHOHQMeshElement(e, curved_faces, curved_edge_coordinates))
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