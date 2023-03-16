struct ISM end
struct ISM_V2 end

struct HOHQMeshData{TV, TE, TC, TB}
    VXYZ::TV
    EToV::TE
    curved_elements::TC
    boundary_tags::TB
end

function read_HOHQMesh(filename::String)
    f = open(filename)
    lines = readlines(f)

    if contains(lines[1], "ISM-V2")
        format = ISM_V2()
        deleteat!(lines, 1)
    else
        format = ISM()
    end

    return HOHQMeshData(_read_HOHQMesh(lines, format)...)
end

# container for curved data for a HOHQMesh element
struct CurvedHOHQMeshElement
    element::Int
    curved_edges::Vector{Int}
    curved_edge_coordinates::Matrix{Float64}
end

function _read_HOHQMesh(lines, ::ISM_V2)    
    num_nodes, num_edges, nelements, polydeg = parse.(Int, split(popat!(lines, 1)))
    VX, VY, VZ = ntuple(_ -> zeros(num_nodes), 3)
    for i in 1:num_nodes
        VX[i], VY[i], VZ[i] = parse.(Float64, split(lines[i]))
    end
    deleteat!(lines, 1:num_nodes)

    # ignore edges for now 
    deleteat!(lines, 1:num_edges)

    curved_elements = CurvedHOHQMeshElement[]
    nvertices = length(split(lines[1]))    
    num_faces = nvertices == 4 ? 4 : 6 # in 2D, 4 faces. In 3D, 6 faces
    EToV = zeros(nelements, nvertices)
    boundary_tags = Matrix{String}(undef, nelements, num_faces)
    for e in 1:nelements
        EToV[e, :] .= parse.(Int, split(lines[1]))
        curved_edges = parse.(Int, split(lines[2]))

        deleteat!(lines, 1:2) # move onto next lines

        if all(curved_edges .== 0)
            # do nothing
        else
            curved_edge_coordinates = mapreduce(vcat, lines[1:polydeg+1]) do s
                (parse.(Float64, split(s)))'
            end
            push!(curved_elements, 
                  CurvedHOHQMeshElement(e, curved_edges, curved_edge_coordinates))
            deleteat!(lines, 1:polydeg+1) # move on to next lines
        end

        tags = split(popat!(lines, 1))
        for i in eachindex(tags)
            boundary_tags[e, i] = tags[i]
        end        
    end
    return (VX, VY, VZ), EToV, curved_elements, boundary_tags
end

function _read_HOHQMesh(lines, ::ISM)    
    num_nodes, nelements, polydeg = parse.(Int, split(popat!(lines, 1)))
    VX, VY, VZ = ntuple(_ -> zeros(num_nodes), 3)
    for i in 1:num_nodes
        VX[i], VY[i], VZ[i] = parse.(Float64, split(lines[i]))
    end
    deleteat!(lines, 1:num_nodes)

    curved_elements = CurvedHOHQMeshElement[]
    nvertices = length(split(lines[1]))    
    num_faces = nvertices == 4 ? 4 : 6 # in 2D, 4 faces. In 3D, 6 faces
    EToV = zeros(nelements, nvertices)
    boundary_tags = Matrix{String}(undef, nelements, num_faces)
    for e in 1:nelements
        EToV[e, :] .= parse.(Int, split(lines[1]))
        curved_edges = parse.(Int, split(lines[2]))

        deleteat!(lines, 1:2) # move onto next lines

        if all(curved_edges .== 0)
            # do nothing
        else
            curved_edge_coordinates = mapreduce(vcat, lines[1:polydeg+1]) do s
                (parse.(Float64, split(s)))'
            end
            push!(curved_elements, 
                  CurvedHOHQMeshElement(e, curved_edges, curved_edge_coordinates))
            deleteat!(lines, 1:polydeg+1) # move on to next lines
        end

        tags = split(popat!(lines, 1))
        for i in eachindex(tags)
            boundary_tags[e, i] = tags[i]
        end        
    end
    return (VX, VY, VZ), EToV, curved_elements, boundary_tags
end

        