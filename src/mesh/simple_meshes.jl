"""
findline(word::String,lines)
## Ouput
line number of given word in file
## Notes
expexted that the word exists once in the file
"""
function findline(name::String, lines::Vector{String})
    for (i, line) in enumerate(lines)
        if line == name
            return i
        end
    end
end

"""
    MeshImportOptions
This struct allows the user to opt for supported features when importing
a gmsh 4.1 .msh file.
## Support
- grouping::Bool | On import would you like to include physical group assignements of 2D elements?
- remap\\_group\\_name::Bool | On import would you like to maintain or remap physical group ID? Remap results in groupIds in the range 1:number\\_group\\_ids.
"""
struct MeshImportOptions
    grouping::Bool
    remap_group_name::Bool
end

"""
 returns the number of elements in a .msh file of a specified dimension
 ## Notes: Gmsh includes elements in a .msh file of multiple dimensions. We want a count of how many
 2D elements are in our file. This corisponds to the number of elements in our tri mesh.
"""
function get_num_elements(lines::Vector{String}, Dim=2)
    elem_start = findline("\$Elements", lines) + 1
    temp_num_blocks = split(lines[elem_start])[1]
    num_blocks = parse(Int, temp_num_blocks)
    num_elements = 0
    block_data_line = elem_start + 1
    for i in 1:num_blocks
        temp_blockDim, _, _, temp_num_elem_in_block = split(lines[block_data_line])
        blockDim = parse(Int, temp_blockDim)
        num_elem_in_block = parse(Int, temp_num_elem_in_block)
        if blockDim == Dim
            num_elements = num_elements + num_elem_in_block
        end
        block_data_line = block_data_line + num_elem_in_block + 1
    end
    return num_elements
end

"""
remap_element_grouping!(eg::Vector{Int})
GMSH uses integers for naming conventions. This function remaps the gmsh ids to
a list of ids 1:numGroups. This just cleans up a little after Gmsh
## Example output
remap_element_grouping([16,16,17,17]) -> [1,1,2,2]
"""
function remap_element_grouping(eg::Vector{Int})
    groupids = unique(eg)
    newids = 1:length(groupids)
    map = Dict{Int,Int}()
    for (i, id) in enumerate(groupids)
        map[id] = newids[i]
    end
    @info "New group names: $map"
    for (i, element) in enumerate(eg)
        eg[i] = map[element]
    end
    return eg
end

"""
 function readGmsh2D_v4(filename)

reads triangular GMSH 2D .msh files.

# Output
This depends on if grouping is opted for or not
- returns: (VX,VY), EToV
- return:(VX,VY),EToV,grouping

# Supported formats and features:
- version 4.1
    'physical group support
    'remap group ids

## grouping application
When modeling the wave equation you might want wave speeds to vary accross your domain. By assigning Physical groups
in gmsh we can maintain such groupings upon importing the .msh file. Each imported element will be a member of a phyical group.

```julia
VXY, EToV = readGmsh2D_v4("eulerSquareCylinder2D.msh")
VXY, EToV = readGmsh2D_v4("eulerSquareCylinder2D.msh",false)
VXY, EToV, grouping = readGmsh2D_v4("eulerSquareCylinder2D.msh", true)

option = MeshOption(true)
VXY, EToV, grouping = readGmsh2D_v4("eulerSquareCylinder2D.msh", option)
```
https://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format

Notes:The version 4 format has a more detailed block data format
this leads to more complicated parser.
"""
function readGmsh2D_v4(filename::String, options::MeshImportOptions)
    @unpack grouping, remap_group_name = options

    if !isfile(filename)
        throw(ArgumentError("file $filename does not exist"))
    end

    f = open(filename)
    lines = readlines(f)

    format_line = findline("\$MeshFormat", lines) + 1
    version, _, dataSize = split(lines[format_line])
    gmsh_version = parse(Float64, version)
    group_requested_but_none_in_file = false


    if gmsh_version == 4.1
        @info "reading gmsh file with legacy ($gmsh_version) format"
    else
        @warn "Gmsh file version is: $gmsh_version consider using a different parsing fuction for this file format"
    end

    # grouping may be requested yet not be present in the file
    if isnothing(findline("\$PhysicalNames", lines)) && grouping == true
        @warn "No grouping data in file. Setting grouping option to false"
        group_requested_but_none_in_file = true
        grouping = false
        remap_group_name = false
    end

    if grouping
        # surface_data_dict = get_surface_tag(lines::Vector{String})

        # For 2D Physical group data we are only interested in groupings in surface entities
        # This will extract relevent grouping data for each surface
        # only supports one physical group per Element
        entities_start_line = findline("\$Entities", lines) + 1
        entities_block_data = split(lines[entities_start_line])
        entities_block_data = [parse(Int, c) for c in entities_block_data[1:3]]
        numPoints, numCurves, numSurfaces = entities_block_data[1:3]
        @info "Entities| points:$numPoints| curves:$numCurves| surfaces:$numSurfaces"

        surface_start_line = entities_start_line + numPoints + numCurves
        @info "Expecting Surface data to start on line $surface_start_line of .msh file"
        surfaceData = Dict{Int,Int}()
        for i in 1:numSurfaces
            surfaceInfo = split(lines[i+surface_start_line])
            surfTag, numPhysTag, PhysTag = [parse(Int, c) for c in surfaceInfo[[1, 8, 9]]]
            @assert numPhysTag == 1 "Surfaces must have one physical tag associated you have $numPhysTag"
            surfaceData[surfTag] = PhysTag
        end
        @info "Tag grouping: $surfaceData"
    end

    # VX,VY,VZ = getnodes(lines::Vector{String})
    node_start = findline("\$Nodes", lines) + 1
    temp_num_blocks, temp_Nv, _, _ = split(lines[node_start])
    num_blocks = parse(Int, temp_num_blocks)
    Nv = parse(Int, temp_Nv)
    VX, VY, VZ = ntuple(x -> zeros(Float64, Nv), 3)
    block_line = node_start + 1
    for block_idx in 1:num_blocks
        temp_num_nodes = split(lines[block_line])[4]
        num_nodes = parse(Int, temp_num_nodes)
        for block_node_idx in 1:num_nodes
            node_num_line = block_line + block_node_idx
            node_cords_line = node_num_line + num_nodes
            node_idx_global = parse(Int, lines[node_num_line])
            vals = [parse(Float64, c) for c in split(lines[node_cords_line])]
            VX[node_idx_global] = vals[1]
            VY[node_idx_global] = vals[2]
        end
        block_line = block_line + 2 * num_nodes + 1
    end

    elem_start = findline("\$Elements", lines) + 1
    temp_num_blocks, _, _, _ = split(lines[elem_start])
    num_elements = get_num_elements(lines) # Some elements in the .msh file are not our 2d mesh elements see docstring of function
    @info "File has $num_elements 2D Elements"
    num_blocks = parse(Int, temp_num_blocks)

    if grouping
        element_grouping = zeros(Int, num_elements) #element physical group assignment
    end

    EToV = zeros(Int64, num_elements, 3)
    # EToV = get_etov(lines::Vector{String})
    block_line_start = elem_start + 1
    elem_global_idx = 0
    for block in 1:num_blocks
        temp_elem_block_dim, temp_tag, _, temp_elem_in_block = split(lines[block_line_start])
        elem_block_dim = parse(Int, temp_elem_block_dim)
        num_elem_in_block = parse(Int, temp_elem_in_block)
        surface_tag = parse(Int, temp_tag)
        if elem_block_dim == 2 # only interesed in 2d triangle elements
            for e_idx in 1:num_elem_in_block
                elem_global_idx = elem_global_idx + 1
                vals = [parse(Int, c) for c in split(lines[e_idx+block_line_start])]
                _, nodeA, nodeB, nodeC = vals
                EToV[elem_global_idx, :] .= [nodeA, nodeB, nodeC]
                if grouping
                    element_grouping[elem_global_idx] = surfaceData[surface_tag]
                end
            end
        end
        block_line_start = block_line_start + num_elem_in_block + 1
    end

    EToV = EToV[:, vec([1 3 2])]  # permute for Gmsh ordering
    EToV = correct_negative_Jacobians!((VX, VY), EToV)

    if grouping
        if remap_group_name
            @info "Remapping Group Names"
            element_grouping = remap_element_grouping(element_grouping)
        end
        return (VX, VY), EToV, element_grouping
    end

    if group_requested_but_none_in_file
        return (VX, VY), EToV, zeros(Int, num_elements)
    end

    return (VX, VY), EToV
end

"""
For brevity while grouping is the only supported feature this allows
for simpler code
    example: VXY, EToV, grouping = readGmsh2D_v4("file.msh",true)
    example: VXY, EToV = readGmsh2D_v4("file.msh",false)
"""
function readGmsh2D_v4(filename::String, groupOpt::Bool=false, remap_group_name::Bool=false)
    options = MeshImportOptions(groupOpt, remap_group_name)
    return readGmsh2D_v4(filename, options)
end

"""
function readGmsh2D(filename)
reads triangular GMSH 2D file format 2.2 0 8. returns (VX, VY), EToV
# Examples
```julia
VXY, EToV = readGmsh2D("eulerSquareCylinder2D.msh")
```

https://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format-version-2-_0028Legacy_0029
"""
function readGmsh2D(filename::String)
    f = open(filename)
    lines = readlines(f)

    format_line = findline("\$MeshFormat", lines) + 1
    version, _, dataSize = split(lines[format_line])
    gmsh_version = parse(Float64, version)
    if gmsh_version == 2.2
        @info "reading gmsh file with legacy ($gmsh_version) format"
    elseif gmsh_version == 4.1
        @assert "This parsing function is not compatable with gmsh version $gmsh_version"
    else
        @warn "Gmsh file version is: $gmsh_version consider using a different parsing fuction for this file format"
    end

    node_start = findline("\$Nodes", lines) + 1
    Nv = parse(Int64, lines[node_start])
    VX, VY, VZ = ntuple(x -> zeros(Float64, Nv), 3)
    for i = 1:Nv
        vals = [parse(Float64, c) for c in split(lines[i+node_start])]
        # first entry =
        VX[i] = vals[2]
        VY[i] = vals[3]
    end

    elem_start = findline("\$Elements", lines) + 1
    K_all = parse(Int64, lines[elem_start])
    K = 0
    for e = 1:K_all
        if length(split(lines[e+elem_start])) == 8
            K = K + 1
        end
    end
    EToV = zeros(Int64, K, 3)
    sk = 1
    for e = 1:K_all
        if length(split(lines[e+elem_start])) == 8
            vals = [parse(Int64, c) for c in split(lines[e+elem_start])]
            EToV[sk, :] .= vals[6:8]
            sk = sk + 1
        end
    end

    EToV = EToV[:, vec([1 3 2])] # permute for Gmsh ordering
    EToV = correct_negative_Jacobians!((VX, VY), EToV)

    return (VX, VY), EToV
end

#     compute_triangle_area(tri)
#
# Computes the area of a triangle given `tri`, which is a tuple of three points (vectors),
# using the [Shoelace_formula](https://en.wikipedia.org/wiki/Shoelace_formula).
function compute_triangle_area(tri)
    B, A, C = tri
    return 0.5 * (A[1] * (B[2] - C[2]) + B[1] * (C[2] - A[2]) + C[1] * (A[2] - B[2]))
end

function correct_negative_Jacobians!((VX, VY), EToV)
    # detect negative Jacobians
    for e in 1:size(EToV, 1)
        v_ids = view(EToV, e, :)
        A, B, C = ((VX[v_ids[i]], VY[v_ids[i]]) for i in 1:3)
        tri = SVector{3}(A, B, C)
        area = compute_triangle_area(tri)
        # if triangle area is negative, permute the vertices
        # so the area is positive
        if area < 0
            view(EToV, e, :) .= (v_ids[2], v_ids[1], v_ids[3])
        end
    end
    return EToV
end

"""
    uniform_mesh(elem::Line,Kx)
    uniform_mesh(elem::Tri,Kx,Ky)
    uniform_mesh(elem::Quad,Kx,Ky)
    uniform_mesh(elem::Hex,Kx,Ky,Kz)
    uniform_mesh(elem, K)

Uniform `Kx` (by `Ky` by `Kz`) mesh on ``[-1,1]^d``, where `d` is the spatial dimension.
Returns `(VX,VY,VZ)`, `EToV`. When only one `K` is specified, it assumes a uniform mesh with
`K` elements in each coordinate direction.

`K` can also be specified using a keyword argument `K1D`, e.g., `uniform_mesh(elem; K1D = 16)`.
"""
function uniform_mesh(elem::Line, K1D)
    VX = collect(LinRange(-1, 1, K1D + 1))
    EToV = transpose(reshape(sort([1:K1D; 2:K1D+1]), 2, K1D))
    return (VX,), Matrix(EToV)
end

function uniform_mesh(elem::Tri, Kx, Ky)

    (VY, VX) = meshgrid(LinRange(-1, 1, Ky + 1), LinRange(-1, 1, Kx + 1))
    sk = 1
    EToV = zeros(Int, 2 * Kx * Ky, 3)
    for ey = 1:Ky
        for ex = 1:Kx
            id(ex, ey) = ex + (ey - 1) * (Kx + 1) # index function
            id1 = id(ex, ey)
            id2 = id(ex + 1, ey)
            id3 = id(ex + 1, ey + 1)
            id4 = id(ex, ey + 1)
            EToV[2*sk-1, :] = [id1 id3 id2]
            EToV[2*sk, :] = [id3 id1 id4]
            sk += 1
        end
    end
    return (VX[:], VY[:]), EToV
end

function uniform_mesh(elem::Wedge, Kx, Ky, Kz)
    (VY, VX, VZ) = meshgrid(LinRange(-1, 1, Ky + 1), LinRange(-1, 1, Kx + 1), LinRange(-1, 1, Kz + 1))
    sk = 1
    shift = (Kx+1)*(Ky+1)
    EToV = zeros(Int, 2 * Kx * Ky * Kz, 6)
    for ey = 1:Ky
        for ex = 1:Kx
            for ez = 1:Kz
                id(ex, ey, ez) = ex + (ey - 1) * (Kx + 1)  + (ez - 1) * (shift)# index function
                id1 = id(ex, ey, ez)
                id2 = id(ex + 1, ey, ez)
                id3 = id(ex, ey + 1, ez)
                id4 = id(ex + 1, ey + 1, ez)
                id5 = id(ex, ey, ez + 1)
                id6 = id(ex + 1, ey, ez + 1)
                id7 = id(ex, ey + 1, ez + 1)
                id8 = id(ex + 1, ey + 1, ez + 1)
                EToV[2*sk-1, :] = [id1 id5 id4 id8 id2 id6]
                EToV[2*sk, :] = [id1 id5 id3 id7 id4 id8]
                sk += 1
            end
        end
    end
    EToV .= EToV[:, SVector(1, 3, 5, 2, 4, 6)]
    return (VX[:], VY[:], VZ[:]), EToV
end

uniform_mesh(elem::Union{Tri,Quad}, Kx) = uniform_mesh(elem, Kx, Kx)
uniform_mesh(elem::Union{Hex, Wedge}, Kx) = uniform_mesh(elem, Kx, Kx, Kx)


function uniform_mesh(elem::Quad, Nx, Ny)

    Nxp = Nx + 1
    Nyp = Ny + 1
    Nv = convert(Int, Nxp * Nyp)
    K = convert(Int, Nx * Ny)

    x1D = LinRange(-1, 1, Nxp)
    y1D = LinRange(-1, 1, Nyp)
    x, y = meshgrid(x1D, y1D)
    I, J = meshgrid(collect(1:Nxp), collect(1:Nyp))
    inds = @. (I - 1) * Ny + (J + I - 1)
    EToV = zeros(Int, K, 4)
    k = 1
    for i = 1:Ny
        for j = 1:Nx
            EToV[k, :] = [inds[i, j] inds[i, j+1] inds[i+1, j] inds[i+1, j+1]]
            k += 1
        end
    end

    VX = x[:]
    VY = y[:]

    return (VX[:], VY[:]), EToV
end


# hack to work with old vertex ordering
function mymeshgrid(x1D, y1D, z1D)
    Np = length(x1D) * length(y1D) * length(z1D)
    x = zeros(Np)
    y = zeros(Np)
    z = zeros(Np)
    sk = 1
    for k = 1:length(z1D)
        for j = 1:length(y1D)
            for i = 1:length(x1D)
                x[sk] = x1D[i]
                y[sk] = y1D[j]
                z[sk] = z1D[k]
                sk += 1
            end
        end
    end
    return x, y, z
end

function uniform_mesh(elem::Hex, Nx, Ny, Nz)
    Nxp = Nx + 1
    Nyp = Ny + 1
    Nzp = Nz + 1
    Nv = convert(Int, Nxp * Nyp * Nzp)
    K = convert(Int, Nx * Ny * Nz)

    x1D = LinRange(-1, 1, Nxp)
    y1D = LinRange(-1, 1, Nyp)
    z1D = LinRange(-1, 1, Nzp)

    x, y, z = mymeshgrid(x1D, y1D, z1D)

    EToV = zeros(Int, K, 8)
    k = 1
    for e = 1:K
        em = e - 1
        k = div(em, (Nx * Ny))
        j = div(em - k * Nx * Ny, Nx)
        i = em % Nx

        EToV[e, 1] = i + Nxp * j + Nxp * Nyp * k
        EToV[e, 2] = i + Nxp * (j + 1) + Nxp * Nyp * k
        EToV[e, 3] = (i + 1) + Nxp * j + Nxp * Nyp * k
        EToV[e, 4] = (i + 1) + Nxp * (j + 1) + Nxp * Nyp * k
        EToV[e, 5] = i + Nxp * j + Nxp * Nyp * (k + 1)
        EToV[e, 6] = i + Nxp * (j + 1) + Nxp * Nyp * (k + 1)
        EToV[e, 7] = (i + 1) + Nxp * j + Nxp * Nyp * (k + 1)
        EToV[e, 8] = (i + 1) + Nxp * (j + 1) + Nxp * Nyp * (k + 1)
    end
    EToV = @. EToV + 1 # re-index to 1 index

    return (x[:], y[:], z[:]), EToV
end

uniform_mesh(elem::Tet, N) = uniform_mesh(elem, N, N, N)

function uniform_mesh(elem::Tet, Nx, Ny, Nz)

    Nxp = Nx + 1
    Nyp = Ny + 1
    Nzp = Nz + 1
    Nv = Nxp * Nyp * Nzp
    K = Nx * Ny * Nz

    id(i, j, k) = i + Nxp * j + Nxp * Nyp * k

    x1D = LinRange(-1, 1, Nxp)
    y1D = LinRange(-1, 1, Nyp)
    z1D = LinRange(-1, 1, Nzp)

    x, y, z = mymeshgrid(x1D, y1D, z1D)

    # 6 tets per cube
    EToV = zeros(Int, 6 * K, 4)
    ee = 1
    for e = 1:K
        em = e - 1
        k = div(em, (Nx * Ny))
        j = div(em - k * Nx * Ny, Nx)
        i = em % Nx

        EToV[ee, 1] = id(i, j, k)
        EToV[ee, 2] = id(i + 1, j, k)
        EToV[ee, 3] = id(i + 1, j + 1, k)
        EToV[ee, 4] = id(i, j, k + 1)

        EToV[ee+1, 1] = id(i + 1, j, k)
        EToV[ee+1, 2] = id(i + 1, j, k + 1)
        EToV[ee+1, 3] = id(i + 1, j + 1, k)
        EToV[ee+1, 4] = id(i, j, k + 1)

        EToV[ee+2, 1] = id(i, j, k + 1)
        EToV[ee+2, 2] = id(i + 1, j + 1, k + 1)
        EToV[ee+2, 3] = id(i + 1, j, k + 1)
        EToV[ee+2, 4] = id(i + 1, j + 1, k)

        EToV[ee+3, 1] = id(i, j + 1, k + 1)
        EToV[ee+3, 2] = id(i + 1, j + 1, k + 1)
        EToV[ee+3, 3] = id(i, j, k + 1)
        EToV[ee+3, 4] = id(i + 1, j + 1, k)

        EToV[ee+4, 1] = id(i, j + 1, k)
        EToV[ee+4, 2] = id(i, j + 1, k + 1)
        EToV[ee+4, 3] = id(i, j, k + 1)
        EToV[ee+4, 4] = id(i + 1, j + 1, k)

        EToV[ee+5, 1] = id(i, j, k)
        EToV[ee+5, 2] = id(i + 1, j + 1, k)
        EToV[ee+5, 3] = id(i, j + 1, k)
        EToV[ee+5, 4] = id(i, j, k + 1)

        ee += 6
    end

    EToV = @. EToV + 1 # re-index to 1 index
    return (x[:], y[:], z[:]), EToV
end

# keyword argument version
uniform_mesh(elem::AbstractElemShape; K1D) = uniform_mesh(elem, K1D)
