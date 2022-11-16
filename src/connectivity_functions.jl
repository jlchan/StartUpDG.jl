"""
    connect_mesh(EToV,fv)

Initialize element connectivity matrices, element to element and element to face
connectivity.

Inputs:
- `EToV` is a `K` by `Nv` matrix whose rows identify the `Nv` vertices
which make up an element.
- `fv` (an array of arrays containing unordered indices of face vertices).

Output: `FToF`, `length(fv)` by `K` index array containing face-to-face connectivity.
"""
function connect_mesh(EToV, fv)
    Nfaces = length(fv)
    K = size(EToV, 1)
    elem_types = [Tri(), Quad(), Tet(), Hex(), Wedge(), Pyr()]
    # sort and find matches
    fnodes = Vector{eltype(first(EToV))}[]
    for e in 1:K
        vertex_ids = EToV[e, :]
        element_type = element_type_from_num_vertices(elem_types, length(vertex_ids))
        for face in 1:num_faces(element_type)
            push!(fnodes, EToV[e, fv[face]])
        end
    end
    
    #sort and find matches
    p = sortperm(fnodes) # sorts by lexicographic ordering by default
    fnodes = fnodes[p, :]

    FToF = reshape(collect(1:Nfaces * K), Nfaces, K)
    for f = 1:size(fnodes, 1) - 1
        if fnodes[f, :]==fnodes[f + 1, :]
            f1 = FToF[p[f]]
            f2 = FToF[p[f + 1]]
            FToF[p[f]] = f2
            FToF[p[f + 1]] = f1
        end
    end
    return FToF
end

# returns back `p` such that `u[p] == v` or false
# u = tuple of vectors containing coordinates
function match_coordinate_vectors(u, v; tol = 100 * eps())
    p = zeros(Int, length(first(u)))
    return match_coordinate_vectors!(p, u, v)
end
function match_coordinate_vectors!(p, u, v; tol = 100 * eps())
    for (i, u_i) in enumerate(zip(u...))
        for (j, v_i) in enumerate(zip(v...))
            if norm(u_i .- v_i) < tol 
                p[i] = j
            end
        end
    end
    return p 
end

@inline function n_face_nodes(rd::RefElemData, face)
    return first(rd.element_type.node_ids_by_face[face]), last(rd.element_type.node_ids_by_face[face])
end

"""
    build_node_maps(FToF, Xf)

Intialize the connectivity table along all edges and boundary node tables of all
elements. `mapM` - map minus (interior). `mapP` - map plus (exterior).

`Xf = (xf, yf, zf)` and `FToF` is size `(Nfaces * K)` and `FToF[face]` = face neighbor

`mapM`, `mapP` are Vectors of Vectors, of the length Nfaces, where each Vector
is of length Nfp of the corresponding

# Examples
```julia
julia> mapM, mapP, mapB = build_node_maps(FToF, (xf, yf))
```
"""
function build_node_maps(FToF, Xf; tol = 1e-12)

    # total number of faces 
    num_faces_total = length(FToF)    
    
    # num face points: assumes all faces have the same number of nodes
    num_nodes_per_face = length(Xf[1]) ÷ num_faces_total 

    # number nodes consecutively
    mapM = reshape(collect(1:length(Xf[1])), num_nodes_per_face, num_faces_total);
    mapP = copy(mapM);

    # reshape to be a face-first array
    Xf = reshape.(Xf, num_nodes_per_face, num_faces_total)

    p = zeros(Int, num_nodes_per_face)
    for (f1, f2) in enumerate(FToF)
        face_1_coordinates = view.(Xf, :, f1)
        face_2_coordinates = view.(Xf, :, f2)        
        match_coordinate_vectors!(p, face_1_coordinates, face_2_coordinates, tol=tol)
        @. mapP[:, f1] = mapM[p, f2]
    end
    mapB = findall(vec(mapM) .== vec(mapP))
    return mapM[:], mapP[:], mapB[:]
end

function build_node_maps(FToF, EToV, rd, Xf...; tol = 1e-12)
    # If the element has uniform surface elements, one can use the
    # original implementation
    # TODO: This can be done more beautiful/better with some Julia-magic
    if rd.element_type ∈ [Line(), Tri(), Quad(), Tet(), Hex()]
        return build_node_maps(FToF, Xf)
    end

    # Store the coordinates of the nodes of every vertex for each face seperatly
    # Each dimension gets its own vector
    Xf_sorted = [Vector{Float64}[], Vector{Float64}[], Vector{Float64}[]]
    dims = length(Xf)

    # Possible types of elements used in a mesh
    elem_types = [Tri(), Quad(), Tet(), Hex(), Wedge(), Pyr()]

    # Degree of the polynomial
    @unpack N = rd

    # Get the number of elements
    num_elements = size(EToV, 1)
   
    # Construct the maps that will be filled with the ids of the connected face-vertices
    mapM = Vector{eltype(first(FToF))}[]

    # Helper variable to keep track of number of nodes.
    number_of_passed_face_nodes = 0
    # Initialize with identity-mapping
    for elem in 1:num_elements
        element_type = element_type_from_num_vertices(elem_types, length(EToV[elem, :]))
        n_faces = num_faces(element_type)
        for face in 1:n_faces
            first_face_node, last_face_node = n_face_nodes(rd, face)
            for i in 1:dims
                push!(Xf_sorted[i], Xf[i][first_face_node: last_face_node, elem])
            end
            num_face_nodes = last_face_node - first_face_node + 1
            # Initialize maps as identity maps
            push!(mapM, collect((number_of_passed_face_nodes+1):(number_of_passed_face_nodes + num_face_nodes)))
            # Update to number of passed nodes
            number_of_passed_face_nodes += num_face_nodes
        end
    end

    mapP = copy(mapM);

    # Iterate over the Face to face connectivity and find out which
    # vertex of face f1 corresponds to a vertex of face2 (face1•p = face2)
    for (f1, f2) in enumerate(FToF)
        face_1_coords = Vector{Float64}[]
        face_2_coords = Vector{Float64}[]
        # Get the coordinates of the faces
        for i in 1:dims
            push!(face_1_coords, Xf_sorted[i][f1])
            push!(face_2_coords, Xf_sorted[i][f2])
        end
        # Compute the permutation of nodes (face )
        p = match_coordinate_vectors(face_1_coords, face_2_coords)
        mapP[f1] = mapM[f2][p]
    end
    mapB = map(x -> x[1], findall(@. mapM[:]==mapP[:]))
    mapM_flat = []
    mapP_flat = []
    for i in eachindex(mapM)
        append!(mapM_flat, mapM[i])
        append!(mapP_flat, mapP[i])
    end
    return mapM_flat, mapP_flat, mapB
end

"""
    make_periodic(md::MeshData{Dim}, is_periodic...) where {Dim}
    make_periodic(md::MeshData{Dim}, is_periodic = ntuple(x->true,Dim)) where {Dim}
    make_periodic(md::MeshData, is_periodic = true)

Returns new MeshData such that the node maps `mapP` and face maps `FToF` are now periodic.
Here, `is_periodic` is a tuple of `Bool` indicating whether or not to impose periodic
BCs in the `x`,`y`, or `z` coordinate.
"""
make_periodic(md::MeshData{Dim}, rd::RefElemData, is_periodic::Bool = true) where {Dim} = 
    make_periodic(md, rd, ntuple(_->is_periodic, Dim)) 

function make_periodic(md::MeshData{Dim}, rd::RefElemData, is_periodic::NTuple{Dim, Bool}) where {Dim, Bool}
    @unpack mapM, mapP, mapB, xyzf, FToF = md
    NfacesTotal = length(FToF)
    FToF_periodic = copy(FToF)
    mapP_periodic, mapB_periodic = build_periodic_boundary_maps!(xyzf...,is_periodic...,NfacesTotal,
                                          mapM, mapP, mapB, FToF_periodic, rd)
    #display(mapP_periodic)
    #mapP_periodic = copy(mapP)
    #mapP_periodic[mapB] = mapPB
    #mapB_periodic = mapB[mapPB .== mapP[mapB]] # keep only non-periodic boundary nodes
    #display(mapB_periodic)
    return setproperties(md, (; mapB=mapB_periodic, mapP = mapP_periodic, 
                                FToF = FToF_periodic, is_periodic = is_periodic)) # from Setfield.jl    
end

# specializes to 1D - periodic = find min/max indices of xf and reverse their order
function make_periodic(md::MeshData{1, Tv, Ti}, rd::RefElemData, is_periodic::Bool = true) where {Tv, Ti}

    if is_periodic == true
        @unpack mapP, mapB, xf, FToF = md
        mapPB = argmax(vec(xf)),argmin(vec(xf))
        mapP_periodic = copy(mapP)
        mapP_periodic[mapB] .= mapPB
        FToF_periodic = copy(FToF)
        FToF_periodic[[1,length(FToF)]] .= mapPB
        mapB_periodic = Ti[] 
        return setproperties(md, (; mapB = mapB_periodic, mapP = mapP_periodic,
                                    FToF = FToF_periodic, is_periodic = (true,))) # from Setfield.jl
    end
end

# Helper functions for `make_nodemaps_periodic!`, 2D version which modifies FToF.
function build_periodic_boundary_maps!(xf, yf, is_periodic_x, is_periodic_y,
                                       NfacesTotal, mapM, mapP, mapB, FToF, rd)

    # find boundary faces (e.g., when FToF[f] = f)
    Flist = 1:length(FToF)
    Bfaces = findall(vec(FToF) .== Flist)

    xb = xf[mapB]

    xb, yb = xf[mapB], yf[mapB]
    Nfp = convert(Int, length(xf) / NfacesTotal)
    Nbfaces = convert(Int, length(xb) / Nfp)
    xb, yb = reshape(xb, Nfp, Nbfaces), reshape(yb, Nfp, Nbfaces)

    # compute centroids of faces
    xc = vec(sum(xb, dims=1) / Nfp)
    yc = vec(sum(yb, dims=1) / Nfp)
    mapMB = reshape(mapM[mapB], Nfp, Nbfaces)
    mapPB = reshape(mapP[mapB], Nfp, Nbfaces)

    xmin, xmax = extrema(xc)
    ymin, ymax = extrema(yc)

    LX, LY = map((x -> x[2] - x[1]) ∘ extrema, (xf, yf))
    NODETOL = 100 * max(eps.((LX, LY))...)
    if abs(abs(xmax - xmin) - LX) > NODETOL && is_periodic_x
        error("periodicity requested in x, but LX = $LX while abs(xmax-xmin) = $(abs(xmax-xmin))")
    end
    if abs(abs(ymax - ymin) - LY) > NODETOL && is_periodic_y
        error("periodicity requested in y, but LX = $LY while abs(ymax-ymin) = $(abs(ymax-ymin))")
    end

    # determine which faces lie on x and y boundaries
    Xscale, Yscale = max(1, LX), max(1, LY)        
    yfaces = map(x -> x[1], findall(@. (@. abs(yc - ymax) < NODETOL * Yscale) | (@. abs(yc - ymin) < NODETOL * Yscale)))
    xfaces = map(x -> x[1], findall(@. (@. abs(xc - xmax) < NODETOL * Xscale) | (@. abs(xc - xmin) < NODETOL * Xscale)))

    D = zeros(eltype(xb), size(xb,1), size(xb,1))
    ids = zeros(Int,size(xb,1))    
    if is_periodic_y # find matches in y faces
        for i in yfaces, j in yfaces
            if i!=j
                if abs(xc[i] - xc[j]) < NODETOL * Xscale && abs(abs(yc[i] - yc[j]) - LY) < NODETOL * Yscale
                    # create distance matrix
                    @. D = abs(xb[:, i] - xb[:, j]')
                    map!(x -> x[1], ids, findall(@. D < NODETOL * Xscale))
                    @. mapPB[:, i] = mapMB[ids, j]

                    FToF[Bfaces[i]] = Bfaces[j]
                end
            end
        end
    end

    if is_periodic_x # find matches in x faces
        for i in xfaces, j in xfaces
            if i!=j
                if abs(yc[i] - yc[j]) < NODETOL * Yscale && abs(abs(xc[i] - xc[j]) - LX) < NODETOL * Xscale
                    @. D = abs(yb[:,i] - yb[:,j]')
                    map!(x->x[1], ids, findall(@. D < NODETOL * Yscale))
                    mapPB[:, i] = mapMB[ids, j]

                    FToF[Bfaces[i]] = Bfaces[j]
                end
            end
        end
    end

    return mapPB[:]
end


# 3D version of build_periodic_boundary_maps, modifies FToF
function build_periodic_boundary_maps!(xf, yf, zf,
                                       is_periodic_x, is_periodic_y, is_periodic_z,
                                       NfacesTotal, mapM, mapP, mapB, FToF, rd)

    # find boundary faces (e.g., when FToF[f] = f)
    Flist = 1:length(FToF)
    Bfaces = findall(vec(FToF) .== Flist)

    # Get the local node ids for each face of a prism
    (; node_ids_by_face) = rd.element_type
    # Total number of prisms in the mesh
    # length mapM = total number of face quadrature points
    # rd.Nfq number of face quadrature points per elemet
    num_elem = Int(length(mapM)/rd.Nfq)

    # map each face to an element
    face_to_elem = Int.(ceil.(mapB/5))

    # reshape mapM and mapP into matrices to make it more easy to access the 
    # face quadrature points of each element
    mat_mapM = reshape(mapM, rd.Nfq, num_elem)
    mat_mapP = reshape(mapP, rd.Nfq, num_elem)

    # map the globale face_id to the face number of each prism
    local_face_ids = [((mapB[i]-1)%5 + 1) for i in eachindex(mapB)]
    # for each boundary face store the face quadrature nodes
    Bface_node_ids_local = [node_ids_by_face[local_face_ids[i]] for i in eachindex(mapB)]
    # store the coords of each quadrature point on the boundary faces
    xb = [xf[Bface_node_ids_local[i], face_to_elem[i]] for i in eachindex(mapB)]
    yb = [yf[Bface_node_ids_local[i], face_to_elem[i]] for i in eachindex(mapB)]
    zb = [zf[Bface_node_ids_local[i], face_to_elem[i]] for i in eachindex(mapB)]

    mapMB = [mat_mapM[Bface_node_ids_local[i], face_to_elem[i]] for i in eachindex(mapB)]
    mapPB = [mat_mapP[Bface_node_ids_local[i], face_to_elem[i]] for i in eachindex(mapB)]


    # centroids of each boundary face
    xc = sum.(xb) ./ length.(xb)
    yc = sum.(yb) ./ length.(yb)
    zc = sum.(zb) ./ length.(zb)

    # max and min coordiantes of the centroids
    xmin, xmax = extrema(xc)
    ymin, ymax = extrema(yc)
    zmin, zmax = extrema(zc)

    # distance between the extrema
    LX, LY, LZ = map((x -> x[2] - x[1]) ∘ extrema, (xf, yf, zf))
    NODETOL = 100 * max(eps.((LX, LY, LZ))...)

    # Check if the construction of the periodicity is possible.
    if abs(abs(xmax - xmin) - LX) > NODETOL && is_periodic_x
        error("periodicity requested in x, but LX = $LX while abs(xmax-xmin) for centroids = $(abs(xmax-xmin))")
    end
    if abs(abs(ymax - ymin) - LY) > NODETOL && is_periodic_y
        error("periodicity requested in y, but LY = $LY while abs(xmax-xmin) for centroids = $(abs(ymax-ymin))")
    end
    if abs(abs(zmax - zmin) - LZ) > NODETOL && is_periodic_z
        error("periodicity requested in z, but LZ = $LZ while abs(xmax-xmin) for centroids = $(abs(zmax-zmin))")
    end

    # determine which faces lie on x and y boundaries
    xfaces = map(x -> x[1], findall(@. (@. abs(xc - xmax) < NODETOL * LX) | (@. abs(xc - xmin) < NODETOL * LX)))
    yfaces = map(x -> x[1], findall(@. (@. abs(yc - ymax) < NODETOL * LY) | (@. abs(yc - ymin) < NODETOL * LY)))
    zfaces = map(x -> x[1], findall(@. (@. abs(zc - zmax) < NODETOL * LZ) | (@. abs(zc - zmin) < NODETOL * LZ)))

    if is_periodic_x # find matches in x faces
        for i in xfaces, j in xfaces
            if i!=j
                if abs(yc[i] - yc[j]) < NODETOL * LY && abs(zc[i] - zc[j]) < NODETOL * LZ && abs(abs(xc[i] - xc[j]) - LX) < NODETOL * LX
                    # Compute the permutation of the points
                    f1_coords = [yb[i], zb[i]]
                    f2_coords = [yb[j], zb[j]]
                    p = match_coordinate_vectors(f1_coords, f2_coords, tol = NODETOL * LY)
                    mapPB[i] = mapMB[j][p]
                    FToF[Bfaces[i]] = Bfaces[j]
                end
            end
        end
    end

    # find matches in y faces
    if is_periodic_y
        for i in yfaces, j = yfaces
            if i!=j
                if abs(xc[i] - xc[j]) < NODETOL * LX && abs(zc[i] - zc[j]) < NODETOL * LZ && abs(abs(yc[i] - yc[j]) - LY) < NODETOL * LY
                    # Compute the permutation of the points
                    f1_coords = [xb[i], zb[i]]
                    f2_coords = [xb[j], zb[j]]
                    p = match_coordinate_vectors(f1_coords, f2_coords, tol = NODETOL * LX)
                    mapPB[i] = mapMB[j][p]
                    FToF[Bfaces[i]] = Bfaces[j]
                end
            end
        end
    end

    # find matches in z faces
    if is_periodic_z
        for i in zfaces, j in zfaces
            if i!=j
                if abs(xc[i] - xc[j]) < NODETOL * LX && abs(yc[i] - yc[j]) < NODETOL * LY && abs(abs(zc[i] - zc[j]) - LZ) < NODETOL * LZ
                    # Compute the permutation of the points
                    f1_coords = [xb[i], yb[i]]
                    f2_coords = [xb[j], yb[j]]
                    p = match_coordinate_vectors(f1_coords, f2_coords, tol = NODETOL * LX)
                    mapPB[i] = mapMB[j][p]
                    FToF[Bfaces[i]] = Bfaces[j]
                end
            end
        end 
    end

    

    for i in eachindex(mapPB)
        elem = face_to_elem[i]
        mat_mapP[Bface_node_ids_local[i], face_to_elem[i]] = mapPB[i]
    end

    mapB = map(x -> x[1], findall(@. mat_mapP[:]==mat_mapM[:]))

    return mat_mapP[:], mapB
end

"""
# 3D version of build_periodic_boundary_maps, modifies FToF
function build_periodic_boundary_maps!(xf, yf, zf,
                                       is_periodic_x, is_periodic_y, is_periodic_z,
                                       NfacesTotal, mapM, mapP, mapB, FToF, rd)

    # find boundary faces (e.g., when FToF[f] = f)
    Flist = 1:length(FToF)
    Bfaces = findall(vec(FToF) .== Flist)


    xb,yb,zb = xf[mapB],yf[mapB],zf[mapB]
    Nbfaces = length(mapB)
    Nfp = length(xf) ÷ NfacesTotal
    Nbfaces = length(xb) ÷ Nfp
    xb, yb, zb = (x->reshape(x, Nfp, Nbfaces)).((xb, yb, zb))

    # compute centroids of faces
    xc = vec(sum(xb, dims=1) / Nfp)
    yc = vec(sum(yb, dims=1) / Nfp)
    zc = vec(sum(zb, dims=1) / Nfp)

    mapMB = reshape(mapM[mapB], Nfp, Nbfaces)
    mapPB = reshape(mapP[mapB], Nfp, Nbfaces)
    
    xmin, xmax = extrema(xc)
    ymin, ymax = extrema(yc)
    zmin, zmax = extrema(zc)

    LX, LY, LZ = map((x -> x[2] - x[1]) ∘ extrema, (xf, yf, zf))
    NODETOL = 100 * max(eps.((LX, LY, LZ))...)
    if abs(abs(xmax - xmin) - LX) > NODETOL && is_periodic_x
        error("periodicity requested in x, but LX = $LX while abs(xmax-xmin) for centroids = $(abs(xmax-xmin))")
    end
    if abs(abs(ymax - ymin) - LY) > NODETOL && is_periodic_y
        error("periodicity requested in y, but LY = $LY while abs(xmax-xmin) for centroids = $(abs(ymax-ymin))")
    end
    if abs(abs(zmax - zmin) - LZ) > NODETOL && is_periodic_z
        error("periodicity requested in z, but LZ = $LZ while abs(xmax-xmin) for centroids = $(abs(zmax-zmin))")
    end

    # determine which faces lie on x and y boundaries
    xfaces = map(x -> x[1], findall(@. (@. abs(xc - xmax) < NODETOL * LX) | (@. abs(xc - xmin) < NODETOL * LX)))
    yfaces = map(x -> x[1], findall(@. (@. abs(yc - ymax) < NODETOL * LY) | (@. abs(yc - ymin) < NODETOL * LY)))
    zfaces = map(x -> x[1], findall(@. (@. abs(zc - zmax) < NODETOL * LZ) | (@. abs(zc - zmin) < NODETOL * LZ)))

    D = zeros(eltype(xb), size(xb,1), size(xb,1))
    ids = zeros(Int, size(xb, 1))
    if is_periodic_x # find matches in x faces
        for i in xfaces, j in xfaces
            if i!=j
                if abs(yc[i] - yc[j]) < NODETOL * LY && abs(zc[i] - zc[j]) < NODETOL * LZ && abs(abs(xc[i] - xc[j]) - LX) < NODETOL * LX
                    # create distance matrix
                    @. D = abs(yb[:, i] - yb[:,j]') + abs(zb[:, i] - zb[j]')
                    map!(x->x[1], ids[i], findall(@. D < NODETOL * LY))
                    @. mapPB[:,i] = mapMB[ids, j]
                    FToF[Bfaces[i]] = Bfaces[j]
                end
            end
        end
    end

    # find matches in y faces
    if is_periodic_y
        for i in yfaces, j = yfaces
            if i!=j
                if abs(xc[i] - xc[j]) < NODETOL * LX && abs(zc[i] - zc[j]) < NODETOL * LZ && abs(abs(yc[i] - yc[j]) - LY) < NODETOL * LY
                    @. D = abs(xb[:, i] - xb[:, j]') + abs(zb[:,i] - zb[:,j]')
                    map!(x->x[1], ids[i], findall(@. D < NODETOL * LX))
                    mapPB[:,i] = mapMB[ids, j]
                    FToF[Bfaces[i]] = Bfaces[j]
                end
            end
        end
    end

    # find matches in z faces
    if is_periodic_z
        for i in zfaces, j in zfaces
            if i!=j
                if abs(xc[i] - xc[j]) < NODETOL * LX && abs(yc[i] - yc[j]) < NODETOL * LY && abs(abs(zc[i] - zc[j]) - LZ) < NODETOL * LZ
                    @. D = abs(xb[:,i] - xb[:,j]') + abs(yb[:,i] - yb[:,j]')
                    map!(x->x[1], ids[i], findall(@. D < NODETOL * LX))
                    mapPB[:,i] = mapMB[ids, j]
                    FToF[Bfaces[i]] = Bfaces[j]
                end
            end
        end
    end
    return mapPB[:]
end
"""