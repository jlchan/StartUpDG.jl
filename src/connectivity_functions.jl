"""
    connect_mesh(EToV,fv)

Inputs:
- `EToV` is a `num_elements` by `Nv` matrix whose rows identify the `Nv` vertices
which make up one of the `num_elements` elements.
- `fv` (an array of arrays containing unordered indices of face vertices).

Output: `FToF`, an `length(fv)` by `num_elements` index array containing 
face-to-face connectivity.
"""
function connect_mesh(EToV, fv)
    num_faces = length(fv)
    num_elements = size(EToV, 1) 

    # create list of faces. Each face is identified by the 
    # vertex nodes indices associated with that face.
    fnodes = Vector{eltype(first(EToV))}[]
    for e in 1:num_elements
        for face_indices in fv
            push!(fnodes, sort(EToV[e, face_indices]))
        end
    end
    
    #sort and find matches
    p = sortperm(fnodes) # sorts by lexicographic ordering by default
    fnodes = fnodes[p, :]

    FToF = reshape(collect(1:num_faces * num_elements), num_faces, num_elements)
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
function match_coordinate_vectors(u, v; tol = 10 * eps())
    p = zeros(Int, length(first(u)))
    return match_coordinate_vectors!(p, u, v; tol)
end
function match_coordinate_vectors!(p, u, v; tol = 10 * eps())    
    for (i, u_i) in enumerate(zip(u...))
        for (j, v_i) in enumerate(zip(v...))
            # Checks if the points are close relative to the magnitude of the coordinates.
            # Scaling by length(first(u)) relaxes this bound for high degree polynomials, 
            # which incur slightly more roundoff error. 
            if norm(u_i .- v_i) < tol * length(first(u)) * max(one(eltype(u_i)), norm(u_i), norm(v_i))
                p[i] = j
            end
        end
    end
    return p 
end

"""
    build_node_maps(FToF, Xf)

Intialize the connectivity table along all edges and boundary node tables of all
elements. `mapM` - map minus (interior). `mapP` - map plus (exterior).

`Xf = (xf, yf, zf)` and `FToF` is size `(Nfaces * K)` and `FToF[face]` = face neighbor

`mapM`, `mapP` are size `Nfp` x `(Nfaces*K)`

# Examples
```julia
julia> mapM, mapP, mapB = build_node_maps(FToF, (xf, yf))
```
"""
function build_node_maps(FToF, Xf; tol = 100 * eps())

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
    return mapM, mapP, mapB
end

# Calls the simpler version of `build_node_maps` for element types that 
# have only one type of face (e.g., Tet has only Tri faces). 
build_node_maps(rd::RefElemData{3, <:Union{Tet, Hex}}, FToF, Xf; kwargs...) = 
    build_node_maps(FToF, Xf; kwargs...) 

# Specialized for elements with multiple types of faces such as the Wedge and Pyramid. 
# We need to pass in `rd::RefElemData` because it contains information necessary to 
# distinguish what nodes lie on a triangular face vs a quadrilateral face (contained 
# in rd.element_type.node_ids_by_face.
# !!! note: this version infers `num_elements` from the dimensions of `FToF::Matrix`. 
# !!! this will not work for adaptive meshes, where `FToF` will be a `Vector`.
function build_node_maps(rd::RefElemData{3, <:Union{Wedge, Pyr}}, FToF, Xf; tol = 100 * eps())    

    _, num_elements = size(FToF)
    @unpack node_ids_by_face = rd.element_type

    mapM = Vector{eltype(FToF)}[]
    offset = zero(eltype(FToF))
    for e in 1:num_elements
        for f in 1:rd.num_faces
            push!(mapM, node_ids_by_face[f] .+ offset)
        end
        offset += rd.Nfq
    end
    mapP = copy(mapM)    

    # create list of face indices
    for (f, fnbr) in enumerate(FToF)
        face_indices = mapM[f]
        nbr_face_indices = mapM[fnbr]

        # TODO: can preallocate `face_coordinates`, etc.
        face_coordinates = map(x->getindex(x, face_indices), Xf)
        nbr_face_coordinates = map(x->getindex(x, nbr_face_indices), Xf)

        p = match_coordinate_vectors(face_coordinates, nbr_face_coordinates; tol)
        mapP[f][p] .= nbr_face_indices
    end

    mapM = vcat(mapM...)
    mapP = vcat(mapP...)
    mapB = findall(vec(mapM) .== vec(mapP))
    return mapM, mapP, mapB
end

"""
    make_periodic(md::MeshData{Dim}, is_periodic...) where {Dim}
    make_periodic(md::MeshData{Dim}, is_periodic = ntuple(x->true,Dim)) where {Dim}
    make_periodic(md::MeshData, is_periodic = true)

Returns new MeshData such that the node maps `mapP` and face maps `FToF` are now periodic.
Here, `is_periodic` is a tuple of `Bool` indicating whether or not to impose periodic
BCs in the `x`,`y`, or `z` coordinate.
"""
make_periodic(md::MeshData{Dim}, is_periodic::Bool = true; kwargs...) where {Dim} = 
    make_periodic(md, ntuple(_->is_periodic, Dim); kwargs...) 

function make_periodic(md::MeshData{Dim}, is_periodic::NTuple{Dim, Bool}; kwargs...) where {Dim, Bool}

    @unpack mapM, mapP, mapB, xyzf, FToF = md
    NfacesTotal = length(FToF)
    FToF_periodic = copy(FToF)
    mapPB = build_periodic_boundary_maps!(xyzf...,is_periodic...,NfacesTotal,
                                          mapM, mapP, mapB, FToF_periodic; kwargs...)
    mapP_periodic = copy(mapP)
    mapP_periodic[mapB] = mapPB
    mapB_periodic = mapB[mapPB .== mapP[mapB]] # keep only non-periodic boundary nodes
    return setproperties(md, (; mapB=mapB_periodic, mapP = mapP_periodic, 
                                FToF = FToF_periodic, is_periodic = is_periodic)) # from Setfield.jl    
end

# specializes to 1D - periodic = find min/max indices of xf and reverse their order
function make_periodic(md::MeshData{1, Tv, Ti}, is_periodic::Bool = true; kwargs...) where {Tv, Ti}

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
                                       NfacesTotal, mapM, mapP, mapB, FToF; 
                                       tol = 100 * eps())

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
    #NODETOL = 100 * max(eps.((LX, LY))...)
    NODETOL = tol * max(LX, LY)
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
                                       NfacesTotal, mapM, mapP, mapB, FToF;
                                       tol = 100 * eps())

    # find boundary faces (e.g., when FToF[f] = f)
    Flist = 1:length(FToF)
    Bfaces = findall(vec(FToF) .== Flist)

    xb,yb,zb = xf[mapB],yf[mapB],zf[mapB]
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
    NODETOL = tol * max(LX, LY, LZ)
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
                    @. D = abs(yb[:,i] - yb[:,j]') + abs(zb[:,i] - zb[:,j]')
                    map!(x->x[1], ids, findall(@. D < NODETOL * LY))
                    @. mapPB[:,i] = mapMB[ids,j]

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
                    @. D = abs(xb[:,i] - xb[:,j]') + abs(zb[:,i] - zb[:,j]')
                    map!(x->x[1], ids, findall(@. D < NODETOL * LX))
                    @. mapPB[:,i] = mapMB[ids,j]

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
                    map!(x->x[1], ids, findall(@. D < NODETOL * LX))
                    @. mapPB[:,i] = mapMB[ids,j]

                    FToF[Bfaces[i]] = Bfaces[j]
                end
            end
        end
    end

    return mapPB[:] 
end