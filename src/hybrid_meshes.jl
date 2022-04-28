import Base: *
function (*)(A::ArrayPartition, B::ArrayPartition)
    if all(size.(A.x, 2) .!== size.(B.x, 1))
        throw(DimensionMismatch("size.(A, 2), $(size.(A.x, 2)), does not match size.(B, 1), $(size.(B.x, 1))"))
    end
    C = ArrayPartition(zeros.(eltype(A), size.(A.x, 1), size.(B.x, 2)))
    return LinearAlgebra.mul!(C, A, B)
end

# Given a tuple of element types, find which element type has 
# `num_vertices_of_target` vertices. 
function element_type_from_num_vertices(elem_types, num_vertices_of_target)
    for elem_type in elem_types
        if num_vertices(elem_type) == num_vertices_of_target
            return elem_type
        end
    end
end

# if EToV is an array of arrays, treat it as a "ragged" index array for a hybrid mesh.
function connect_mesh(EToV::AbstractVector{<:AbstractArray}, 
                      face_vertex_indices::Dict{AbstractElemShape})

    elem_types = (keys(face_vertex_indices)...,)

    # EToV = vector of index vectors
    K = length(EToV)    

    # create fnodes
    fnodes = Vector{eltype(first(EToV))}[]
    for e in 1:K
        vertex_ids = EToV[e]
        element_type = element_type_from_num_vertices(elem_types, length(vertex_ids))
        for ids in face_vertex_indices[element_type]
            push!(fnodes, sort(EToV[e][ids]))
        end
    end

    num_vertices_per_elem = length.(EToV)
    Nfaces_per_elem = [num_faces(element_type_from_num_vertices(elem_types, nv)) 
                       for nv in num_vertices_per_elem]
    NfacesTotal = sum(Nfaces_per_elem)

    # sort and find matches
    p = sortperm(fnodes) # sorts by lexicographic ordering by default
    fnodes = fnodes[p]

    FToF = collect(1:NfacesTotal)
    for f = 1:size(fnodes, 1) - 1
        if fnodes[f]==fnodes[f + 1]
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
function match_coordinate_vectors(u, v; tol = 1e-12)
    p = zeros(Int, length(first(u)))
    for i in eachindex(first(u)), j in eachindex((first(v)))
        if norm(getindex.(u, i) .- getindex.(v, j)) < tol
            p[i] = j
        end
    end
    return p # [findfirst(abs.(u[i] .- v) .< tol) for i in eachindex(u)]
end

# returns element type of global element `global_e`
element_type(global_e, element_types, EToV) = 
    element_types[findfirst(length(EToV[global_e]) .== num_vertices.(element_types))]

# construct node connectivity arrays for hybrid meshes. 
# note that `rds` (the container of `RefElemData`s) must have the 
# same ordering as the `ArrayPartition` `Xf`.
# Here `element_ids` is a Dict{AbstractElemShape, "indices of elements"}.
function build_node_maps(rds::Dict{AbstractElemShape, <:RefElemData{2}}, 
                         EToV, FToF, Xf::NTuple{2, ArrayPartition}; tol = 1e-12)

    # TODO: fix, this is repeated code                     
    element_types = (keys(rds)..., ) # convert to tuple for indexing
    element_ids = Dict((Pair(elem, findall(length.(EToV) .== num_vertices(elem))) for elem in element_types))
                     
    # NOTE: this part assumes all faces have the same number of points (valid in 2D)
    rd = first(values(rds))
    num_points_per_face = rd.Nfq ÷ num_faces(rd.element_type)
    fids(f) = (1:num_points_per_face) .+ (f-1) * num_points_per_face    
    xf, yf = ntuple(_ -> [zeros(num_points_per_face) for _ in 1:length(FToF)], 2)

    # create list of element types
    global_faces = UnitRange{Int}[]
    face_offset = 0
    for e in 1:length(EToV)        
        elem_type = element_type(e, element_types, EToV)        
        push!(global_faces, (1:num_faces(elem_type)) .+ face_offset)
        face_offset += num_faces(elem_type)
    end

    # create xf, yf = vector of vectors::{points on each face}
    for (elem_type_id, rd) in enumerate(values(rds))
        elem_type = rd.element_type
        for (e_local, e) in enumerate(element_ids[elem_type])
            for f in 1:num_faces(elem_type)                
                xf[global_faces[e][f]] .= Xf[1].x[elem_type_id][fids(f), e_local]
                yf[global_faces[e][f]] .= Xf[2].x[elem_type_id][fids(f), e_local]
            end
        end
    end

    mapM = [collect(fids(f)) for f in 1:length(FToF)]
    mapP = copy(mapM)
    for (f1, f2) in enumerate(FToF)
        if f1 != f2
            x1, y1 = xf[f1], yf[f1]
            x2, y2 = xf[f2], yf[f2]
            p = match_coordinate_vectors((x1, y1), (x2, y2))
            mapP[f1] .= mapM[f2][p]
        end
    end

    mapB = nothing

    # TODO: finish
    return mapM, mapP, mapB
end

# computes geometric terms from nodal coordinates
function compute_geometric_data(xyz, rd::RefElemData{2})
    x, y = xyz

    @unpack Dr, Ds = rd
    rxJ, sxJ, ryJ, syJ, J = geometric_factors(x, y, Dr, Ds)
    rstxyzJ = SMatrix{2, 2}(rxJ, ryJ, sxJ, syJ)

    @unpack Vq, Vf, wq = rd
    xyzf = map(x -> Vf * x, (x, y))
    xyzq = map(x -> Vq * x, (x, y))
    wJq = Diagonal(wq) * (Vq * J)

    nxJ, nyJ, Jf = compute_normals(rstxyzJ, rd.Vf, rd.nrstJ...)
    nxyzJ = (nxJ, nyJ)
    return (; xyzf, xyzq, wJq, rstxyzJ, J, nxyzJ, Jf)
end

# constructs MeshData for a hybrid mesh given a Dict of `RefElemData` 
# with element type keys (e.g., `Tri()` or `Quad`). 
function MeshData(VX, VY, EToV, rds::Dict{AbstractElemShape, <:RefElemData};
                  is_periodic = (false, false))

    # # sort EToV so that elements of the same type are contiguous
    # p = sortperm(length.(EToV_unsorted))
    # EToV = EToV_unsorted[p]

    fvs = Dict(Pair(getproperty(rd, :element_type), getproperty(rd, :fv)) for rd in values(rds))
    FToF = StartUpDG.connect_mesh(EToV, fvs)

    # Dict between element type and element_ids of that type, e.g., element_ids[Tri()] = ...
    # We distinguish between different elements by the number of vertices. 
    # This should work in 3D too, might have issues for 2D/3D though.
    element_types = keys(rds)
    element_ids = Dict((Pair(elem, findall(length.(EToV) .== num_vertices(elem))) for elem in element_types))
    num_elements_of_type(elem) = length(element_ids[elem])

    # make node arrays 
    allocate_node_arrays(num_rows, element_type) = ntuple(_ -> zeros(num_rows, num_elements_of_type(element_type)), 
                                                          ndims(element_type))
    xyz_hybrid = Dict((rd.element_type => allocate_node_arrays(size(rd.V1, 1), rd.element_type) for rd in values(rds)))
    for elem_type in element_types
        eids = element_ids[elem_type]
        x, y = xyz_hybrid[elem_type]
        @unpack V1 = rds[elem_type]
        for (e_local, e) in enumerate(eids)
            etov = EToV[e]        
            x[:, e_local] .= V1 * VX[etov']
            y[:, e_local] .= V1 * VY[etov']
        end
    end

    # returns tuple of NamedTuples containing geometric fields
    geo = compute_geometric_data.(values(xyz_hybrid), values(rds))

    # create array partitions for all geometric quantities
    xyz = ArrayPartition.(values(xyz_hybrid)...)    
    xyzf = ArrayPartition.(getproperty.(geo, :xyzf)...)
    xyzq = ArrayPartition.(getproperty.(geo, :xyzq)...)
    wJq = ArrayPartition(getproperty.(geo, :wJq)...)
    rstxyzJ = ArrayPartition.(getproperty.(geo, :rstxyzJ)...)
    J = ArrayPartition(getproperty.(geo, :J)...)
    nxyzJ = ArrayPartition.(getproperty.(geo, :nxyzJ)...)
    Jf = ArrayPartition(getproperty.(geo, :Jf)...)    

    # TODO fix this
    mapM, mapP, mapB = build_node_maps(rds, EToV, FToF, xyzf)

    return MeshData((VX, VY), EToV, FToF, 
                    xyz, xyzf, xyzq, wJq,
                    mapM, mapP, mapB,
                    rstxyzJ, J, nxyzJ, Jf, 
                    is_periodic)
end