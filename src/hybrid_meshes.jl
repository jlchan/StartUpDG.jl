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

# construct node connectivity arrays for hybrid meshes. 
# note that `rds` (the container of `RefElemData`s) must have the 
# same ordering as the `ArrayPartition` `Xf`.
function build_node_maps(rds::Dict{AbstractElemShape, <:RefElemData}, FToF, 
                         Xf::NTuple{2, ArrayPartition}; tol = 1e-12)
    xf, yf = Xf
    for rd in values(rds)
        Nfq = size(rd.Vf, 1)
        
    end

    # TODO: finish
    return nothing, nothing, nothing    

    # NfacesK = length(FToF)
    # dims = length(Xf)
    # # number nodes consecutively
    # Nfp  = length(Xf[1]) ÷ NfacesK
    # mapM = reshape(collect(1:length(Xf[1])), Nfp, NfacesK)
    # mapP = copy(mapM)
    # D = zeros(Nfp, Nfp)
    # idM, idP = zeros(Int, Nfp), zeros(Int, Nfp)
    # for (f1, f2) in enumerate(FToF)
    #     fill!(D, zero(eltype(D)))
    #     # find volume node numbers of left and right nodes
    #     for i in 1:dims
    #         Xfi = reshape(Xf[i], Nfp, NfacesK)
    #         for j in 1:Nfp, k in 1:Nfp
    #             D[j, k] += abs(Xfi[j, f1] - Xfi[k, f2])
    #         end
    #     end
    #     refd = maximum(D[:])
    #     map!(id -> id[1], idM, findall(@. D < tol * refd))
    #     map!(id -> id[2], idP, findall(@. D < tol * refd))        
    #     @. mapP[idM, f1] = idP + (f2 - 1) * Nfp
    # end
    # mapB = map(x -> x[1], findall(@. mapM[:]==mapP[:]))
    # return mapM, mapP, mapB
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

    fvs = Dict(Pair(getproperty(rd, :element_type), getproperty(rd, :fv)) for rd in values(rds))
    FToF = StartUpDG.connect_mesh(EToV, fvs)

    # Dict between element type and element_ids of that type, e.g., element_ids[Tri()] = ...
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

    mapM, mapP, mapB = build_node_maps(rds, FToF, xyzf)

    return MeshData((VX, VY), EToV, FToF, 
                    xyz, xyzf, xyzq, wJq,
                    mapM, mapP, mapB,
                    rstxyzJ, J, nxyzJ, Jf, 
                    is_periodic)
end