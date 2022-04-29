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

    element_types = (keys(rds)..., ) # convert to tuple for indexing
                
    # TODO: this part assumes all faces have the same number of points (valid in 2D)
    rd = first(values(rds))
    num_points_per_face = rd.Nfq ÷ num_faces(rd.element_type)

    # TODO: fix, repeated code
    element_ids = Dict((Pair(elem, findall(length.(EToV) .== num_vertices(elem))) for elem in element_types))
    num_elements_of_type(elem) = length(element_ids[elem])
    
    # TODO: this part assumes all faces have the same number of points (valid in 2D)
    xf = hcat(reshape.(Xf[1].x, num_points_per_face, 
                       num_faces.(element_types) .* num_elements_of_type.(element_types))...)
    yf = hcat(reshape.(Xf[2].x, num_points_per_face, 
                       num_faces.(element_types) .* num_elements_of_type.(element_types))...)

    return build_node_maps(FToF, xf, yf)
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

# returns a Dict{element_type, RefElemData} when specifying multiple element types in 2D
function RefElemData(element_types::NTuple{N, <:Union{Tri, Quad}}, args...; kwargs...) where {N} 
    rds = Dict((elem => RefElemData(elem, args...; kwargs...) for elem in element_types))

    # check if number of face nodes 
    # TODO: this only works in 2D
    num_face_nodes = length.(getproperty.(values(rds), :rf)) .÷ num_faces.(element_types)
    allequal(x) = all(y->y==x[1],x)
    if !allequal(num_face_nodes)
        Base.@warn "Number of nodes per face for each element should be the same, but instead is:" num_face_nodes
    end
    return rds
end

# constructs MeshData for a hybrid mesh given a Dict of `RefElemData` 
# with element type keys (e.g., `Tri()` or `Quad`). 
function MeshData(VX, VY, EToV_unsorted, rds::Dict{AbstractElemShape, <:RefElemData};
                  is_periodic = (false, false))

    # sort EToV so that elements of the same type are contiguous
    # order by number of vertices, e.g., Quad(), then Tri()
    p = sortperm(length.(EToV_unsorted), order=Base.Order.Reverse)
    EToV = EToV_unsorted[p]

    # connect faces together 
    fvs = Dict(Pair(getproperty(rd, :element_type), getproperty(rd, :fv)) for rd in values(rds))
    FToF = StartUpDG.connect_mesh(EToV, fvs)

    # Dict between element type and element_ids of that type, e.g., element_ids[Tri()] = ...
    # We distinguish between different elements by the number of vertices. 
    # This should work in 3D too (might have issues if we ever do mixed 2D/3D meshes).
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
    # for each reference element in `rds`
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

    mapM, mapP, mapB = build_node_maps(rds, EToV, FToF, xyzf)

    return MeshData((VX, VY), EToV, FToF, 
                    xyz, xyzf, xyzq, wJq,
                    mapM, mapP, mapB,
                    rstxyzJ, J, nxyzJ, Jf, 
                    is_periodic)
end