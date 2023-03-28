# mesh_type identifier for hybrid meshes
struct HybridMesh{T, TV, TE}
    element_types::T
    VXYZ::TV
    EToV::TE
    # This is an inner constructor for HybridMesh. It sorts tuples of element types alphabetically 
    # so that `HybridMesh.element_types` always yields the same ordering of element types.
    function HybridMesh(element_types::T, VXYZ, EToV) where {T <: Tuple}
        p = sortperm(collect(typeof.(element_types)); by=string)
        return new{typeof(Tuple(element_types[p])), typeof(VXYZ), typeof(EToV)}(Tuple(element_types[p]), VXYZ, EToV)
    end
end

function Base.getproperty(x::MeshData{Dim, <:HybridMesh}, s::Symbol) where {Dim}

    if s===:VX
        return getfield(getfield(x, :mesh_type), :VXYZ)[1]
    elseif s===:VY
        return getfield(getfield(x, :mesh_type), :VXYZ)[2]
    elseif s===:VZ
        return getfield(getfield(x, :mesh_type), :VXYZ)[3]
    elseif s===:EToV
        return getfield(getfield(x, :mesh_type), s)
    else
        meshdata_getproperty(x, s)
    end
end


function HybridMeshExample()
    # Simple hybrid mesh for testing
    #   1  7______8______9
    #      |      | 3  / |
    #      |   4  |  / 5 |
    #   0  4 ---- 5 ---- 6 
    #      |      |      |
    #      |   1  |   2  |
    #   -1 1 ---- 2 ---- 3
    #     -1      0      1
    VX = [-1; 0; 1; -1; 0; 1; -1; 0; 1]
    VY = [-1; -1; -1; 0; 0; 0; 1; 1; 1]
    EToV = [[1 2 4 5], [2 3 5 6], [5 8 9], [4 5 7 8], [9 6 5]]
    return (VX, VY), EToV
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
                      face_vertex_indices::LittleDict{<:AbstractElemShape})

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

# returns element type of global element `global_e`
element_type(global_e, element_types, EToV) = 
    element_types[findfirst(length(EToV[global_e]) .== num_vertices.(element_types))]

# computes geometric terms from nodal coordinates
function compute_geometric_data(xyz, rd::RefElemData{2})
    x, y = xyz

    (; Dr, Ds ) = rd
    rxJ, sxJ, ryJ, syJ, J = geometric_factors(x, y, Dr, Ds)
    rstxyzJ = SMatrix{2, 2}(rxJ, ryJ, sxJ, syJ)

    (; Vq, Vf, wq ) = rd
    xyzf = map(x -> Vf * x, (x, y))
    xyzq = map(x -> Vq * x, (x, y))
    wJq = Diagonal(wq) * (Vq * J)

    nxJ, nyJ, Jf = compute_normals(rstxyzJ, rd.Vf, rd.nrstJ...)
    nxyzJ = (nxJ, nyJ)
    return (; xyzf, xyzq, wJq, rstxyzJ, J, nxyzJ, Jf)
end

# TODO: switch LittleDict to NamedTuple-based type for consistency (e.g., `rd.Quad.Dr` instead of `rd[Quad()].Dr`)
function RefElemData(element_types::NTuple{N, Union{Tri, Quad}}, args...; kwargs...) where {N} 
    # returns a LittleDict{element_type, RefElemData} when specifying multiple element types in 2D
    rds = LittleDict((elem => RefElemData(elem, args...; kwargs...) for elem in element_types)...)

    # check if number of face nodes is the same 
    num_face_nodes = length.(getproperty.(values(rds), :rf)) .÷ num_faces.(keys(rds))
    allequal(x) = all(y->y==x[1],x)
    if !allequal(num_face_nodes)
        Base.@warn "Number of nodes per face for each element should be the same, but instead is:" num_face_nodes
    end
    return rds
end

typename(x) = typeof(x).name.name

# constructs MeshData for a hybrid mesh given a LittleDict of `RefElemData` 
# with element type keys (e.g., `Tri()` or `Quad`). 
function MeshData(VX, VY, EToV_unsorted, rds::LittleDict{<:AbstractElemShape, <:RefElemData};
                  is_periodic = (false, false))

    # sort EToV so that elements of the same type are contiguous
    # order by number of vertices, e.g., Quad(), then Tri()
    p = sortperm(length.(EToV_unsorted))
    EToV = EToV_unsorted[p]

    # connect faces together 
    fvs = LittleDict((Pair(getproperty(rd, :element_type), 
                           getproperty(rd, :fv)) for rd in values(rds))...)
    FToF = StartUpDG.connect_mesh(EToV, fvs)

    # LittleDict between element type and element_ids of that type, e.g., element_ids[Tri()] = ...
    # We distinguish between different elements by the number of vertices. 
    # This should work in 3D too (but might have issues if we ever do mixed 2D/3D meshes).
    element_types = keys(rds)
    element_ids = LittleDict((Pair(elem, findall(length.(EToV) .== num_vertices(elem))) for elem in element_types))
    num_elements_of_type(elem) = length(element_ids[elem])

    # make node arrays 
    allocate_node_arrays(num_rows, element_type) = 
        ntuple(_ -> zeros(num_rows, num_elements_of_type(element_type)), ndims(element_type))
    xyz_hybrid = LittleDict((rd.element_type => allocate_node_arrays(size(rd.V1, 1), rd.element_type) for rd in values(rds)))
    for elem_type in element_types
        eids = element_ids[elem_type]
        x, y = xyz_hybrid[elem_type]
        (; V1 ) = rds[elem_type]
        for (e_local, e) in enumerate(eids)
            etov = EToV[e]        
            x[:, e_local] .= vec(V1 * VX[etov'])
            y[:, e_local] .= vec(V1 * VY[etov'])
        end
    end    

    # returns tuple of NamedTuples containing geometric fields 
    # for each reference element in `rds`
    geo = compute_geometric_data.(values(xyz_hybrid), values(rds))

    n_dims = 2
    xyz = ntuple(i -> NamedArrayPartition(NamedTuple(Pair.(typename.(keys(rds)), getindex.(values(xyz_hybrid), i)))), n_dims)
    xyzf = ntuple(i -> NamedArrayPartition(NamedTuple(Pair.(typename.(keys(rds)), getindex.(getproperty.(geo, :xyzf), i)))), n_dims)
    xyzq = ntuple(i -> NamedArrayPartition(NamedTuple(Pair.(typename.(keys(rds)), getindex.(getproperty.(geo, :xyzq), i)))), n_dims)

    # 4 entries in the geometric term matrix for 2D hybrid meshes
    rstxyzJ = ntuple(i -> NamedArrayPartition(NamedTuple(Pair.(typename.(keys(rds)), getindex.(getproperty.(geo, :rstxyzJ), i)))), n_dims * n_dims)
    rstxyzJ = SMatrix{2, 2}(rstxyzJ...)

    nxyzJ = ntuple(i -> NamedArrayPartition(NamedTuple(Pair.(typename.(keys(rds)), getindex.(getproperty.(geo, :nxyzJ), i)))), n_dims)
    wJq = NamedArrayPartition(NamedTuple(Pair.(typename.(keys(rds)), getproperty.(geo, :wJq))))
    J = NamedArrayPartition(NamedTuple(Pair.(typename.(keys(rds)), getproperty.(geo, :J))))
    Jf = NamedArrayPartition(NamedTuple(Pair.(typename.(keys(rds)), getproperty.(geo, :Jf))))

    mapM, mapP, mapB = vec.(build_node_maps(FToF, xyzf))

    return MeshData(HybridMesh(Tuple(element_types), (VX, VY), EToV), FToF, 
                    xyz, xyzf, xyzq, wJq,
                    mapM, mapP, mapB,
                    rstxyzJ, J, nxyzJ, Jf, 
                    is_periodic)
end

function MeshData(rds::LittleDict{<:AbstractElemShape, <:RefElemData{Dim}}, 
                  md::MeshData{Dim}, xyz_curved...) where {Dim}

    # TODO: can this be made type stable?
    tuple_fields = LittleDict{AbstractElemShape, NamedTuple}()
    scalar_fields = LittleDict{AbstractElemShape, NamedTuple}()
    for element_type in keys(rds)
        rd = rds[element_type]

        # compute curved geometric properties for each element type
        xyz = getproperty.(xyz_curved, typename(element_type))
        xyzf, xyzq, rstxyzJ, J, wJq, nxyzJ, Jf = recompute_geometry(rd, xyz)

        tuple_fields[element_type] = (; xyz, xyzq, xyzf, rstxyzJ, nxyzJ)
        scalar_fields[element_type] = (; J, wJq, Jf)
    end

    element_type_names = Symbol.(typename.(keys(rds)))

    xyz   = Tuple(NamedArrayPartition(NamedTuple(zip(element_type_names, (getproperty(tuple_fields[element_type], :xyz)[dim] for element_type in keys(rds))))) for dim in 1:Dim)
    xyzf  = Tuple(NamedArrayPartition(NamedTuple(zip(element_type_names, (getproperty(tuple_fields[element_type], :xyzf)[dim] for element_type in keys(rds))))) for dim in 1:Dim)
    xyzq  = Tuple(NamedArrayPartition(NamedTuple(zip(element_type_names, (getproperty(tuple_fields[element_type], :xyzq)[dim] for element_type in keys(rds))))) for dim in 1:Dim)
    nxyzJ = Tuple(NamedArrayPartition(NamedTuple(zip(element_type_names, (getproperty(tuple_fields[element_type], :nxyzJ)[dim] for element_type in keys(rds))))) for dim in 1:Dim)
    rstxyzJ = SMatrix{Dim, Dim}(NamedArrayPartition(NamedTuple(zip(element_type_names, (getproperty(tuple_fields[element_type], :rstxyzJ)[dim] for element_type in keys(rds))))) for dim in 1:Dim * Dim)

    J = NamedArrayPartition(NamedTuple(zip(element_type_names, (getproperty(scalar_fields[element_type], :J) for element_type in keys(rds)))))
    wJq = NamedArrayPartition(NamedTuple(zip(element_type_names, (getproperty(scalar_fields[element_type], :wJq) for element_type in keys(rds)))))
    Jf = NamedArrayPartition(NamedTuple(zip(element_type_names, (getproperty(scalar_fields[element_type], :Jf) for element_type in keys(rds)))))

    return setproperties(md, (; xyz, xyzq, xyzf, rstxyzJ, J, wJq, nxyzJ, Jf))
end
