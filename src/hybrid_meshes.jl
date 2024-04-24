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

"""
    struct MultipleRefElemData{T <: NamedTuple}
        data::T
    end
        
Holds multiple `RefElemData` objects in `data` where `typeof(data) <: NamedTuple`.

Individual `RefElemData` can be accessed via `getproperty`, for example `rds.Tri`. 
""" 
struct MultipleRefElemData{T <: NamedTuple}
    data::T
end

function Base.show(io::IO, ::MIME"text/plain", rds::MultipleRefElemData)
    @nospecialize rds
    print(io, "MultipleRefElemData: ")
    print(io, "\n")
    for rd in values(rds)
        print(io, "⋅ ")
        Base.show(io, MIME("text/plain"), rd)
        print(io, "\n")
    end
end

import Base: keys, values, getproperty
@inline keys(rds::MultipleRefElemData) = keys(getfield(rds, :data))
@inline values(rds::MultipleRefElemData) = values(getfield(rds, :data))
@inline getproperty(rds::MultipleRefElemData, s::Symbol) = getproperty(getfield(rds, :data), s)

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
                      face_vertex_indices, rds::MultipleRefElemData)

    elem_types = Tuple(getproperty(rds, key).element_type for key in keys(face_vertex_indices))

    # EToV = vector of index vectors
    K = length(EToV)    

    # create fnodes
    fnodes = Vector{eltype(first(EToV))}[]
    for e in 1:K
        vertex_ids = EToV[e]
        element_type = element_type_from_num_vertices(elem_types, length(vertex_ids))
        for ids in getproperty(face_vertex_indices, typename(element_type))
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

typename(x) = typeof(x).name.name

function RefElemData(element_types::NTuple{NT, Union{Tri, Quad}}, N::Int; 
                     quad_rule_face=gauss_quad(0, 0, N), kwargs...) where {NT} 

    rds = NamedTuple((typename(elem) => 
        RefElemData(elem, N; quad_rule_face, kwargs...) for elem in element_types))

    # check if number of face nodes is the same 
    num_face_nodes = length.(getproperty.(values(rds), :rf)) .÷ num_faces.(getproperty.(values(rds), :element_type))
    allequal(x) = all(y -> y == x[1], x)
    if !allequal(num_face_nodes)
        Base.@warn "Number of nodes per face for each element should be the same, but instead is:" num_face_nodes
    end

    return MultipleRefElemData(rds)
end

MeshData(mesh::Tuple{<:Tuple, Vector{Matrix{Int64}}}, rd::MultipleRefElemData, other_args...) = 
    MeshData(mesh..., rd; other_args...)
MeshData(VXYZ::Tuple, EToV, rd::MultipleRefElemData, other_args...) = 
    MeshData(VXYZ..., EToV, rd; other_args...) # splats VXYZ 

# constructs MeshData for a hybrid mesh given a NamedTuple of `RefElemData` 
# with element type keys (e.g., `:Tri` or `:Quad`). 
function MeshData(VX, VY, EToV_unsorted, 
                  rds::MultipleRefElemData;
                  is_periodic = (false, false))

    # sort EToV so that elements of the same type are contiguous
    # order by number of vertices, e.g., Quad(), then Tri()
    p = sortperm(length.(EToV_unsorted))
    EToV = EToV_unsorted[p]

    # connect faces together 
    face_vertex_indices = NamedTuple((Pair(typename(getproperty(rd, :element_type)), 
                                           getproperty(rd, :fv)) for rd in values(rds)))

    FToF = StartUpDG.connect_mesh(EToV, face_vertex_indices, rds)

    # We distinguish between different elements by the number of vertices. 
    # This should work in 3D too (but might have issues if we ever do mixed 2D/3D meshes).
    element_types = Tuple(getproperty(rd, :element_type) for rd in values(rds))
    element_ids = NamedTuple((Pair(typename(elem), findall(length.(EToV) .== num_vertices(elem))) for elem in element_types))
    num_elements_of_type(elem) = length(getproperty(element_ids, typename(elem)))

    # make node arrays 
    allocate_node_arrays(num_rows, element_type) = 
        ntuple(_ -> zeros(num_rows, num_elements_of_type(element_type)), ndims(element_type))

    xyz_hybrid = NamedTuple((typename(rd.element_type) => allocate_node_arrays(size(rd.V1, 1), rd.element_type) for rd in values(rds)))
    for elem_type in element_types
        eids = getproperty(element_ids, typename(elem_type))
        x, y = getproperty(xyz_hybrid, typename(elem_type))

        (; V1 ) = getproperty(rds, typename(elem_type))
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
    xyz = ntuple(i -> NamedArrayPartition(NamedTuple(Pair.(keys(rds), getindex.(values(xyz_hybrid), i)))), n_dims)
    xyzf = ntuple(i -> NamedArrayPartition(NamedTuple(Pair.(keys(rds), getindex.(getproperty.(geo, :xyzf), i)))), n_dims)
    xyzq = ntuple(i -> NamedArrayPartition(NamedTuple(Pair.(keys(rds), getindex.(getproperty.(geo, :xyzq), i)))), n_dims)

    # 4 entries in the geometric term matrix for 2D hybrid meshes
    rstxyzJ = ntuple(i -> NamedArrayPartition(NamedTuple(Pair.(keys(rds), getindex.(getproperty.(geo, :rstxyzJ), i)))), n_dims * n_dims)
    rstxyzJ = SMatrix{2, 2}(rstxyzJ...)

    nxyzJ = ntuple(i -> NamedArrayPartition(NamedTuple(Pair.(keys(rds), getindex.(getproperty.(geo, :nxyzJ), i)))), n_dims)
    wJq = NamedArrayPartition(NamedTuple(Pair.(keys(rds), getproperty.(geo, :wJq))))
    J = NamedArrayPartition(NamedTuple(Pair.(keys(rds), getproperty.(geo, :J))))
    Jf = NamedArrayPartition(NamedTuple(Pair.(keys(rds), getproperty.(geo, :Jf))))

    mapM, mapP, mapB = vec.(build_node_maps(FToF, xyzf))

    periodicity = (false, false)
    md = MeshData(HybridMesh(Tuple(element_types), (VX, VY), EToV), FToF, 
                  xyz, xyzf, xyzq, wJq,
                  mapM, mapP, mapB,
                  rstxyzJ, J, nxyzJ, Jf, 
                  periodicity)

    if any(is_periodic)     
        @warn "Periodic boundary conditions not yet implemented for hybrid meshes."
        # md = make_periodic(md, is_periodic)
    end

    return md
end

function MeshData(rds::MultipleRefElemData, 
                  md::MeshData{Dim}, xyz_curved...) where {Dim}

    # TODO: can this be made type stable?
    tuple_fields = Dict{AbstractElemShape, NamedTuple}()
    scalar_fields = Dict{AbstractElemShape, NamedTuple}()
    for rd in values(rds)
        # compute curved geometric properties for each element type
        xyz = getproperty.(xyz_curved, typename(rd.element_type))
        xyzf, xyzq, rstxyzJ, J, wJq, nxyzJ, Jf = recompute_geometry(rd, xyz)

        tuple_fields[rd.element_type] = (; xyz, xyzq, xyzf, rstxyzJ, nxyzJ)
        scalar_fields[rd.element_type] = (; J, wJq, Jf)
    end

    xyz   = Tuple(NamedArrayPartition(NamedTuple(zip(keys(rds), (getproperty(tuple_fields[rd.element_type], :xyz)[dim] for rd in values(rds))))) for dim in 1:Dim)
    xyzf  = Tuple(NamedArrayPartition(NamedTuple(zip(keys(rds), (getproperty(tuple_fields[rd.element_type], :xyzf)[dim] for rd in values(rds))))) for dim in 1:Dim)
    xyzq  = Tuple(NamedArrayPartition(NamedTuple(zip(keys(rds), (getproperty(tuple_fields[rd.element_type], :xyzq)[dim] for rd in values(rds))))) for dim in 1:Dim)
    nxyzJ = Tuple(NamedArrayPartition(NamedTuple(zip(keys(rds), (getproperty(tuple_fields[rd.element_type], :nxyzJ)[dim] for rd in values(rds))))) for dim in 1:Dim)
    rstxyzJ = SMatrix{Dim, Dim}(NamedArrayPartition(NamedTuple(zip(keys(rds), (getproperty(tuple_fields[rd.element_type], :rstxyzJ)[dim] for rd in values(rds))))) for dim in 1:Dim * Dim)

    J = NamedArrayPartition(NamedTuple(zip(keys(rds), (getproperty(scalar_fields[rd.element_type], :J) for rd in values(rds)))))
    wJq = NamedArrayPartition(NamedTuple(zip(keys(rds), (getproperty(scalar_fields[rd.element_type], :wJq) for rd in values(rds)))))
    Jf = NamedArrayPartition(NamedTuple(zip(keys(rds), (getproperty(scalar_fields[rd.element_type], :Jf) for rd in values(rds)))))

    return setproperties(md, (; xyz, xyzq, xyzf, rstxyzJ, J, wJq, nxyzJ, Jf))
end
