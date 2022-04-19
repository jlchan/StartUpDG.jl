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
                      face_vertex_indices::Dict{AbstractElemShape}) where {N}

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
    fnodes = fnodes[p, :]

    FToF = collect(1:NfacesTotal)
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

function build_node_maps(FToF, Xf...; tol = 1e-12)
    
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