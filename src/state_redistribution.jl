# convenience extractors for `MeshData{2, CutCellMesh}` fields, which are `ComponentArrays`. 
struct Cut end
struct Cartesian end

# create a typed index, e.g., CellIndex{Cut}(e)
struct CellIndex{CellT}
    index::Int
end

import Base: getindex
property_name(::CellIndex{Cut}) = :cut
property_name(::CellIndex{Cartesian}) = :cartesian
getindex(x::ComponentArray, i::CellIndex) = getindex(getproperty(x, property_name(i)), i.index)
getcolumns(x::ComponentArray, i::CellIndex) = view(getproperty(x, property_name(i)), :, i.index)

get_face_nodes(x::ComponentArray, i::CellIndex{Cartesian}, args...) = view(x.cartesian, :, i.index)
get_face_nodes(x::ComponentArray, i::CellIndex{Cut}, md::MeshData) = view(x.cut, md.mesh_type.cut_face_nodes[i.index])

getcolumns(x::ComponentArray, indices::AbstractVector{<:CellIndex}) = (getcolumns(x, i.index) for i in indices)
vcat_columns(x::ComponentArray, list::AbstractVector{<:CellIndex}) = vcat((vec(getcolumns(x, i)) for i in list)...)

struct StateRedistribution{TP, TN, TO, TU}
    projection_operators::TP
    projection_indices::TN
    overlap_counts::TO
    u_tmp::TU # temporary storage for operations
end

function StateRedistribution(rd::RefElemData{2, Quad}, md::MeshData{2, <:CutCellMesh})
    (; physical_frame_elements) = md.mesh_type
    # (; cut_cell_data) = md.mesh_type
    # (; region_flags, cells_per_dimension, element_indices) = cut_cell_data

    # compute neighbor list
    neighbor_list = Vector{CellIndex}[[CellIndex{Cut}(e)] for e in 1:num_cut_elements(md)]
    push!(neighbor_list[4], CellIndex{Cut}(1))
    push!(neighbor_list[4], CellIndex{Cartesian}(1))

    # indexing by elements is a little tricky. for consistency, we store overlap counts
    # separately for cut and cartesian cells. 
    overlap_counts = ComponentArray(cartesian=ones(Int, num_cartesian_elements(md)), 
                                    cut=zeros(Int, num_cut_elements(md)))
    for neighbors in neighbor_list
        for e in neighbors
            # equivalent to `overlap_counts.cut[e] +=1` (similarly for `cartesian`)
            getproperty(overlap_counts, property_name(e))[e.index] += 1
        end
    end

    # cartesian first, then cut     
    cut_indices = (1:length(md.x.cut)) .+ length(md.x.cartesian)
    indices = ComponentArray(cartesian=reshape(1:length(md.x.cartesian), size(md.x.cartesian)),
                             cut=reshape(cut_indices, size(md.x.cut)))

    # scale weights by overlap counts
    wJq = copy(md.wJq)
    wJq.cut *= Diagonal(1 ./ overlap_counts.cut)
    wJq.cartesian *= Diagonal(1 ./ overlap_counts.cartesian)

    # compute projection operators
    projection_operators = Matrix{eltype(wJq)}[]
    projection_indices = Vector{eltype(indices)}[]
    for neighbors in neighbor_list

        wq = vcat_columns(wJq, neighbors)
        xq, yq = vcat_columns(md.xq, neighbors), vcat_columns(md.yq, neighbors)        

        xf = vcat((get_face_nodes(md.xf, e, md) for e in neighbors)...)
        yf = vcat((get_face_nodes(md.yf, e, md) for e in neighbors)...)
        merged_elem = PhysicalFrame(xf, yf)

        # "merged" Vandermonde matrix
        Vq = vandermonde(merged_elem, rd.N, xq, yq)
        M = Vq' * Diagonal(vcat(wq...)) * Vq

        # individual cell Vandermonde matrices
        Vq_list = typeof(Vq)[]
        for e in neighbors
            if property_name(e) == :cut
                VDM = vandermonde(physical_frame_elements[e.index], rd.N, 
                                view(md.x.cut, :, e.index), view(md.y.cut, :, e.index))
                Vq_e = vandermonde(physical_frame_elements[e.index], rd.N, 
                                view(md.xq.cut, :, e.index), view(md.yq.cut, :, e.index))                                  
                push!(Vq_list, Vq_e / VDM)
            else # vandermonde matrix for a nodal basis is identity
                push!(Vq_list, Matrix(I(rd.Np)))
            end
        end

        # blockdiag is only defined for sparse matrices
        B = Vq' * Diagonal(wq) * blockdiag(sparse.(Vq_list)...)

        interp_to_nodal_points = 
            vandermonde(merged_elem, rd.N, vcat_columns(md.x, neighbors), vcat_columns(md.y, neighbors))
        projection_operator = interp_to_nodal_points * (M \ B)

        push!(projection_operators, projection_operator)        
        push!(projection_indices, vcat_columns(indices, neighbors))
    end 

    # temporary storage for state redistribution - stores output of projection operators
    u_tmp = [zeros(size(projection_operators[e], 2)) for e in eachindex(projection_operators)]

    return StateRedistribution(projection_operators, projection_indices, overlap_counts, u_tmp)
end

function (srd::StateRedistribution)(u)
    (; projection_operators, projection_indices, overlap_counts, u_tmp)= srd
    for i in eachindex(projection_indices, projection_operators)
        mul!(u_tmp[i], projection_operators[i], view(u, projection_indices[i]))
    end
end