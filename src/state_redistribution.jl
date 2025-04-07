# create a typed index, e.g., CellIndex{Cut}(e)
struct CellIndex{CellT}
    index::Int
end

# convenience extractors for `MeshData{2, CutCellMesh}` fields, which have `cartesian` and `cut` properties. 
struct Cut end
struct Cartesian end
_propertyname(::CellIndex{Cut}) = :cut
_propertyname(::CellIndex{Cartesian}) = :cartesian

getcolumns(x, i::CellIndex) = view(getproperty(x, _propertyname(i)), :, i.index)
getcolumns(x, indices::AbstractVector{<:CellIndex}) = (getcolumns(x, i.index) for i in indices)

vcat_columns(x, list::AbstractVector{<:CellIndex}) = vcat((vec(getcolumns(x, i)) for i in list)...)

get_face_nodes(x, i::CellIndex{Cartesian}, args...) = view(x.cartesian, :, i.index)
get_face_nodes(x, i::CellIndex{Cut}, md::MeshData) = view(x.cut, md.mesh_type.cut_face_nodes[i.index])

# ================== neighborhood computation code ================
struct VolumeScore end

# Note: Scores should be constructed such that a larger score is considered a better match
function compute_nbhd_score(::VolumeScore, neighbor_list, md)
    total_volume = zero(eltype(md.wJq))
    for e in neighbor_list
        total_volume += sum(getcolumns(md.wJq, e))
    end

    return total_volume
end

# merge until a cell is some ratio of the volume of a uniform Cartesian cell.
# the default ratio is 1/2. 
function default_threshold_score(::VolumeScore, mesh_data; volume_ratio = 0.5)
    vx, vy = mesh_data.mesh_type.cut_cell_data.vxyz
    cartesian_cell_volume = maximum(diff(vx)) * maximum(diff(vy))
    return volume_ratio * cartesian_cell_volume
end

# Helper functions for constructing merge neighborhoods ---------------------------------

function create_cell_index(e, flag)
    if is_Cartesian(flag)
        return CellIndex{Cartesian}(e)
    elseif is_cut(flag)
        return CellIndex{Cut}(e)
    end
end

function get_cartesian_nbhrs(e::CellIndex{T}, md) where T
    if T == Cut
        ex, ey = md.mesh_type.cut_cell_data.linear_to_cartesian_element_indices.cut[e.index];
    elseif T == Cartesian
        ex, ey = md.mesh_type.cut_cell_data.linear_to_cartesian_element_indices.cartesian[e.index];
    end
    neighbors = CellIndex[]

    if is_inside_domain(ex+1, ey, md.mesh_type.cut_cell_data.region_flags)
        e_nbhr = md.mesh_type.cut_cell_data.cartesian_to_linear_element_indices[ex+1,ey]
        flag = md.mesh_type.cut_cell_data.region_flags[ex+1, ey]
        push!(neighbors, create_cell_index(e_nbhr, flag))
    end
    if is_inside_domain(ex-1, ey, md.mesh_type.cut_cell_data.region_flags)
        e_nbhr = md.mesh_type.cut_cell_data.cartesian_to_linear_element_indices[ex-1,ey]
        flag = md.mesh_type.cut_cell_data.region_flags[ex-1, ey]
        push!(neighbors, create_cell_index(e_nbhr, flag))
    end
    if is_inside_domain(ex, ey+1, md.mesh_type.cut_cell_data.region_flags)
        e_nbhr = md.mesh_type.cut_cell_data.cartesian_to_linear_element_indices[ex,ey+1]
        flag = md.mesh_type.cut_cell_data.region_flags[ex, ey+1]
        push!(neighbors, create_cell_index(e_nbhr, flag))
    end
    if is_inside_domain(ex, ey-1, md.mesh_type.cut_cell_data.region_flags)
        e_nbhr = md.mesh_type.cut_cell_data.cartesian_to_linear_element_indices[ex,ey-1]
        flag = md.mesh_type.cut_cell_data.region_flags[ex, ey-1]
        push!(neighbors, create_cell_index(e_nbhr, flag))
    end

    return neighbors
end


function compute_neighbor_list(md; threshold_score=-Inf, score_type=VolumeScore(), score_params=md)

    if threshold_score == -Inf
        threshold_score = default_threshold_score(score_type, score_params)
    end

    merge_nbhds = Vector{CellIndex}[[CellIndex{Cut}(e)] for e in 1:num_cut_elements(md)]
    for e in 1:num_cut_elements(md)
        # Initialize the merge nbhd with the its featured cut cell
        merged_score = compute_nbhd_score(score_type, merge_nbhds[e], score_params)

        # Find the neighbors of the initial cut cell
        neighbors = get_cartesian_nbhrs(CellIndex{Cut}(e), md)
        
        while length(neighbors) > 0 && merged_score < threshold_score
            best_nbhr = neighbors[1]
            best_score = compute_nbhd_score(score_type, vcat(merge_nbhds[e], neighbors[1]), score_params)
            
            for e_nbhr in neighbors[2:end]
                new_score = compute_nbhd_score(score_type, vcat(merge_nbhds[e], e_nbhr), score_params)
                if new_score >= best_score
                    best_score = new_score
                    best_nbhr = e_nbhr
                end
            end

            # Update the merge nbhd
            push!(merge_nbhds[e], best_nbhr)
            merged_score = best_score

            # Remove the just-merged element from the list of neighbors
            filter!((e)-> e != best_nbhr, neighbors)

            # Add the just-merged cell's neighbors to the list of neighbors
            new_neighbors = get_cartesian_nbhrs(best_nbhr, md)
            for nbhr in new_neighbors
                if !(nbhr in merge_nbhds[e])
                    push!(neighbors, nbhr)
                end
            end
            unique!(neighbors)

        end # while score is not good enough
    end # for each cut cell

    return merge_nbhds
end

# end neighborhood computation code ================

struct StateRedistribution{TP, TN, TE, TO, TU}
    projection_operators::TP
    projection_indices::TN
    projection_indices_by_element::TE
    overlap_counts::TO
    u_tmp::TU # temporary storage for operations
end

"""
    function StateRedistribution(rd::RefElemData{2, Quad}, 
        md::MeshData{2, <:CutCellMesh}; solution_eltype=Float64)

Constructs a state redistribution instance `srd`, which can be applied to a DG solution via
`apply!(u, srd::StateRedistribution)` or `srd(u)`. 

Note that, by default, the constructor assumes that the solution `u` has eltype `Float64`; if 
this is not the case, then `solution_eltype` should be specified in the constructor.
"""
function StateRedistribution(rd::RefElemData{2, Quad}, md::MeshData{2, <:CutCellMesh}; solution_eltype=Float64)
    (; physical_frame_elements) = md.mesh_type

    neighbor_list = compute_neighbor_list(md)    

    # indexing by elements is a little tricky. for consistency, we store overlap counts
    # separately for cut and cartesian cells. 
    overlap_counts = NamedArrayPartition(cartesian=ones(Int, num_cartesian_elements(md)), 
                                         cut=zeros(Int, num_cut_elements(md)))
    for neighbors in neighbor_list
        for e in neighbors
            # equivalent to `overlap_counts.cut[e] +=1` (similarly for `cartesian`)
            getproperty(overlap_counts, _propertyname(e))[e.index] += 1
        end
    end

    # cartesian first, then cut     
    cut_indices = (1:length(md.x.cut)) .+ length(md.x.cartesian)
    indices = NamedArrayPartition(cartesian=reshape(1:length(md.x.cartesian), size(md.x.cartesian)),
                                  cut=reshape(cut_indices, size(md.x.cut)))

    # scale weights by overlap counts
    wJq = copy(md.wJq)
    wJq.cut *= Diagonal(1 ./ overlap_counts.cut)
    wJq.cartesian *= Diagonal(1 ./ overlap_counts.cartesian)

    # compute projection operators
    projection_operators = Matrix{eltype(wJq)}[]
    projection_indices_by_element = Vector{Vector{eltype(indices)}}[]
    projection_indices = Vector{eltype(indices)}[]
    for neighbors in neighbor_list

        xf = vcat((get_face_nodes(md.xf, e, md) for e in neighbors)...)
        yf = vcat((get_face_nodes(md.yf, e, md) for e in neighbors)...)
        merged_elem = PhysicalFrame(xf, yf)

        # "merged" Vandermonde matrix
        wq = vcat_columns(wJq, neighbors)
        Vq = vandermonde(merged_elem, rd.N, vcat_columns(md.xq, neighbors), vcat_columns(md.yq, neighbors))
        M = Vq' * Diagonal(wq) * Vq

        # individual cell Vandermonde matrices
        Vq_list = typeof(Vq)[]
        for e in neighbors
            if _propertyname(e) == :cut
                VDM = vandermonde(physical_frame_elements[e.index], rd.N, 
                                view(md.x.cut, :, e.index), view(md.y.cut, :, e.index))
                Vq_e = vandermonde(physical_frame_elements[e.index], rd.N, 
                                view(md.xq.cut, :, e.index), view(md.yq.cut, :, e.index))                                  
                push!(Vq_list, Vq_e / VDM)
            else # vandermonde matrix for a nodal basis
                Vq_cartesian = vandermonde(rd.element_type, rd.N, rd.rstq...)
                push!(Vq_list, Vq_cartesian / rd.VDM)
            end
        end

        # blockdiag is only defined for sparse matrices
        B = Vq' * Diagonal(wq) * blockdiag(sparse.(Vq_list)...)

        eval_at_nodal_points = 
            vandermonde(merged_elem, rd.N, vcat_columns(md.x, neighbors), vcat_columns(md.y, neighbors))
        projection_operator = eval_at_nodal_points * (M \ B)

        push!(projection_operators, projection_operator)        
        push!(projection_indices_by_element, [getcolumns(indices, e) for e in neighbors])
        push!(projection_indices, vcat_columns(indices, neighbors))
    end 

    # temporary storage for state redistribution - stores output of projection operators
    u_tmp = [zeros(solution_eltype, size(projection_operators[e], 2)) 
        for e in eachindex(projection_operators)]

    return StateRedistribution(projection_operators, projection_indices, 
                               projection_indices_by_element, overlap_counts, u_tmp)
end

# functor syntax, e.g., `srd(u)` applies state redistribution
(srd::StateRedistribution)(u) = apply!(u, srd::StateRedistribution)

function apply!(u, srd::StateRedistribution)
    (; projection_operators, overlap_counts, u_tmp) = srd
    (; projection_indices, projection_indices_by_element) = srd

    # cell merging
    for i in eachindex(projection_indices, projection_operators)
        mul!(u_tmp[i], projection_operators[i], view(u, projection_indices[i]))
    end    

    # zero out cut cell solutions (these will be overwritten by projections)
    fill!(u.cut, zero(eltype(u)))

    # sum up overlapping solutions
    for i in eachindex(projection_indices, projection_operators)
        offset = zero(eltype(first(projection_indices)))
        for (id, element_indices) in enumerate(projection_indices_by_element[i])
            local_ids = (1:length(projection_indices_by_element[i][id])) .+ offset
            view(u, element_indices) .+= u_tmp[i][local_ids] 
            offset += length(projection_indices_by_element[i][id])
        end
    end

    # averaging of summed solutions
    for (e, overlap) in enumerate(overlap_counts.cut)
        view(u.cut, :, e) ./= overlap
    end
    for (e, overlap) in enumerate(overlap_counts.cartesian)
        view(u.cartesian, :, e) ./= overlap
    end
end