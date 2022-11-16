# Volume Score functions -----------------------------------------------------------------
struct VolumeScore end

# Note: Scores should be constructed such that a larger score is considered a better match
function compute_nbhd_score(::VolumeScore, neighbor_list, md)
    total_volume = zero(eltype(md.wJq))
    for e in neighbor_list
        total_volume += sum(getcolumns(md.wJq, e))
    end

    return total_volume
end

function default_threshold_score(::VolumeScore, mesh_data)
    vx, vy = mesh_data.mesh_type.cut_cell_data.vxyz
    cartesian_cell_volume = maximum(diff(vx)) * maximum(diff(vy))
    return cartesian_cell_volume
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
            filter!((e)-> e == best_nbhr, neighbors)

            # Add the just-merged cell's neighbors to the list of neighbors
            new_neighbors = get_cartesian_nbhrs(best_nbhr, md)
            for nbhr in new_neighbors
                if nbhr in merge_nbhds[e]
                    push!(neighbors, nbhr)
                end
            end
            unique!(neighbors)

        end # while score is not good enough
    end # for each cut cell

    return merge_nbhds
end