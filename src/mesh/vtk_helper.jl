
"""
    n_verts_between(n, from, to, dim)

Compute the coordinates of n equally distributed points between the points
given by `from` and `to`. `dim` is the dimension of `from` and `to`. 
Inspired by: https://github.com/ju-kreber/paraview-scripts/blob/master/node_ordering.py
"""
function n_verts_between(n, from, to, dim)
    if n <= 0
        return
    end
    edge_verts = [LinRange(from[i], to[i], n+2)[2: 1+n] for i in 1:dim]
    return edge_verts
end


"""
    triangle_vtk_order(corner_verts, order, dim, skip = false)

Compute the coordinates of a `VTK_LAGRANGE_TRIANGLE` of a triangle or order `order` 
defined by the coordinates of the vertices given in `corner_verts`. `dim` is the
dimension of the coordinates given. If `skip` is set to true, the coordinates
of the vertex- and edge-points aren't computed, which can be used to compute
points of a `VTK_LAGRANGE_WEDGE`
Inspired by: https://github.com/ju-kreber/paraview-scripts/blob/master/node_ordering.py
"""
function triangle_vtk_order(corner_verts, order, dim, skip = false)
    if order < 0
        return
    end
    if order == 0
        # For recursion
        return corner_verts[:,1]
    end
    #Corner vertices
    coords = Matrix{Float64}(undef, 0, dim)
    if skip == false
        coords = corner_verts
    end
    #edge vertices
    num_verts_on_edge = order -1 
    edges = [(1,2), (2,3), (3,1)]
    for (frm, to) in edges
        if skip == false
            tmp = n_verts_between(num_verts_on_edge, corner_verts[:, frm], corner_verts[:, to], dim)
            tmp_vec = Vector{Float64}(undef, dim)
            for i in 1:num_verts_on_edge
                for j in 1:dim
                    tmp_vec[j] = tmp[j][i]
                end
                coords = hcat(coords, tmp_vec)
            end
        end
    end
    if order == 2
        return coords
    end
    #faces
    e_x = (corner_verts[:,2] .- corner_verts[:,1]) ./ order
    e_y = (corner_verts[:,3] .- corner_verts[:,1]) ./ order
    inc = [(e_x + e_y) (-2*e_x + e_y) (e_x - 2*e_y)]
    if order >= 3
        coords = hcat(coords, triangle_vtk_order(corner_verts .+ inc, order - 3, dim, false))
    end
    return coords
end

"""
    vtk_order(elem::Tri, order)

Construct all node-points of a VTK_LAGRANGE_TRIANGLE of order order. The corner-nodes are
given by the reference-triangle used by StartUpDG
"""
function vtk_order(elem::Tri, order)
    tri_sud_vertices = [-1.0 1.0 -1.0; -1.0 -1.0 1.0]
    return triangle_vtk_order(tri_sud_vertices, order, 2)
end

"""
    SUD_to_vtk_order(rd::RefElemData, dim)

Compute the permutation of the nodes between StartUpDG and VTK
"""
function SUD_to_vtk_order(rd::RefElemData, dim)
    #nodes in vtk order
    vtk_nodes = vtk_order(rd.element_type, rd.N)
    vtk_formatted = Tuple(vtk_nodes[i,:] for i in 1:dim)
    
    #nodes in StartUpDG order
    interpolate = vandermonde(rd.element_type, rd.N, equi_nodes(rd.element_type, rd.N)...) / rd.VDM
    equi_dist_vertices = Tuple(interpolate * rd.rst[i] for i in 1:dim)

    #permutation
    return match_coordinate_vectors(vtk_formatted, equi_dist_vertices)
end

"""
    type_to_vtk(elem::Tri)
    return the VTK-type
"""
function type_to_vtk(elem::Tri)
    return VTKCellTypes.VTK_LAGRANGE_TRIANGLE
end
