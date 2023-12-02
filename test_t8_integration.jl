using MPI
using T8code
using T8code.Libt8: sc_init, sc_free, sc_finalize, sc_array_new_data, sc_array_destroy
using T8code.Libt8: SC_LP_ESSENTIAL, SC_LP_PRODUCTION

mpiret = MPI.Init() # Initialize MPI. This has to happen before we initialize sc or t8code.
comm = MPI.COMM_WORLD # We will use MPI_COMM_WORLD as a communicator.
sc_init(comm, 1, 1, C_NULL, SC_LP_ESSENTIAL) # Initialize the sc library, has to happen before we initialize t8code.
t8_init(SC_LP_PRODUCTION) # Initialize t8code with log level SC_LP_PRODUCTION. See sc.h for more info on the log levels.

using LinearAlgebra
using StartUpDG

using Plots
using StaticArrays
using Printf

using Logging

# Suppress some anoying warnings when using `scatter` plots.
Logging.disable_logging(Logging.Warn)

function adapt_callback(forest,
        forest_from,
        which_tree,
        lelement_id,
        ts,
        is_family,
        num_elements,
        elements_ptr)::Cint
    centroid = Vector{Cdouble}(undef, 3) # Will hold the element midpoint.
    elements = unsafe_wrap(Array, elements_ptr, num_elements)
    t8_forest_element_centroid(forest_from, which_tree, elements[1], pointer(centroid))

    level = t8_element_level(ts, elements[1])

    # Maximum refinement level
    if level >= 1
        return 0
    end

    if centroid[2] < centroid[1]
        return 1
    end

    return 0
end

function build_forest_hypercube(element_type, comm, level; do_adapt = false, dim = 2, is_periodic=false)

    if element_type isa StartUpDG.Tri
        cmesh = is_periodic ? t8_cmesh_new_periodic_tri(comm) :
                t8_cmesh_new_hypercube(T8_ECLASS_TRIANGLE, comm, 0, 0, 0)
    elseif element_type isa StartUpDG.Quad
        cmesh = is_periodic ? t8_cmesh_new_periodic(comm, dim) :
                t8_cmesh_new_hypercube(T8_ECLASS_QUAD, comm, 0, 0, 0)
    else
        @error element_type
    end

    # t8_cmesh_vtk_write_file(cmesh, "cmesh", 1.0)

    # Start with a uniform forest.
    scheme = t8_scheme_new_default_cxx()
    forest = t8_forest_new_uniform(cmesh, scheme, level, 0, comm)

    # t8_forest_write_vtk_ext(forest, "forest", 1, 1, 1, 1, 0, 1, 0, 0, C_NULL)
    # t8_forest_write_vtk(forest, "forest")

    if !do_adapt
        return forest
    end

    forest_apbg_ref = Ref(t8_forest_t())
    t8_forest_init(forest_apbg_ref)
    forest_apbg = forest_apbg_ref[]

    # Adapt, partition, balance and create ghost elements all in one go.
    # t8_forest_set_user_data(forest_apbg, Ref(adapt_data))
    t8_forest_set_adapt(forest_apbg, forest, @t8_adapt_callback(adapt_callback), 1)
    t8_forest_set_partition(forest_apbg, C_NULL, 0)
    t8_forest_set_balance(forest_apbg, C_NULL, 0)
    t8_forest_set_ghost(forest_apbg, 1, T8_GHOST_FACES)
    t8_forest_commit(forest_apbg)

    # t8_forest_write_vtk_ext(forest_apbg, "forest", 1, 1, 1, 1, 0, 1, 0, 0, C_NULL)

    return forest_apbg
end

# TODO: 
# make face ordering consistent with t8code; StartUpDG is ccw from bottom
# fix md.xf not being the right dimension in nonconforming.jl
# 
# 
# ? Aspirational mapP construction: 
# mapM = reshape(1:length(FToF) * rd.Nfq, rd.Nfq, :)
# mapP[:, f] = mapM[:, nbr_face] # without accounting for orientation 
# mapP[:, f] = mapP[p, f] # accounting for orientation 
# mapP[i, f] = mapM[p[f][i], fnbr[f]]
# 
# struct T8MappingArray <: AbstractMatrix{Int}
#     forest::...
#     mapM
# end
# getindex(x::T8MappingArray, i, f) = mapM[p[f][i], fnbr[f]]
# 
# * run update! after every refinement of StartUpDG
# update!(x::T8MappingArray)

# Note: `iface` must be one-indexed.
@inline function map_iface(::Tri, iface)
    mapping = [2, 3, 1]
    return mapping[iface]
    # return iface
end

# Note: `iface` must be one-indexed.
@inline function map_iface(::Quad, iface)
    return iface
end

# Note: `iface` must be one-indexed.
@inline function flip_faces(iface)
    mapping = [1, 3, 2, 4]
    return mapping[iface]
end

# Note: `iface` must be zero-indexed.
@inline function t8_flip_faces(iface)
    mapping = [0, 2, 1, 3]
    return mapping[iface + 1]
end

# Note: `itree` must be zero-indexed.
function t8_forest_element_vertices(forest, itree, eclass_scheme, element)
    num_corners = t8_element_num_corners(eclass_scheme, element)

    vertices = Matrix{Cdouble}(undef, 3, num_corners)

    for corner_number in 1:num_corners
        t8_forest_element_coordinate(forest,
            itree,
            element,
            corner_number - 1,
            @view(vertices[:, corner_number]))
    end

    u = [vertices[1, 2] - vertices[1, 1], vertices[2, 2] - vertices[2, 1], 0.0]
    v = [vertices[1, 3] - vertices[1, 1], vertices[2, 3] - vertices[2, 1], 0.0]
    w = [0.0, 0.0, 1.0]

    vol = dot(cross(u, v), w)

    # In case the volume is negative, we flip order of the vertices.
    if vol < 0.0
        # @error "Discovered negative volumes in `cmesh`: vol = $vol"
        vertices[:, 3], vertices[:, 2] = vertices[:, 2], vertices[:, 3]
    end

    return vertices, vol < 0.0
end

# Note: `itree` and `iface` must be zero-indexed.
function t8_forest_canonical_leaf_face_neighbors(forest,
        itree,
        eclass_scheme,
        element,
        iface)
    neighids_ref = Ref{Ptr{t8_locidx_t}}()
    neighbors_ref = Ref{Ptr{Ptr{t8_element}}}()
    neigh_scheme_ref = Ref{Ptr{t8_eclass_scheme}}()

    dual_faces_ref = Ref{Ptr{Cint}}()
    num_neighbors_ref = Ref{Cint}()
    forest_is_balanced = Cint(1)

    t8_forest_leaf_face_neighbors(forest,
        itree,
        element,
        neighbors_ref,
        iface,
        dual_faces_ref,
        num_neighbors_ref,
        neighids_ref,
        neigh_scheme_ref,
        forest_is_balanced)

    num_neighbors = num_neighbors_ref[]
    dual_faces = copy(unsafe_wrap(Array, dual_faces_ref[], num_neighbors))
    neighids = copy(unsafe_wrap(Array, neighids_ref[], num_neighbors))
    neighbors = copy(unsafe_wrap(Array, neighbors_ref[], num_neighbors))

    # Free allocated memory.
    sc_free(t8_get_package_id(), neighbors_ref[])
    sc_free(t8_get_package_id(), dual_faces_ref[])
    sc_free(t8_get_package_id(), neighids_ref[])

    # Compute the `orientation` of the touching faces.
    if t8_element_is_root_boundary(eclass_scheme, element, iface) == 1
        cmesh = t8_forest_get_cmesh(forest)
        # itree_in_cmesh = t8_forest_itree_to_cmesh_itree(forest, itree)
        itree_in_cmesh = itree
        iface_in_tree = t8_element_tree_face(eclass_scheme, element, iface)
        orientation_ref = Ref{Cint}(0)

        t8_cmesh_get_face_neighbor(cmesh,
            itree_in_cmesh,
            iface_in_tree,
            C_NULL,
            orientation_ref)
        orientation = orientation_ref[]
    else
        orientation = zero(Cint)
    end

    for ineighbor in 1:num_neighbors
        # if element has negative volume 
        _, flipped = t8_forest_element_vertices(forest,
            itree,
            eclass_scheme,
            neighbors[ineighbor])
        if flipped
            dual_faces[ineighbor] = t8_flip_faces(dual_faces[ineighbor])
        end
    end

    return num_neighbors, neighids, neighbors, dual_faces, orientation
end

function compute_connectivity(forest, rd)
    num_elements = t8_forest_get_local_num_elements(forest)

    # coordinates
    VX = [zeros(StartUpDG.num_vertices(rd.element_type)) for _ in 1:num_elements]
    VY = [zeros(StartUpDG.num_vertices(rd.element_type)) for _ in 1:num_elements]

    # count faces per element and get coordinates
    faces_per_element = zeros(Int, num_elements)
    current_element = 1
    num_mortar_faces = 0
    num_local_trees = t8_forest_get_num_local_trees(forest)

    for itree in 0:(num_local_trees - 1)
        tree_class = t8_forest_get_tree_class(forest, itree)
        eclass_scheme = t8_forest_get_eclass_scheme(forest, tree_class)
        num_elements_in_tree = t8_forest_get_tree_num_elements(forest, itree)

        for ielement in 0:(num_elements_in_tree - 1)
            element = t8_forest_get_element_in_tree(forest, itree, ielement)
            num_faces = t8_element_num_faces(eclass_scheme, element)

            num_corners = t8_element_num_corners(eclass_scheme, element)

            vertices, _ = t8_forest_element_vertices(forest, itree, eclass_scheme, element)

            for corner_number in 1:num_corners
                VX[current_element][corner_number] = vertices[1, corner_number]
                VY[current_element][corner_number] = vertices[2, corner_number]
            end

            for iface in 1:num_faces
                neighids_ref = Ref{Ptr{t8_locidx_t}}()
                neighbors_ref = Ref{Ptr{Ptr{t8_element}}}()
                neigh_scheme_ref = Ref{Ptr{t8_eclass_scheme}}()
                dual_faces_ref = Ref{Ptr{Cint}}()
                num_neighbors_ref = Ref{Cint}()
                forest_is_balanced = Cint(1)

                t8_forest_leaf_face_neighbors(forest,
                    itree,
                    element,
                    neighbors_ref,
                    iface - 1,
                    dual_faces_ref,
                    num_neighbors_ref,
                    neighids_ref,
                    neigh_scheme_ref,
                    forest_is_balanced)

                num_neighbors = num_neighbors_ref[]

                if num_neighbors > 1
                    num_mortar_faces += num_neighbors
                end

                sc_free(t8_get_package_id(), neighbors_ref[])
                sc_free(t8_get_package_id(), dual_faces_ref[])
                sc_free(t8_get_package_id(), neighids_ref[])
            end

            faces_per_element[current_element] = num_faces
            current_element += 1
        end
    end

    face_offsets = cumsum(faces_per_element) .- faces_per_element[1]

    num_element_faces = sum(faces_per_element)
    FToF = collect(1:(num_element_faces + num_mortar_faces))
    nonconforming_faces = Int[]

    orientations = zeros(Int8, length(FToF))

    # Loop over all local trees in the forest. 
    current_element = 1
    nonconforming_face_offset = 0
    for itree in 0:(num_local_trees - 1)
        tree_class = t8_forest_get_tree_class(forest, itree)
        eclass_scheme = t8_forest_get_eclass_scheme(forest, tree_class)
        num_elements_in_tree = t8_forest_get_tree_num_elements(forest, itree)
        for ielement in 0:(num_elements_in_tree - 1)
            element = t8_forest_get_element_in_tree(forest, itree, ielement)
            num_faces = t8_element_num_faces(eclass_scheme, element)
            for iface in 1:num_faces
                num_neighbors, neighids, neighbors, dual_faces, orientation = t8_forest_canonical_leaf_face_neighbors(forest,
                    itree,
                    eclass_scheme,
                    element,
                    iface - 1)

                # Convert to 1-based indexing.
                dual_faces .+= 1
                neighids .+= 1

                # Flipping node ordering is necessary for triangle meshes.
                _, flipped = t8_forest_element_vertices(forest,
                    itree,
                    eclass_scheme,
                    element)
                if flipped
                    iface = flip_faces(iface)
                end

                # Convert to StartUpDG face ordering convention.
                iface = map_iface(rd.element_type, iface)
                dual_faces = map_iface(rd.element_type, dual_faces)

                current_face = iface + face_offsets[current_element]
                if num_neighbors == 1 # then it's a conforming face
                    neighbor_face = dual_faces[1] + face_offsets[neighids[1]]

                    FToF[current_face] = neighbor_face
                    orientations[current_face] = orientation

                    t8_productionf("itree = %d, ielement = (%d,%d), iface = (%d,%d), orientation = %d\n",
                        itree + 1,
                        current_element,
                        neighids[1],
                        iface,
                        dual_faces[1],
                        orientation)

                    # t8_productionf("itree = %d, ielement = %d, iface = %d, neighid = %d, dual_face = %d, curr_face = %d, neigh_face = %d\n",
                    # itree + 1, ielement + 1, iface, neighids[1], dual_faces[1], current_face, neighbor_face)

                elseif num_neighbors > 1 # then it's a non-conforming face

                    # add the current face index to the list of non-conforming faces (to split)
                    push!(nonconforming_faces, current_face)

                    # if it's a non-conforming face with 2:1 balance, the neighboring faces
                    # are conforming (e.g., non-split) faces. 
                    neighbor_faces = dual_faces .+ face_offsets[neighids]

                    # split faces are ordered after un-split faces, so we
                    # track the total number of conforming faces. 
                    split_faces_indices = @. num_element_faces + nonconforming_face_offset +
                                             (1:num_neighbors)

                    nonconforming_face_offset += num_neighbors

                    # make connections between mortar faces
                    FToF[split_faces_indices] .= neighbor_faces
                    FToF[neighbor_faces] .= split_faces_indices
                    orientations[split_faces_indices] .= orientation

                else
                    t8_productionf("itree = %d, ielement = (%d,%d), iface = (%d,%d), orientation = %d\n",
                        itree + 1,
                        current_element,
                        current_element,
                        iface,
                        iface,
                        0)
                end
            end

            current_element += 1

            println("")
        end
    end

    return VX, VY, FToF, nonconforming_faces, orientations
end

function compute_coordinates(element_type, forest, rd, md)
    num_local_trees = t8_forest_get_num_local_trees(forest)

    if element_type isa StartUpDG.Tri
        # This is the strange reference basis for triangles in t8code.
        e_r = [0.0, 1.0, 0.0]
        e_s = [1.0, 1.0, 0.0]
    elseif element_type isa StartUpDG.Quad
        e_r = [1.0, 0.0, 0.0]
        e_s = [0.0, 1.0, 0.0]
    else
        @error element_type
    end

    x = similar(md.x)
    y = similar(md.y)

    current_element = 1
    for itree in 0:(num_local_trees - 1)

        num_elements_in_tree = t8_forest_get_tree_num_elements(forest, itree)
        for ielement in 0:(num_elements_in_tree - 1)
            element = t8_forest_get_element_in_tree(forest, itree, ielement)

            for i in 1:length(rd.r)
                r = 0.5 * (1 + rd.r[i])
                s = 0.5 * (1 + rd.s[i])

                ref_coords = r * e_r + s * e_s
                out_coords = Vector{Cdouble}(undef, 3)

                t8_forest_element_from_ref_coords(forest,
                    itree,
                    element,
                    pointer(ref_coords),
                    1,
                    pointer(out_coords))

                x[i, current_element] = out_coords[1]
                y[i, current_element] = out_coords[2]
            end

            current_element += 1
        end
    end

    return x, y
end

# ========= Main =========== #

# The uniform refinement level of the forest.
dim = 2
level = 0
npoly = 17

nodes, weights = StartUpDG.gauss_lobatto_quad(0, 0, npoly);

# etype = Quad()
etype = Tri()
rd = RefElemData(etype, Polynomial(), npoly)

# Initialize an adapted forest. 
forest = build_forest_hypercube(rd.element_type, comm, level; do_adapt = true, is_periodic=false)

# Compute connectivity arrays.
VX, VY, FToF, nonconforming_faces, orientations = compute_connectivity(forest, rd)

# Clean-up. t8code is not needed from here on.
t8_forest_unref(Ref(forest)) # Destroy the forest.
sc_finalize()

# md = MeshData((VX, VY), FToF, nonconforming_faces, rd)

# # Transfer t8code's geometry information to MeshData.
# # Not supported yet by the registered T8code.jl package.
# # x, y = compute_coordinates(rd.element_type, forest, rd, md)
# # md = MeshData(rd, md, x, y)


# # Plot the mesh.
# scatter(x, y, leg = false, markersize = 1.0 * [1.0, 1.0, 1.0]);
# plot!(rd, md)
# savefig("mesh.png")

# # Plot all interpolation nodes.
# xp = x
# yp = y
# up = (x + y)

# scatter(xp,
#     yp,
#     up,
#     zcolor = up,
#     msw = 0,
#     leg = false,
#     markersize = 3.0 * [1.0, 1.0, 1.0],
#     ratio = 1.0,
#     c = :darkrainbow,
#     cam = (0, 90),
#     size = 4.0 .* (600, 600)
#     # clim = (-0.001, 0.001),
#     # cbar = true,
# )


# Create a StartUpDG mesh.
# construct element nodal coordinates
VX_local = VX
VY_local = VY

(; V1) = rd
x = zeros(size(V1, 1), length(VX_local))
y = zeros(size(V1, 1), length(VX_local))
for e in eachindex(VX_local)
    view(x, :, e) .= V1 * VX_local[e]
    view(y, :, e) .= V1 * VY_local[e]
end

xf, yf = (x -> reshape(rd.Vf * x, rd.Nfq ÷ rd.num_faces, :)).((x, y))

# (; mortar_interpolation_matrix) = md.mesh_type
mortar_interpolation_matrix, mortar_projection_matrix = StartUpDG.compute_mortar_operators(rd)
if length(nonconforming_faces) > 0
    xm, ym = (x -> reshape(mortar_interpolation_matrix * x, :, 2 * length(nonconforming_faces))).((xf[:, nonconforming_faces], yf[:, nonconforming_faces]))
    xM, yM = [xf xm], [yf ym]
else
    xM, yM = xf, yf
end

mapM = reshape(1:length(xM), size(xM))
mapP = copy(mapM)
for (f, fnbr) in enumerate(FToF)
    # if fnbr < 0, it is a split face and we do not compute connectivity
    if fnbr > 0 && f != fnbr 
        if orientations[f] == orientations[fnbr]
            mapP[end:-1:1, f] = mapM[:, fnbr]
        else
            mapP[:, f] = mapM[:, fnbr]
        end
    end
end

# check that all node coordinates match
xy = [[xM[i, j], yM[i, j]] for i in axes(xM, 1), j in axes(xM, 2)]
norm(norm.(xy .- xy[mapP]))

# function number!(xyz...)
#     annotate!(vec.(xyz)..., string.(1:length(first(xyz))))
# end

