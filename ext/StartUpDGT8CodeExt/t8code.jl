using StartUpDG

using MPI: MPI
using T8code
using T8code.Libt8: sc_init, sc_free, sc_finalize, sc_array_new_data, sc_array_destroy
using T8code.Libt8: SC_LP_ESSENTIAL, SC_LP_PRODUCTION

using StaticArrays: SVector
using Static: Static

include("t8code_auxiliary.jl")

# ========= Initialization ===========

mpiret = MPI.Init() # Initialize MPI. This has to happen before we initialize sc or t8code.
comm = MPI.COMM_WORLD # We will use MPI_COMM_WORLD as a communicator.
sc_init(comm, 1, 1, C_NULL, SC_LP_ESSENTIAL) # Initialize the sc library, has to happen before we initialize t8code.
t8_init(SC_LP_PRODUCTION) # Initialize t8code with log level SC_LP_PRODUCTION. See sc.h for more info on the log levels.

# create a mesh
function t8_uniform_forest(element_type, comm, level; is_periodic = false)
    if element_type isa StartUpDG.Tri
        cmesh = is_periodic ? t8_cmesh_new_periodic_tri(comm) :
                t8_cmesh_new_hypercube(T8_ECLASS_TRIANGLE, comm, 0, 0, 0)
    elseif element_type isa StartUpDG.Quad
        cmesh = is_periodic ? t8_cmesh_new_periodic(comm, dim) :
                t8_cmesh_new_hypercube(T8_ECLASS_QUAD, comm, 0, 0, 0)
    else
        error("Uniform t8 mesh for element type $element_type not yet implemented.")
    end

    # Start with a uniform forest.
    forest = t8_forest_new_uniform(cmesh, t8_scheme_new_default_cxx(), level, 0, comm)
    return forest
end

function convert_ref_coords_to_t8(::Tri, r, s)
    # map to biunit
    r = @. 0.5 * (1 + r)
    s = @. 0.5 * (1 + s)

    # used for triangular mappings 
    e_r = [0, 1, 0]
    e_s = [1, 1, 0]

    # t8 has a [r, r + s, ...] coordinate system
    return map((r, s) -> r * e_r + s * e_s, r, s)
end

function convert_ref_coords_to_t8(::Quad, r, s)
    # map to biunit
    r = @. 0.5 * (1 + r)
    s = @. 0.5 * (1 + s)

    return [[r_i, s_i] for (r_i, s_i) in zip(r, s)]
end

function map_reference_nodes_to_physical(mesh::T8codeMesh{NDIMS},
    element_type,
    rst...) where {NDIMS}
    t8_ref_coords = convert_ref_coords_to_t8(element_type, rst...)

    xyz = ntuple(_ -> zeros(length(t8_ref_coords), num_elements(mesh)), 2)

    out_coords = Vector{Cdouble}(undef, 3)
    num_local_trees = t8_forest_get_num_local_trees(forest)
    current_element = 1
    for itree in 0:(num_local_trees - 1)
        num_elements_in_tree = t8_forest_get_tree_num_elements(forest, itree)
        for ielement in 0:(num_elements_in_tree - 1)
            element = t8_forest_get_element_in_tree(forest, itree, ielement)
            for i in eachindex(t8_ref_coords)
                t8_forest_element_from_ref_coords(forest,
                    itree,
                    element,
                    pointer(t8_ref_coords[i]),
                    1,
                    pointer(out_coords),
                    C_NULL)
                for dim in 1:NDIMS
                    xyz[dim][i, current_element] = out_coords[dim]
                end
            end
            current_element += 1
        end
    end
    return xyz
end

rd = RefElemData(Tri(), 1)

# The uniform refinement level of the forest.
dim = 2
level = 1
forest = t8_uniform_forest(rd.element_type, comm, level)

# Check that forest is a committed, that is valid and usable, forest.
@T8_ASSERT(t8_forest_is_committed(forest)==1)
mesh = T8codeMesh{dim}(forest)

# map to unit triangle, permute r/s for t8code
x, y = map_reference_nodes_to_physical(mesh, rd.element_type, rd.rst...)

for e in axes(x, 2)
    rxJ, sxJ, ryJ, syJ, J = StartUpDG.geometric_factors(view(x, :, e), view(y, :, e), rd.Drst...)
    if any(@. J < 0)
        @show e
        # then we should permute x
    end
end

function annotate_scatter(xyz...)
    scatter(xyz...)
    annotate!((xyz..., string.(1:length(xyz[1]))))
end

# face_node_indices[face] = indices of nodes
function get_face_node_indices(rd::RefElemData)
    face_node_ids = reshape(1:(rd.Nfq), :, rd.num_faces)
    return [face_node_ids[:, f] for f in 1:(rd.num_faces)]
end

# routines to convert t8 face indexing to StartUpDG face indexing
@inline function convert_t8_face_index(::Tri, iface)
    # mapping = [2, 1, 3] 
    mapping = [2, 3, 1]
    return mapping[iface]
end

@inline convert_t8_face_index(::Quad, iface) = iface

face_node_indices = get_face_node_indices(rd)

mapM = reshape(1:length(xf), size(xf)...)
mapP = copy(mapM)

if true
    # Loop over all local trees in the forest. 
    current_element = 1
    num_local_trees = t8_forest_get_num_local_trees(forest)
    for itree in 0:(num_local_trees - 1)
        tree_class = t8_forest_get_tree_class(forest, itree)
        eclass_scheme = t8_forest_get_eclass_scheme(forest, tree_class)
        num_elements_in_tree = t8_forest_get_tree_num_elements(forest, itree)
        for ielement in 0:(num_elements_in_tree - 1)
            element = t8_forest_get_element_in_tree(forest, itree, ielement)
            num_faces = t8_element_num_faces(eclass_scheme, element)
            for iface in 1:num_faces
                neighids_ref = Ref{Ptr{t8_locidx_t}}()
                neighbors_ref = Ref{Ptr{Ptr{t8_element}}}()
                neigh_scheme_ref = Ref{Ptr{t8_eclass_scheme}}()

                dual_faces_ref = Ref{Ptr{Cint}}()
                num_neighbors_ref = Ref{Cint}()
                forest_is_balanced = Cint(1)
                t8_forest_leaf_face_neighbors(forest, itree, element,
                    neighbors_ref, iface - 1,
                    dual_faces_ref, num_neighbors_ref,
                    neighids_ref, neigh_scheme_ref,
                    forest_is_balanced)

                num_neighbors = num_neighbors_ref[]
                dual_faces = 1 .+ unsafe_wrap(Array, dual_faces_ref[], num_neighbors)
                neighids = 1 .+ unsafe_wrap(Array, neighids_ref[], num_neighbors)

                @assert num_neighbors <= 1 # assume it's a conforming or boundary face for now

                # neighbor_face = dual_faces[1] + face_offsets[neighids[1]]
                # FToF[current_face] = neighbor_face

                # TODO: map between StartUpDG and t8 face ordering

                if num_neighbors == 1
                    if true # iface==3 && current_element == 6
                        @show iface, current_element
                        @show dual_faces[1], neighids[1]

                        f = convert_t8_face_index(rd.element_type, iface)
                        @show [
                            xf[face_node_indices[f], current_element],
                            yf[face_node_indices[f], current_element],
                        ]

                        f = convert_t8_face_index(rd.element_type, dual_faces[1])
                        @show [
                            xf[face_node_indices[f], neighids[1]],
                            yf[face_node_indices[f], neighids[1]],
                        ]

                        # Compute the `orientation` of the touching faces.
                        if t8_element_is_root_boundary(eclass_scheme, element, iface-1) == 1
                            cmesh = t8_forest_get_cmesh(forest)
                            itree_in_cmesh = t8_forest_ltreeid_to_cmesh_ltreeid(forest, itree)
                            iface_in_tree = t8_element_tree_face(eclass_scheme, element, iface-1)
                            orientation_ref = Ref{Cint}()

                            t8_cmesh_get_face_neighbor(cmesh, itree_in_cmesh, iface_in_tree, C_NULL,
                                                        orientation_ref)
                            orientation = orientation_ref[]
                        else
                            orientation = zero(Cint)
                        end
                        @show orientation

                        println("")
                    end

                end
            end
            current_element += 1
        end
    end
end

gr(leg=false)

function compress_towards_mean(x_in)
    x = copy(x_in)
    for i in axes(x, 2)
        xavg = sum(x[:, i]) / size(x, 1)
        @. x[:, i] = xavg + .8 * (x[:, i] - xavg)
    end
    return vec(x)
end

t8_f, elem = 3, 6
nbr = 4

scatter(xf, yf)
xf_annotate, yf_annotate = compress_towards_mean.((xf, yf))
for i in mapM    
    if abs(xf[i]-0.5) < 0.5 && abs(yf[i]^2-0.5) < .50
        annotate!((xf_annotate[i], yf_annotate[i], string(mapM[i])))
    end
end
plot!()
f = convert_t8_face_index(rd.element_type, t8_f)
scatter!(xf[face_node_indices[f], elem], yf[face_node_indices[f], elem], ms=10)
annotate!((compress_towards_mean.((xf[:,elem], yf[:,elem]))..., string.(vec(mapM[:,elem]))))
annotate!((compress_towards_mean.((xf[:,nbr], yf[:,nbr]))..., string.(vec(mapM[:,nbr]))))

