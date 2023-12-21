using MPI
using T8code
using T8code.Libt8: sc_init, sc_free, sc_finalize, sc_array_new_data, sc_array_destroy
using T8code.Libt8: SC_LP_ESSENTIAL, SC_LP_PRODUCTION

# This is useful in live sessions.
# if T8code.Libt8.sc_is_initialized() == 1
if @isdefined(forest) && forest != C_NULL
  # Clean-up
  t8_forest_unref(Ref(forest)) # Destroy the forest.
  sc_finalize()
  forest = C_NULL
  println("Cleaned up.")
end

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

function adapt_callback(forest, forest_from, which_tree, lelement_id,
                                  ts, is_family, num_elements, elements_ptr) :: Cint
    
  centroid = Vector{Cdouble}(undef,3) # Will hold the element midpoint.
  elements = unsafe_wrap(Array, elements_ptr, num_elements)
  t8_forest_element_centroid(forest_from, which_tree, elements[1], pointer(centroid))

  level = t8_element_level(ts, elements[1])

  # Maximum refinement level
  if level >= max_level
    return 0
  end

  if centroid[2] > centroid[1]
    return 1
  end

  return 0
end

function build_forest_hypercube(element_type, comm, level; do_adapt = false, dim = 2)
  
  if element_type isa StartUpDG.Tri
    # cmesh = t8_cmesh_new_hypercube(T8_ECLASS_TRIANGLE, comm, 0, 0, 0)
    cmesh = t8_cmesh_new_periodic_tri(comm)

  elseif element_type isa StartUpDG.Quad
    cmesh = t8_cmesh_new_hypercube(T8_ECLASS_QUAD, comm, 0, 0, 0)
    # cmesh = t8_cmesh_new_periodic(comm, dim)

  else
    cmesh = t8_cmesh_new_periodic_hybrid(comm)

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

# Note: `iface` must be zero-indexed.
@inline function t8_flip_faces(iface)
  # This function works only for triangles.
  # It won't be reached for quads anyway since they do not have to be flipped.

  if iface > 2
    @error "iface = $iface"
  end

  mapping = [0, 2, 1]
  return mapping[iface+1]
end

# Note: `itree` must be zero-indexed.
function t8_forest_element_vertices_flipped(forest, itree, eclass_scheme, element)
    num_corners = t8_element_num_corners(eclass_scheme, element)

    vertices = Matrix{Cdouble}(undef, 3, num_corners) 
    vertices .= 0.0

    for corner_number in 1:num_corners 
      t8_forest_element_coordinate(forest, itree, element, corner_number-1, @view(vertices[:,corner_number]))
    end

    u = [vertices[1,2] - vertices[1,1], vertices[2,2] - vertices[2,1], 0.0]
    v = [vertices[1,3] - vertices[1,1], vertices[2,3] - vertices[2,1], 0.0]
    w = [0.0, 0.0, 1.0]

    vol = dot(cross(u, v), w)

    do_flipping = vol < 0.0

    # In case the volume is negative, we swap corners 2 and 3.
    if do_flipping
      vertices[:,3], vertices[:,2] = vertices[:,2], vertices[:,3]
    end

    return vertices, do_flipping
end

# Note: `itree` and `iface` must be zero-indexed.
function t8_forest_leaf_face_neighbors_flipped(forest, itree, eclass_scheme, element, iface)
    neighids_ref = Ref{Ptr{t8_locidx_t}}()
    neighbors_ref = Ref{Ptr{Ptr{t8_element}}}()
    neigh_scheme_ref = Ref{Ptr{t8_eclass_scheme}}()

    dual_faces_ref = Ref{Ptr{Cint}}()
    num_neighbors_ref = Ref{Cint}()
    forest_is_balanced = Cint(1)

    t8_forest_leaf_face_neighbors(forest, itree, element,
      neighbors_ref, iface, dual_faces_ref, num_neighbors_ref,
      neighids_ref, neigh_scheme_ref, forest_is_balanced)

    num_neighbors = num_neighbors_ref[]
    dual_faces    = copy(unsafe_wrap(Array, dual_faces_ref[], num_neighbors))
    neighids      = copy(unsafe_wrap(Array, neighids_ref[], num_neighbors))
    neighbors     = copy(unsafe_wrap(Array, neighbors_ref[], num_neighbors))

    # Free allocated memory.
    sc_free(t8_get_package_id(), neighbors_ref[])
    sc_free(t8_get_package_id(), dual_faces_ref[])
    sc_free(t8_get_package_id(), neighids_ref[])

    # Compute the `orientation` of the touching faces.
    if t8_element_is_root_boundary(eclass_scheme, element, iface) == 1
      cmesh = t8_forest_get_cmesh(forest)
      iface_in_tree = t8_element_tree_face(eclass_scheme, element, iface)
      orientation_ref = Ref{Cint}(0)
      t8_cmesh_get_face_neighbor(cmesh, itree, iface_in_tree, C_NULL, orientation_ref)
      orient = orientation_ref[]
    else
      orient = zero(Cint)
    end
    orientation = [orient for _ in 1:num_neighbors]

    for ineighbor = 1:num_neighbors
      neigh_face_ref = Ref{Cint}()
      neigh_itree = t8_forest_element_face_neighbor(forest, itree, element, neighbors[ineighbor], neigh_scheme_ref[], iface, neigh_face_ref)
      _, flipped = t8_forest_element_vertices_flipped(forest, neigh_itree, eclass_scheme, neighbors[ineighbor])
      if flipped
        # Only count flip when dual face is `0`. Flipping changes only face orientation for face `0`.
        if dual_faces[ineighbor] == 0
          orientation[ineighbor] = 1 - orientation[ineighbor]
        end
        dual_faces[ineighbor] = t8_flip_faces(dual_faces[ineighbor])
      end
    end

    return num_neighbors, neighids, neighbors, dual_faces, orientation
end

function compute_connectivity(forest)

  num_elements = t8_forest_get_local_num_elements(forest)

  vertices = []
  levels = []
  neighbor_ids = []
  dual_faces = []
  orientations = []

  current_element = 1

  # Loop over trees.
  num_local_trees = t8_forest_get_num_local_trees(forest)
  for itree = 0:num_local_trees-1
    tree_class = t8_forest_get_tree_class(forest, itree)
    eclass_scheme = t8_forest_get_eclass_scheme(forest, tree_class)
    num_elements_in_tree = t8_forest_get_tree_num_elements(forest, itree)

    # Loop over all elements in tree.
    for ielement = 0:num_elements_in_tree-1
      element = t8_forest_get_element_in_tree(forest, itree, ielement)
      level = t8_element_level(eclass_scheme, element)
      num_faces = t8_element_num_faces(eclass_scheme, element)
      num_corners = t8_element_num_corners(eclass_scheme, element)
        
      push!(levels, level)

      # Compute the neighbors at face and push them to `local_*` arrays.
      local_neighbor_ids = []
      local_dual_faces = []
      local_orientations = []
      for iface = 1:num_faces
        face_num_neighbors, face_neighbor_ids, face_neighbors, face_dual_faces, face_orientation =
          t8_forest_leaf_face_neighbors_flipped(forest, itree, eclass_scheme, element, iface-1)

        # Convert to 1-based indexing.
        face_dual_faces .+= 1
        face_neighbor_ids .+= 1

        push!(local_neighbor_ids, face_neighbor_ids)
        push!(local_dual_faces, face_dual_faces)
        push!(local_orientations, face_orientation)
      end

      # Compute current element vertices.
      element_vertices, element_flipped = t8_forest_element_vertices_flipped(forest, itree, eclass_scheme, element)
      
      # Just copy `element_vertices` to `loca_vertices`.
      local_vertices = []
      for corner_number in 1:num_corners 
        push!(local_vertices, element_vertices[:,corner_number])
      end
      push!(vertices, local_vertices)
    
      # Apply flipping of current element.
      if element_flipped
        local_neighbor_ids[3], local_neighbor_ids[2] = local_neighbor_ids[2], local_neighbor_ids[3]
        local_dual_faces[3], local_dual_faces[2] = local_dual_faces[2], local_dual_faces[3]
        local_orientations[3], local_orientations[2] = local_orientations[2], local_orientations[3]
        for (i,_) in enumerate(local_orientations[1])
          local_orientations[1][i] = 1-local_orientations[1][i]
        end
      end

      push!(neighbor_ids, local_neighbor_ids)
      push!(dual_faces, local_dual_faces)
      push!(orientations, local_orientations)

      current_element += 1
    end
  end

  return vertices, levels, neighbor_ids, dual_faces, orientations
end

# ========= Main =========== #

# The uniform refinement level of the forest.
ini_level = 0
max_level = 1
figsize = 1000 .*(1,1) # increase if things are getting crowded

# etype = Quad()
etype = Tri()
# etype = :hybrid # This actually works! :)

# Initialize an adapted forest. 
forest = build_forest_hypercube(etype, comm, ini_level; do_adapt = true)

# Compute connectivity arrays.
vertices, levels, neighbor_ids, dual_faces, orientations = compute_connectivity(forest)

# Clean-up. t8code is not needed from here on.
t8_forest_unref(Ref(forest)) # Destroy the forest.
forest = C_NULL
sc_finalize()

# Mapping of inner face ids to inner corner ids.
t8_iface_to_corners = Dict(
  :tri => [ [2,3], [1,3], [1,2] ],
  :quad => [ [1,3], [2,4], [1,2], [3,4] ],
)

t8_num_corners_to_etype = [ :point, :line, :tri, :quad ]

# local vertices
VX = map(x->getindex.(x, 1), vertices)
VY = map(x->getindex.(x, 2), vertices)