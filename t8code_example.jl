using MPI
using T8code
using T8code.Libt8: sc_init, sc_free, sc_finalize, sc_array_new_data, sc_array_destroy
using T8code.Libt8: SC_LP_ESSENTIAL, SC_LP_PRODUCTION

# This is our own defined data that we will pass on to the
# adaptation callback.
mutable struct t8_step3_adapt_data_t
  midpoint                    :: NTuple{3,Cdouble}
  refine_if_inside_radius     :: Cdouble
  coarsen_if_outside_radius   :: Cdouble 
end

# The adaptation callback function. This function will be called once for each element
# and the return value decides whether this element should be refined or not.
#   return > 0 -> This element should get refined.
#   return = 0 -> This element should not get refined.
# If the current element is the first element of a family (= all level l elements that arise from refining
# the same level l-1 element) then this function is called with the whole family of elements
# as input and the return value additionally decides whether the whole family should get coarsened.
#   return > 0 -> The first element should get refined.
#   return = 0 -> The first element should not get refined.
#   return < 0 -> The whole family should get coarsened.
#  
# \param [in] forest       The current forest that is in construction.
# \param [in] forest_from  The forest from which we adapt the current forest (in our case, the uniform forest)
# \param [in] which_tree   The process local id of the current tree.
# \param [in] lelement_id  The tree local index of the current element (or the first of the family).
# \param [in] ts           The refinement scheme for this tree's element class.
# \param [in] is_family    if 1, the first \a num_elements entries in \a elements form a family. If 0, they do not.
# \param [in] num_elements The number of entries in \a elements elements that are defined.
# \param [in] elements     The element or family of elements to consider for refinement/coarsening.
function t8_step3_adapt_callback(forest, forest_from, which_tree, lelement_id,
                                  ts, is_family, num_elements, elements_ptr) :: Cint
  # Our adaptation criterion is to look at the midpoint coordinates of the current element and if
  # they are inside a sphere around a given midpoint we refine, if they are outside, we coarsen.

  centroid = Vector{Cdouble}(undef,3) # Will hold the element midpoint.
  # In t8_step3_adapt_forest we pass a t8_step3_adapt_data pointer as user data to the
  # t8_forest_new_adapt function. This pointer is stored as the used data of the new forest
  # and we can now access it with t8_forest_get_user_data (forest).
  adapt_data_ptr = Ptr{t8_step3_adapt_data_t}(t8_forest_get_user_data(forest))

  # You can use assert for assertions that are active in debug mode (when configured with --enable-debug).
  # If the condition is not true, then the code will abort.
  # In this case, we want to make sure that we actually did set a user pointer to forest and thus
  # did not get the NULL pointer from t8_forest_get_user_data.
  @T8_ASSERT(adapt_data_ptr != C_NULL)

  adapt_data = unsafe_load(adapt_data_ptr)

  elements = unsafe_wrap(Array, elements_ptr, num_elements)

  # Compute the element's centroid coordinates.
  t8_forest_element_centroid(forest_from, which_tree, elements[1], pointer(centroid))

  # Compute the distance to our sphere midpoint.
  dist = t8_vec_dist(centroid, Ref(adapt_data.midpoint[1]))
  if dist < adapt_data.refine_if_inside_radius
    # Refine this element.
    return 1
  elseif is_family == 1 && dist > adapt_data.coarsen_if_outside_radius
    # Coarsen this family. Note that we check for is_family before, since returning < 0
    # if we do not have a family as input is illegal. 
    return -1
  end

  # Do not change this element.
  return 0
end

# Adapt a forest according to our t8_step3_adapt_callback function.
# This will create a new forest and return it.
function t8_step3_adapt_forest(forest)
  adapt_data = t8_step3_adapt_data_t(
    (0.5, 0.5, 1.0),      # Midpoints of the sphere.
    0.2,                  # Refine if inside this radius.
    0.4                   # Coarsen if outside this radius.
  )

  # Check that forest is a committed, that is valid and usable, forest.
  @T8_ASSERT(t8_forest_is_committed(forest) == 1)

  # Create a new forest that is adapted from \a forest with our adaptation callback.
  # We provide the adapt_data as user data that is stored as the used_data pointer of the
  # new forest (see also t8_forest_set_user_data).
  # The 0, 0 arguments are flags that control
  #   recursive  -    If non-zero adaptation is recursive, thus if an element is adapted the children
  #                   or parents are plugged into the callback again recursively until the forest does not
  #                   change any more. If you use this you should ensure that refinement will stop eventually.
  #                   One way is to check the element's level against a given maximum level.
  #   do_face_ghost - If non-zero additionally a layer of ghost elements is created for the forest.
  #                   We will discuss ghost in later steps of the tutorial.
  forest_adapt = t8_forest_new_adapt(forest, @t8_adapt_callback(t8_step3_adapt_callback), 0, 0, Ref(adapt_data))

  return forest_adapt
end

# Print the local and global number of elements of a forest.
function t8_step3_print_forest_information(forest)
  # Check that forest is a committed, that is valid and usable, forest.
  @T8_ASSERT(t8_forest_is_committed(forest) == 1)

  # Get the local number of elements.
  local_num_elements = t8_forest_get_local_num_elements(forest)
  # Get the global number of elements.
  global_num_elements = t8_forest_get_global_num_elements(forest)

  t8_global_productionf(" [step3] Local number of elements:\t\t%i\n", local_num_elements)
  t8_global_productionf(" [step3] Global number of elements:\t%li\n", global_num_elements)
end


# In this function we first allocate a new uniformly refined forest at given
# refinement level. Then a second forest is created, where user data for the
# adaption call (cf. step 3) is registered.  The second forest inherts all
# properties of the first ("root") forest and deallocates it. The final
# adapted and commited forest is returned back to the calling scope.
function t8_step6_build_forest(element_type, comm, level)
  
  if element_type isa StartUpDG.Tri
    cmesh = t8_cmesh_new_hypercube(T8_ECLASS_TRIANGLE, comm, 0, 0, 0)

    # cmesh = t8_cmesh_new_periodic_tri(comm)

  elseif element_type isa StartUpDG.Quad
    cmesh = t8_cmesh_new_hypercube(T8_ECLASS_QUAD, comm, 0, 0, 0)

    # dim = 2
    # cmesh = t8_cmesh_new_periodic(comm, dim)  
  else
    @show element_type
  end

  scheme = t8_scheme_new_default_cxx()

  adapt_data = t8_step3_adapt_data_t(
    (0.0, 0.0, 0.0),      # Midpoints of the sphere.
    0.5,                  # Refine if inside this radius.
    0.7                   # Coarsen if outside this radius.
  )

  # Start with a uniform forest.
  forest = t8_forest_new_uniform(cmesh, scheme, level, 0, comm)

  forest_apbg_ref = Ref(t8_forest_t())
  t8_forest_init(forest_apbg_ref)
  forest_apbg = forest_apbg_ref[]

  # Adapt, partition, balance and create ghost elements all in one go.
  # See steps 3 and 4 for more details.
  t8_forest_set_user_data(forest_apbg, Ref(adapt_data))
  t8_forest_set_adapt(forest_apbg, forest, @t8_adapt_callback(t8_step3_adapt_callback), 0)
  t8_forest_set_partition(forest_apbg, C_NULL, 0)
  t8_forest_set_balance(forest_apbg, C_NULL, 0)
  t8_forest_set_ghost(forest_apbg, 1, T8_GHOST_FACES)
  t8_forest_commit(forest_apbg)

  return forest_apbg
end


# The uniform refinement level of the forest.
dim = 2
level = 2

using StartUpDG
rd = RefElemData(Quad(), 3)


# ========= Initialization ===========

mpiret = MPI.Init() # Initialize MPI. This has to happen before we initialize sc or t8code.
comm = MPI.COMM_WORLD # We will use MPI_COMM_WORLD as a communicator.
sc_init(comm, 1, 1, C_NULL, SC_LP_ESSENTIAL) # Initialize the sc library, has to happen before we initialize t8code.
t8_init(SC_LP_PRODUCTION) # Initialize t8code with log level SC_LP_PRODUCTION. See sc.h for more info on the log levels.

# Initialize an adapted forest 
forest = t8_step6_build_forest(rd.element_type, comm, level)

# Check that forest is a committed, that is valid and usable, forest.
@T8_ASSERT(t8_forest_is_committed(forest) == 1)

# ======== Extract StartUpDG quantities ========= 

function count_num_elements(forest)
  # Get the number of trees that have elements of this process. 
  num_local_trees = t8_forest_get_num_local_trees(forest)

  # count total number of elements
  num_elements = 0
  for itree = 0:num_local_trees-1
    num_elements_in_tree = t8_forest_get_tree_num_elements(forest, itree)
    num_elements += num_elements_in_tree
  end
  return num_elements
end

num_elements = count_num_elements(forest)
faces_per_element = zeros(Int, num_elements)

# coordinates
VX = [zeros(StartUpDG.num_vertices(rd.element_type)) for _ in 1:num_elements]
VY = [zeros(StartUpDG.num_vertices(rd.element_type)) for _ in 1:num_elements]

# count faces per element and get coordinates
current_element = 1
num_mortar_faces = 0
num_local_trees = t8_forest_get_num_local_trees(forest)
for itree = 0:num_local_trees-1
  tree_class = t8_forest_get_tree_class(forest, itree)
  eclass_scheme = t8_forest_get_eclass_scheme(forest, tree_class)
  num_elements_in_tree = t8_forest_get_tree_num_elements(forest, itree)
  for ielement = 0:num_elements_in_tree-1
    element = t8_forest_get_element_in_tree(forest, itree, ielement)
    num_faces = t8_element_num_faces(eclass_scheme, element)

    num_corners = t8_element_num_corners(eclass_scheme, element)
    for corner_number in 1:num_corners 
      coordinates = Vector{Cdouble}(undef, 3) 
      t8_forest_element_coordinate(forest, itree, element, corner_number-1, pointer(coordinates))
      VX[current_element][corner_number] = coordinates[1]
      VY[current_element][corner_number] = coordinates[2]
    end

    for iface = 1:num_faces
      neighids_ref = Ref{Ptr{t8_locidx_t}}()
      neighbors_ref = Ref{Ptr{Ptr{t8_element}}}()
      neigh_scheme_ref = Ref{Ptr{t8_eclass_scheme}}()
      dual_faces_ref = Ref{Ptr{Cint}}()
      num_neighbors_ref = Ref{Cint}()
      forest_is_balanced = Cint(1)

      t8_forest_leaf_face_neighbors(forest, itree, element,
          neighbors_ref, iface-1, dual_faces_ref, num_neighbors_ref,
          neighids_ref, neigh_scheme_ref, forest_is_balanced)

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

# t8code supplies these
EToE = [[[1] for _ in 1:faces_per_element[e]] for e in 1:num_elements]
EToF = [[[1] for _ in 1:faces_per_element[e]] for e in 1:num_elements]

# Loop over all local trees in the forest. 
current_element = 1
nonconforming_face_offset = 0
for itree = 0:num_local_trees-1
  tree_class = t8_forest_get_tree_class(forest, itree)
  eclass_scheme = t8_forest_get_eclass_scheme(forest, tree_class)
  num_elements_in_tree = t8_forest_get_tree_num_elements(forest, itree)
  for ielement = 0:num_elements_in_tree-1
    element = t8_forest_get_element_in_tree(forest, itree, ielement)
    num_faces = t8_element_num_faces(eclass_scheme, element)
    for iface = 1:num_faces
      neighids_ref = Ref{Ptr{t8_locidx_t}}()
      neighbors_ref = Ref{Ptr{Ptr{t8_element}}}()
      neigh_scheme_ref = Ref{Ptr{t8_eclass_scheme}}()

      dual_faces_ref = Ref{Ptr{Cint}}()
      num_neighbors_ref = Ref{Cint}()
      forest_is_balanced = Cint(1)
      t8_forest_leaf_face_neighbors(forest, itree, element,
        neighbors_ref, iface-1, dual_faces_ref, num_neighbors_ref,
        neighids_ref, neigh_scheme_ref, forest_is_balanced)

      num_neighbors = num_neighbors_ref[]
      dual_faces    = 1 .+ unsafe_wrap(Array, dual_faces_ref[], num_neighbors)
      neighids      = 1 .+ unsafe_wrap(Array, neighids_ref[], num_neighbors)

      EToE[current_element][iface] = neighids
      EToF[current_element][iface] = dual_faces

      current_face = iface + face_offsets[current_element]
      if num_neighbors == 1 # then it's a conforming face

        neighbor_face = dual_faces[1] + face_offsets[neighids[1]]

        
        FToF[current_face] = neighbor_face

      elseif num_neighbors > 1 # then it's a non-conforming face

        # add the current face index to the list of non-conforming faces (to split)
        push!(nonconforming_faces, current_face)
        
        # if it's a non-conforming face with 2:1 balance, the neighboring faces
        # are conforming (e.g., non-split) faces. 
        neighbor_faces = dual_faces .+ face_offsets[neighids]

        # split faces are ordered after un-split faces, so we
        # track the total number of conforming faces. 
        split_faces_indices = 
          @. num_element_faces + nonconforming_face_offset + (1:num_neighbors)

        nonconforming_face_offset += num_neighbors 

        # make connections between mortar faces
        FToF[split_faces_indices] .= neighbor_faces
      end

      # Free allocated memory.
      sc_free(t8_get_package_id(), neighbors_ref[])
      sc_free(t8_get_package_id(), dual_faces_ref[])
      sc_free(t8_get_package_id(), neighids_ref[])
    end

    current_element += 1
    
  end
end

# Clean-up
t8_forest_unref(Ref(forest)) # Destroy the forest.
t8_global_productionf(" Destroyed forest.\n")
sc_finalize()


# permute indices for StartUpDG ordering
function permute_t8code_vertices!(::Quad, VXY)
    VX, VY = VXY
    p = [1, 3, 2, 4]
    for e in eachindex(VX, VY)
        VX[e] = VX[e][p]
        VY[e] = VY[e][p]
    end
    return VX, VY
end

function permute_t8code_vertices!(::Tri, VXY)
    VX, VY = VXY

    for e in eachindex(VX, VY)
        vx, vy = VX[e], VY[e]
        area = StartUpDG.compute_triangle_area(zip(vx, vy))
        if area < 0
            VX[e] = vx[[2, 1, 3]]
            VY[e] = vy[[2, 1, 3]]
        end
    end
    return VX, VY
end

md = MeshData(permute_t8code_vertices!(rd.element_type, (VX, VY)), FToF, nonconforming_faces, rd)

using Plots
scatter(md.xyz..., leg=false, ms=3); plot!(rd, md)