using StartUpDG
using MPI
using T8code
using T8code.Libt8: sc_init, sc_free, sc_finalize, sc_array_new_data, sc_array_destroy
using T8code.Libt8: SC_LP_ESSENTIAL, SC_LP_PRODUCTION

using StaticArrays: SVector
using Static

include("t8code_auxiliary.jl")
mutable struct T8codeMesh{NDIMS, RealT <: Real, IsParallel}               
    forest      :: Ptr{t8_forest} # cpointer to forest

    # TODO: add element type(s), is_periodic ??
    # element_type::ElemType - add to type parameters
    
    is_parallel :: IsParallel

    ninterfaces :: Int
    nmortars    :: Int
    nboundaries :: Int

    # nodes: interpolation nodes
    function T8codeMesh{NDIMS}(forest) where {NDIMS}
        is_parallel = False()

        mesh = new{NDIMS, Float64, typeof(is_parallel)}(forest, is_parallel)

        finalizer(mesh) do mesh
            # When finalizing `mesh.forest`, `mesh.scheme` and `mesh.cmesh` are
            # also cleaned up from within `t8code`. The cleanup code for
            # `cmesh` does some MPI calls for deallocating shared memory
            # arrays. Due to garbage collection in Julia the order of shutdown
            # is not deterministic. The following code might happen after MPI
            # is already in finalized state.
            # If the environment variable `TRIXI_T8CODE_SC_FINALIZE` is set the
            # `finalize_hook` of the MPI module takes care of the cleanup. See
            # further down. However, this might cause a pile-up of `mesh`
            # objects during long-running sessions.
            if !MPI.Finalized()
                t8_unref_forest(mesh.forest)
            end
        end

        # This finalizer call is only recommended during development and not for
        # production runs, especially long-running sessions since a reference to
        # the `mesh` object will be kept throughout the lifetime of the session.
        # See comments in `init_t8code()` in file `src/auxiliary/t8code.jl` for
        # more information.
        if haskey(ENV, "T8CODE_SC_FINALIZE")
            MPI.add_finalize_hook!() do
                t8_unref_forest(mesh.forest)
            end
        end

        return mesh
    end
end

const SerialT8codeMesh{NDIMS} = T8codeMesh{NDIMS, <:Real, <:False}
@inline mpi_parallel(mesh::SerialT8codeMesh) = False()

@inline Base.ndims(::T8codeMesh{NDIMS}) where {NDIMS} = NDIMS

@inline ntrees(mesh::T8codeMesh) = Int(t8_forest_get_num_local_trees(mesh.forest))
import StartUpDG: num_elements
@inline StartUpDG.num_elements(mesh::T8codeMesh) = Int(t8_forest_get_local_num_elements(mesh.forest))
@inline ninterfaces(mesh::T8codeMesh) = mesh.ninterfaces
@inline nmortars(mesh::T8codeMesh) = mesh.nmortars
@inline nboundaries(mesh::T8codeMesh) = mesh.nboundaries

function Base.show(io::IO, mesh::T8codeMesh)
    print(io, "T8codeMesh{", ndims(mesh), "}")
end

function Base.show(io::IO, ::MIME"text/plain", mesh::T8codeMesh)
    print(io, "T8codeMesh{", ndims(mesh), "}")
end

# ========= Initialization ===========

mpiret = MPI.Init() # Initialize MPI. This has to happen before we initialize sc or t8code.
comm = MPI.COMM_WORLD # We will use MPI_COMM_WORLD as a communicator.
sc_init(comm, 1, 1, C_NULL, SC_LP_ESSENTIAL) # Initialize the sc library, has to happen before we initialize t8code.
t8_init(SC_LP_PRODUCTION) # Initialize t8code with log level SC_LP_PRODUCTION. See sc.h for more info on the log levels.

# Initialize an adapted forest 
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
  
    # adapt_data = t8_step3_adapt_data_t(
    #   (0.0, 0.0, 0.0),      # Midpoints of the sphere.
    #   0.5,                  # Refine if inside this radius.
    #   0.7                   # Coarsen if outside this radius.
    # )
  
    # Start with a uniform forest.
    forest = t8_forest_new_uniform(cmesh, scheme, level, 0, comm)
    return forest
  
    # forest_apbg_ref = Ref(t8_forest_t())
    # t8_forest_init(forest_apbg_ref)
    # forest_apbg = forest_apbg_ref[]
  
    # # Adapt, partition, balance and create ghost elements all in one go.
    # # See steps 3 and 4 for more details.
    # t8_forest_set_user_data(forest_apbg, Ref(adapt_data))
    # t8_forest_set_adapt(forest_apbg, forest, @t8_adapt_callback(t8_step3_adapt_callback), 0)
    # t8_forest_set_partition(forest_apbg, C_NULL, 0)
    # t8_forest_set_balance(forest_apbg, C_NULL, 0)
    # t8_forest_set_ghost(forest_apbg, 1, T8_GHOST_FACES)
    # t8_forest_commit(forest_apbg)
  
    # return forest_apbg
end

  
# The uniform refinement level of the forest.
dim = 2
level = 2
forest = t8_step6_build_forest(Tri(), comm, level)

# Check that forest is a committed, that is valid and usable, forest.
@T8_ASSERT(t8_forest_is_committed(forest) == 1)

mesh = T8codeMesh{dim}(forest)

rd = RefElemData(Tri(), 3)

# compute physical coordinates

function set_coordinates!(forest, rd, md)
    num_local_trees = t8_forest_get_num_local_trees(forest)
    current_element = 1
    for itree = 0:num_local_trees-1
      cmesh = t8_forest_get_cmesh(forest)
      gitree = t8_forest_global_tree_id(forest, itree)
      num_elements_in_tree = t8_forest_get_tree_num_elements(forest, itree)
      for ielement = 0:num_elements_in_tree-1
        element = t8_forest_get_element_in_tree(forest, itree, ielement)
        for i = 1:length(rd.r)
          ref_coords = [0.5*(1 + rd.r[i]), 0.5*(1 + rd.s[i]),0.0]
          out_coords = Vector{Cdouble}(undef,3)
          t8_forest_element_from_ref_coords(forest, itree, element, pointer(ref_coords), 1, pointer(out_coords), C_NULL)
          x[i, current_element] = out_coords[1]
          y[i, current_element] = out_coords[2]
        end
        current_element += 1
      end
    end
  end
  