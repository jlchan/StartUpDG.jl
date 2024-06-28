module StartUpDG

using Reexport: @reexport

using ConstructionBase: ConstructionBase
using FillArrays: Fill
using HDF5: h5open # used to read in SBP triangular node data
using Kronecker: kronecker # for Hex element matrix manipulations
using LinearAlgebra:
    cond, diagm, eigvals, Diagonal, UniformScaling, I, mul!, norm, qr, ColumnNorm, Symmetric, nullspace, pinv
using NodesAndModes: meshgrid, find_face_nodes, face_vertices
@reexport using NodesAndModes # for basis functions
using PathIntersections: PathIntersections
@reexport using PathIntersections: PresetGeometries
using Printf: @sprintf
using RecipesBase: RecipesBase
using StaticArrays: SVector, SMatrix
using Setfield: setproperties, @set # for "modifying" structs (setproperties)
using SparseArrays: sparse, droptol!, blockdiag, nnz
using Triangulate: Triangulate, TriangulateIO, triangulate
@reexport using WriteVTK

@inline mean(x) = sum(x) / length(x)

include("RefElemData.jl")

include("RefElemData_polynomial.jl")
export RefElemData, Polynomial
export MultidimensionalQuadrature, TensorProductQuadrature
export TensorProductGaussCollocation, Gauss

include("RefElemData_TensorProductWedge.jl")
export TensorProductWedge

include("RefElemData_SBP.jl")
export SBP, DefaultSBPType, TensorProductLobatto, Hicken, Kubatko # types for SBP node dispatch
export LobattoFaceNodes, LegendreFaceNodes # type parameters for SBP{Kubatko{...}}
export hybridized_SBP_operators

include("low_order_sbp.jl")
export sparse_low_order_SBP_operators
export subcell_limiting_operators
export inverse_trace_constant, face_type

include("ref_elem_utils.jl")

include("MeshData.jl")
export MeshData, num_elements

include("geometric_functions.jl")
export geometric_factors, estimate_h

# spatial connectivity routines
include("connectivity_functions.jl")
export make_periodic

# helper array type for cut cell and hybrid meshes
include("named_array_partition.jl")
export NamedArrayPartition

# for tagging faces on boundaries
include("boundary_utils.jl")
export boundary_face_centroids, tag_boundary_faces, tag_boundary_nodes

include("hybrid_meshes.jl")
export num_faces, num_vertices, HybridMeshExample

include("physical_frame_basis.jl")
include("cut_cell_meshes.jl")
export PhysicalFrame, equi_nodes
export Subtriangulation 

# ! this will be deprecated in a future release
include("cut_cell_moment_fitting.jl") 
export MomentFitting

include("state_redistribution.jl")
export StateRedistribution, apply!

include("nonconforming.jl")
export num_mortars_per_face, NonConformingQuadMeshExample

# uniform meshes + face vertex orderings
include("mesh/simple_meshes.jl")
export uniform_mesh
include("mesh/gmsh_utilities.jl")
export read_Gmsh_2D # unifies v2.2.8 and v4.1 mesh reading
export readGmsh2D, readGmsh2D_v4 # TODO: deprecate
export read_Gmsh_2D_v2, read_Gmsh_2D_v4
export MeshImportOptions
include("mesh/hohqmesh_utilities.jl")
export read_HOHQMesh

# Plots.jl recipes for meshes
include("mesh/vtk_helper.jl")
include("mesh/mesh_visualization.jl")
export VertexMeshPlotter, MeshPlotter, MeshData_to_vtk

# Triangulate interfaces and pre-built meshes
include("mesh/triangulate_utils.jl")
export refine, triangulateIO_to_VXYEToV, get_node_boundary_tags
export BoundaryTagPlotter
include("mesh/triangulate_example_meshes.jl")
export triangulate_domain
export Scramjet, SquareDomain, RectangularDomain, RectangularDomainWithHole
export CircularDomain, PartialCircularDomain

# simple explicit time-stepping included for conveniencea
include("explicit_timestep_utils.jl")
export ck45 # LSERK 45

end # module
