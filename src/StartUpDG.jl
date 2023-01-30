module StartUpDG

using Reexport: @reexport

using ConstructionBase: ConstructionBase
using FillArrays: Ones, Zeros, Fill
using HDF5: h5open # used to read in SBP triangular node data
using Kronecker: kronecker # for Hex element matrix manipulations
using LinearAlgebra: cond, diagm, eigvals, Diagonal, I, mul!, norm, qr, ColumnNorm
using NodesAndModes: meshgrid, find_face_nodes, face_vertices
@reexport using NodesAndModes # for basis functions
using OrderedCollections: LittleDict # fast ordered dict for a small number of entries
using PathIntersections
@reexport using PathIntersections: PresetGeometries
using Printf: @sprintf
using RecipesBase: RecipesBase
using StaticArrays: SVector, SMatrix
using Setfield: setproperties, @set # for "modifying" structs (setproperties)
using SparseArrays: sparse, droptol!, blockdiag
using Triangulate: Triangulate, TriangulateIO, triangulate
using WriteVTK
@reexport using UnPack: @unpack  # for getting values in RefElemData and MeshData

# reference element utility functions
include("RefElemData.jl")
include("RefElemData_polynomial.jl")
include("RefElemData_SBP.jl")
include("ref_elem_utils.jl")
export RefElemData, Polynomial
export SBP, DefaultSBPType, TensorProductLobatto, Hicken, Kubatko # types for SBP node dispatch
export LobattoFaceNodes, LegendreFaceNodes # type parameters for SBP{Kubatko{...}}
export hybridized_SBP_operators, inverse_trace_constant, face_type

include("MeshData.jl")
export MeshData, num_elements

include("geometric_functions.jl")
export geometric_factors, estimate_h

# spatial connectivity routines
include("connectivity_functions.jl")
export make_periodic, match_coordinate_vectors, match_coordinate_vectors!

# for tagging faces on boundaries
include("boundary_utils.jl")
export boundary_face_centroids, tag_boundary_faces, tag_boundary_nodes

# helper array type for cut cell and hybrid meshes
include("named_array_partition.jl")
export NamedArrayPartition

include("hybrid_meshes.jl")
export num_faces, num_vertices, HybridMeshExample

include("physical_frame_basis.jl")
include("cut_cell_meshes.jl")
export PhysicalFrame, equi_nodes

include("state_redistribution.jl")
export StateRedistribution, apply!

include("nonconforming.jl")
export num_mortars_per_face, NonConformingQuadMeshExample

# uniform meshes + face vertex orderings
include("mesh/simple_meshes.jl")
export readGmsh2D, uniform_mesh
export readGmsh2D_v4, MeshImportOptions

# Plots.jl recipes for meshes
include("mesh/vtk_helper.jl")
include("mesh/mesh_visualization.jl")
export VertexMeshPlotter, MeshPlotter, Meshdata_to_vtk, type_to_vtk

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


using Requires: @require

function __init__()
  @require Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80" begin
    include("TriangulatePlots.jl")
  end
end


end # module
