"""
Module StartUpDG

Module to aid in setting up reference operators, meshes, and geometric terms
"""

module StartUpDG

# non-DG modules
using LinearAlgebra # for diagm, identity matrix I
using UnPack # for setting/getting values in RefElemData and MeshData
import VectorizedRoutines.Matlab.meshgrid
export meshgrid

using NodesAndModes # for basis functions

# Convenience routines for identity matrices.
export eye
eye(n) = diagm(ones(n))

# containers for reference element and physical data types.
# these are optional and provided for convenience.
export RefElemData
export MeshData
include("DG_types.jl")

# initialization of mesh/reference element data
export init_reference_interval
export init_reference_tri, init_reference_quad
export init_reference_hex
include("ref_elem_functions.jl")

# spatial assembly routines
export connect_mesh
export build_node_maps, geometric_factors
export build_periodic_boundary_maps, build_periodic_boundary_maps!
export init_DG_mesh
include("./mesh_functions.jl")
include("./geometric_mapping_functions.jl")
include("./node_connectivity_functions.jl")

# uniform meshes + face vertex orderings are included
export readGmsh2D
export uniform_quad_mesh, uniform_tri_mesh, uniform_hex_mesh
export tri_face_vertices, quad_face_vertices, hex_face_vertices
include("mesh_utils.jl")

end
