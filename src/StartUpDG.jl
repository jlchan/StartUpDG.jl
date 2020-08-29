"""
Module StartUpDG

Module to aid in setting up reference operators, meshes, and geometric terms
"""

module StartupDG

# non-DG modules
using LinearAlgebra # for diagm, identity matrix I
using UnPack # for easy setting/getting in mutable structs
import VectorizedRoutines.Matlab.meshgrid
export meshgrid

using NodesAndModes # for basis functions
using UniformMeshes # for face vertex orderings

# Convenience routines for identity matrices.
export eye
eye(n) = diagm(ones(n))

# initialization of mesh/reference element data
export RefElemData
export init_reference_interval
export init_reference_tri, init_reference_quad
export init_reference_hex
include("ref_elem_functions.jl")

# spatial assembly routines
export MeshData
export init_mesh
export connect_mesh, build_node_maps, geometric_factors
export build_periodic_boundary_maps, build_periodic_boundary_maps!
include("./mesh_functions.jl")
include("./geometric_mapping_functions.jl")
include("./node_connectivity_functions.jl")

# uniform meshes included
export uniform_quad_mesh, uniform_tri_mesh
export uniform_hex_mesh
export tri_face_vertices, quad_face_vertices
export hex_face_vertices
include("uniform_meshes.jl")


end
