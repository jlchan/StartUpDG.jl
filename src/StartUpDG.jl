"""
Module StartUpDG

Module to aid in setting up reference operators, meshes, and geometric terms
"""

module StartUpDG

import UnPack: @unpack
export @unpack 

using LinearAlgebra # for diagm, identity matrix I
export eye # Convenience routine for identity matrices.
eye(n) = diagm(ones(n))

import VectorizedRoutines.Matlab.meshgrid
export meshgrid

using NodesAndModes # for basis functions

using UnPack # for setting/getting values in RefElemData and MeshData

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

# submodule for explicit time-stepping included for convenience
export ExplicitTimestepUtils
module ExplicitTimestepUtils
import ..@unpack
export bcopy!, bmult
export ck45, dp56 # carpenter/kennedy and dormand/prince
export PIparams, init_PI_controller, compute_adaptive_dt # for embedded RK methods like dp56()
include("explicit_timestep_utils.jl")
end

end
