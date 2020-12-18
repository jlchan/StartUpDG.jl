"""
Module StartUpDG

Module to aid in setting up reference operators, meshes, and geometric terms
"""

module StartUpDG

using NodesAndModes # for basis functions
export Line, Tri, Quad, Hex # element types from NodesAndModes

using UnPack # for setting/getting values in RefElemData and MeshData
using LinearAlgebra # for diagm, identity matrix I

import UnPack: @unpack
#import VectorizedRoutines.Matlab.meshgrid
import NodesAndModes.meshgrid

export @unpack
export meshgrid

# convenience containers for reference element and physical data types.
export RefElemData
export MeshData
include("DG_types.jl")

# initialization of mesh/reference element data
export init_reference_elem
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
export uniform_mesh
export face_vertices
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
