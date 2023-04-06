using Test
using Suppressor
using LinearAlgebra
using StartUpDG

include("write_vtk_tests.jl")
include("named_array_partition_tests.jl")
include("triangulate_tests.jl")
include("reference_elem_tests.jl")
include("sbp_tests.jl")
include("MeshData_tests.jl")
include("boundary_util_tests.jl")
include("hybrid_mesh_tests.jl")
include("noncon_mesh_tests.jl")
include("cut_mesh_tests.jl")
include("misc_tests.jl")
include("gmsh_parse_tests.jl")
include("hohqmesh_tests.jl")
