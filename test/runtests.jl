using Test
using Suppressor
using LinearAlgebra

# we need to load this first before StartUpDG for Julia pre-v1.9 
# since Requires.jl needs us to load this first.
using SummationByPartsOperators

using StartUpDG

include("write_vtk_tests.jl")
include("triangulate_tests.jl")
include("reference_elem_tests.jl")
include("multidim_sbp_tests.jl")
include("sparse_SBP_operator_tests.jl")
include("SummationByPartsOperatorsExt_tests.jl")
include("MeshData_tests.jl")
include("MeshData_wedge_pyr_tests.jl")
include("boundary_util_tests.jl")
include("hybrid_mesh_tests.jl")
include("noncon_mesh_tests.jl")
include("cut_mesh_tests.jl")
include("misc_tests.jl")
include("gmsh_parse_tests.jl")
include("hohqmesh_tests.jl")
