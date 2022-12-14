using Test
using Suppressor
using LinearAlgebra
using RecipesBase
using Triangulate # load before StartUpDG b/c of @require
using StartUpDG

include("triangulate_tests.jl")
# include("geometric_term_tests.jl")
include("reference_elem_tests.jl")
include("sbp_tests.jl")
include("MeshData_tests.jl")
include("boundary_util_tests.jl")
include("hybrid_mesh_tests.jl")
include("noncon_mesh_tests.jl")
include("cut_mesh_tests.jl")
include("misc_tests.jl")
include("gmsh_parse_tests.jl")
