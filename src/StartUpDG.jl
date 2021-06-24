module StartUpDG

using ConstructionBase  # for constructorof(...)
using LinearAlgebra     # for diagm, identity matrix I
using StaticArrays      # for SMatrix
using Setfield          # for "modifying" structs to modify node mappings (setproperties)

using Reexport 
@reexport using UnPack            # for getting values in RefElemData and MeshData
@reexport using NodesAndModes     # for basis functions

using NodesAndModes: meshgrid
using SparseArrays: sparse, droptol!
using MAT # read in .mat files
using RecipesBase, Colors 

# reference element utility functions
include("RefElemData.jl")
include("RefElemData_polynomial.jl")
include("RefElemData_SBP.jl")
include("ref_elem_utils.jl")
include("element_type_utils.jl")
export RefElemData, Polynomial, SBP
export DefaultSBPType, DGSEM, Hicken, Kubatko  # types for SBP node dispatch
export LobattoFaceNodes, LegendreFaceNodes # type parameters for SBP{Kubatko{...}}
export hybridized_SBP_operators, inverse_trace_constant

include("MeshData.jl")
export MeshData 

include("geometric_functions.jl")
export geometric_factors, estimate_h

# spatial connectivity routines
include("connectivity_functions.jl")
export make_periodic

# for tagging faces on boundaries
include("boundary_utils.jl")
export boundary_face_centroids, tag_boundary_faces

# uniform meshes + face vertex orderings are included
include("mesh/simple_meshes.jl")
export readGmsh2D, uniform_mesh

include("mesh/mesh_visualization.jl")
export VertexMeshPlotter, MeshPlotter

using Requires
using Printf
function __init__()                 
    @require Triangulate="f7e6ffb2-c36d-4f8f-a77e-16e897189344" begin
        using .Triangulate: TriangulateIO, triangulate
        include("mesh/triangulate_utils.jl")      
        export refine, triangulateIO_to_VXYEToV, get_node_boundary_tags
        export BoundaryTagPlotter
        include("mesh/triangulate_example_meshes.jl")
        export square_domain, rectangular_domain, square_hole_domain, scramjet
    end
end

# simple explicit time-stepping included for conveniencea
export ck45 # LSERK 45 
include("explicit_timestep_utils.jl")

end
