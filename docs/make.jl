push!(LOAD_PATH,"../src/") # necessary for unregistered packages

using Documenter
using StartUpDG

makedocs(
    sitename = "StartUpDG.jl",
    repo = "https://github.com/jlchan/StartUpDG.jl",
    modules=[StartUpDG],
    pages = [
        "Home" => "index.md",        
        "Background and conventions" => "conventions.md",        
        "`RefElemData`" => "RefElemData.md",
        "`MeshData`" => "MeshData.md",
        "Example: computing DG derivatives" => "ex_dg_deriv.md",
        "Hybrid meshes" => "hybrid_mesh.md",
        "Timestepping" => "tstep_usage.md",
        "Reference" => "index_refs.md",
        "Authors" => "authors.md",
    ]
)

deploydocs(
    repo = "github.com/jlchan/StartUpDG.jl",
)
