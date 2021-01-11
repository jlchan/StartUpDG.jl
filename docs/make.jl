push!(LOAD_PATH,"../src/") # necessary for unregistered packages

using Documenter
using StartUpDG

makedocs(
    sitename = "StartUpDG.jl",
    repo = "https://github.com/jlchan/StartUpDG.jl",
    modules=[StartUpDG],
    pages = [
        "Home" => "index.md",
        "Data structures" => "data_structures.md",
        "Reference elements" => "refelem_usage.md",
        "Meshes" => "mesh_usage.md",
        "Time-stepping" => "tstep_usage.md",
        "Authors" => "authors.md",
        "Reference" => "reference.md"
    ]
)

deploydocs(
    repo = "github.com/jlchan/StartUpDG.jl",
)
