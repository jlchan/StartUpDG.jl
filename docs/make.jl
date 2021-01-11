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
        "DG mesh data" => "mesh_usage.md",
        "Timestepping" => "tstep_usage.md",
        "Reference" => "index_refs.md",
        "Authors" => "authors.md",
    ]
)

deploydocs(
    repo = "github.com/jlchan/StartUpDG.jl",
)
