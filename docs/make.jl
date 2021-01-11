push!(LOAD_PATH,"src/") # necessary for unregistered packages

using Documenter
using StartUpDG

makedocs(
    sitename = "StartUpDG.jl",
    repo = "https://github.com/jlchan/StartUpDG.jl",
    modules=[StartUpDG],
    pages = [
        "Home" => "index.md",
        "Reference elements" => "refelem_usage.md",
        "Meshes" => "refelem_usage.md",
        "Time-stepping" => "tstep_usage.md",
        "Authors" => "authors.md",
        "Reference" => "refs.md"
    ]
)

deploydocs(
    repo = "github.com/jlchan/StartUpDG.jl",
)
