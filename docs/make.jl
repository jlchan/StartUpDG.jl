using Documenter
using StartUpDG

makedocs(
    sitename = "StartUpDG.jl",
    repo = "https://github.com/jlchan/StartUpDG.jl",
    modules=[StartUpDG],
    pages = [
        "Home" => "index.md",
        "Authors" => "authors.md"
    ]
)

deploydocs(
    repo = "github.com/jlchan/StartUpDG.jl",
)
