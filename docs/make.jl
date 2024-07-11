using Documenter

makedocs(
    sitename = "MRI.jl",
    format = Documenter.HTML(prettyurls = false),
    pages = [
        "Introduction" => "index.md"
    ]
)

deploydocs(
    repo = "github.com/control-toolbox/MRI.jl.git",
    devbranch = "main"
)