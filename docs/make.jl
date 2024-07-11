using Documenter

makedocs(
    warnonly = :cross_references,
    sitename = "MRI.jl",
    format = Documenter.HTML(
        prettyurls = false, 
        size_threshold_ignore = ["saturation.md", "bloch-equation.md"]),
    pages = [
        "Introduction"       => "index.md",
        "Bloch equation"     => "bloch-equation.md",
        "Saturation problem" => "saturation.md",
    ]
)

deploydocs(
    repo = "github.com/control-toolbox/MRI.jl.git",
    devbranch = "main"
)