using Documenter

makedocs(
    warnonly = :cross_references,
    sitename = "MRI.jl",
    format = Documenter.HTML(prettyurls = false, size_threshold_ignore = ["application-mri-saturation.md"]),
    pages = [
        "Introduction" => "index.md",
        "MRI: saturation problem" => "application-mri-saturation.md",
    ]
)

deploydocs(
    repo = "github.com/control-toolbox/MRI.jl.git",
    devbranch = "main"
)