using Documenter
using DocumenterInterLinks

#
links = InterLinks(
    "Tutorials" => (
        "https://control-toolbox.org/Tutorials.jl/stable/",
        "https://control-toolbox.org/Tutorials.jl/stable/objects.inv",
        joinpath(@__DIR__, "inventories", "Tutorials.toml"),
    ),
)

# For reproducibility
mkpath(joinpath(@__DIR__, "src", "assets"))
cp(
    joinpath(@__DIR__, "Manifest.toml"),
    joinpath(@__DIR__, "src", "assets", "Manifest.toml");
    force=true,
)
cp(
    joinpath(@__DIR__, "Project.toml"),
    joinpath(@__DIR__, "src", "assets", "Project.toml");
    force=true,
)

repo_url = "github.com/control-toolbox/MedicalResonanceImaging.jl"

# if draft is true, then the julia code from .md is not executed # debug
# to disable the draft mode in a specific markdown file, use the following:
#=
```@meta
Draft = false
```
=#
makedocs(;
    draft=false,
    warnonly=:cross_references,
    sitename="Medical Resonance Imaging",
    format=Documenter.HTML(;
        repolink="https://" * repo_url,
        prettyurls=false,
        size_threshold_ignore=["saturation.md", "bloch-equation.md"],
        assets=[
            asset("https://control-toolbox.org/assets/css/documentation.css"),
            asset("https://control-toolbox.org/assets/js/documentation.js"),
        ],
    ),
    pages=[
        "Introduction" => "index.md",
        "Bloch equation" => "bloch-equation.md",
        "Saturation problem" => "saturation.md",
    ],
    plugins=[links],
)

deploydocs(; repo=repo_url * ".git", devbranch="main", push_preview=true)
