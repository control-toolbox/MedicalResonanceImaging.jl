using Documenter

mkpath("./docs/src/assets")
cp("./docs/Manifest.toml", "./docs/src/assets/Manifest.toml"; force=true)
cp("./docs/Project.toml", "./docs/src/assets/Project.toml"; force=true)

repo_url = "github.com/control-toolbox/MedicalResonanceImaging.jl"

makedocs(;
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
)

deploydocs(; repo=repo_url * ".git", devbranch="main")
