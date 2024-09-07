using Documenter

makedocs(;
    warnonly=:cross_references,
    sitename="Medical Resonance Imaging",
    format=Documenter.HTML(;
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

deploydocs(;
    repo="github.com/control-toolbox/medical_resonance_imaging.git", devbranch="main"
)
