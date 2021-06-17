using Documenter
using LaserTypes

makedocs(
    sitename = "LaserTypes",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
    pages = [
        "index.md",
        "setup.md",
        "Laser profiles" => [
            "Spatial profiles" => [
                "Gaussian profile" => "gauss.md",
                "Laguerre-Gauss profile" => "laguerre-gauss.md"
            ],
            "temporal-profiles.md"
        ],
        "covariant.md",
        "advanced.md"
    ],
    modules = [LaserTypes]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/SebastianM-C/LaserTypes.jl",
    push_preview = true
)
