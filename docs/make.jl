using Documenter
using LaserTypes
using Literate

# create literate versions of the source files
filepath = joinpath(@__DIR__, "..", "src")
files = joinpath.(filepath, readdir(filepath))

Literate.markdown.(files, joinpath(@__DIR__, "src", "generated"); documenter = false)

makedocs(
    sitename = "LaserTypes",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
    pages = [
        "index.md",
        "Laser profiles" => [
            "Spatial profiles" => [
                "Gaussian profile" => "gauss.md",
                "Laguerre-Gauss profile" => "laguerre-gauss.md"
            ],
            "temporal-profiles.md"
        ],
        "Source code" => joinpath.("generated",
                [
                "LaserTypes.md",
                "envelopes.md",
                "electricfield.md",
                "magneticfield.md",
                "potential.md",
                "gauss.md",
                "laguerre-gauss.md"
            ]
        ),
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
