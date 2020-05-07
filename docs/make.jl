using Documenter
using LaserTypes

makedocs(
    sitename = "LaserTypes",
    format = Documenter.HTML(),
    modules = [LaserTypes]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
