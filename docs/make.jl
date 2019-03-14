using Documenter
using DiffEqBase

makedocs(
    sitename = "DiffEqBase",
    format = :html,
    modules = [DiffEqBase]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
