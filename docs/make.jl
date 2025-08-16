using Documenter
using DiffEqBase

DocMeta.setdocmeta!(DiffEqBase, :DocTestSetup, :(using DiffEqBase); recursive = true)

makedocs(;
    modules = [DiffEqBase],
    sitename = "DiffEqBase.jl",
    authors = "Christopher Rackauckas et al.",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://docs.sciml.ai/DiffEqBase/stable/",
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
        "API Reference" => [
            "Callbacks" => "api/callbacks.md",
            "Integrator Interface" => "api/integrator.md",
            "Statistics" => "api/statistics.md",
            "Utilities" => "api/utilities.md",
            "Tableaus" => "api/tableaus.md",
            "Internal Functions" => "api/internal.md",
        ],
        "Interfaces" => [
            "Problem Interface" => "interfaces/problems.md",
            "Solution Interface" => "interfaces/solutions.md", 
            "Algorithm Interface" => "interfaces/algorithms.md",
            "Function Interface" => "interfaces/functions.md",
        ],
    ],
)

deploydocs(;
    repo = "github.com/SciML/DiffEqBase.jl",
    devbranch = "master",
    push_preview = true,
)