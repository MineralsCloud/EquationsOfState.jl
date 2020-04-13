using EquationsOfState
using Documenter

makedocs(;
    modules=[EquationsOfState],
    authors="Qi Zhang <singularitti@outlook.com>",
    repo="https://github.com/MineralsCloud/EquationsOfState.jl/blob/{commit}{path}#L{line}",
    sitename="EquationsOfState.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://MineralsCloud.github.io/EquationsOfState.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Installation" => "Installation.md",
        "Manual" => [
            "Collections" => "Collections.md",
            "Nonlinear fitting" => "NonlinearFitting.md",
            "Find volume" => "Find.md",
            "Portability" => "Portability.md",
            "Interoperability" => "Python.md",
        ],
        "FAQ" => "FAQ.md",
    ],
)

deploydocs(;
    repo="github.com/MineralsCloud/EquationsOfState.jl",
)
