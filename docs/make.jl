using EquationsOfState
using Documenter

DocMeta.setdocmeta!(EquationsOfState, :DocTestSetup, :(using EquationsOfState.Collections, Unitful, UnitfulAtomic); recursive=true)
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
        "Installation" => "installation.md",
        "Manual" => [
            "Collections" => "Collections.md",
            "Nonlinear fitting" => "NonlinearFitting.md",
            "Find volume" => "Find.md",
            "Portability" => "portability.md",
            "Interoperability" => "interoperability.md",
        ],
        "FAQ" => "faq.md",
    ],
)

deploydocs(;
    repo="github.com/MineralsCloud/EquationsOfState.jl",
)
