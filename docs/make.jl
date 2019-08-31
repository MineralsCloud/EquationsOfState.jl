using Documenter, EquationsOfState

makedocs(;
    modules=[EquationsOfState],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
        "Manual" => [
            "Collections" => "Collections.md",
            "Nonlinear fitting" => "NonlinearFitting.md"
        ]
    ],
    repo="https://github.com/MineralsCloud/EquationsOfState.jl/blob/{commit}{path}#L{line}",
    sitename="EquationsOfState.jl",
    authors="Qi Zhang <singularitti@outlook.com>",
    assets=String[],
)

deploydocs(;
    repo="github.com/MineralsCloud/EquationsOfState.jl.git"
)
