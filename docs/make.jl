using Documenter, EquationsOfState

makedocs(;
    modules=[EquationsOfState],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/singularitti/EquationsOfState.jl/blob/{commit}{path}#L{line}",
    sitename="EquationsOfState.jl",
    authors="Qi Zhang <singularitti@outlook.com>",
    assets=String[],
)

deploydocs(;
    repo="github.com/singularitti/EquationsOfState.jl",
)
