using EquationsOfState
using Documenter

DocMeta.setdocmeta!(EquationsOfState, :DocTestSetup, :(using EquationsOfState); recursive=true)

makedocs(;
    modules=[EquationsOfState],
    authors="Qi Zhang <singularitti@outlook.com>",
    repo="https://github.com/MineralsCloud/EquationsOfState.jl/blob/{commit}{path}#{line}",
    sitename="EquationsOfState.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://MineralsCloud.github.io/EquationsOfState.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/MineralsCloud/EquationsOfState.jl",
)
