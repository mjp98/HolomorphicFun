using HolomorphicFun
using Documenter

DocMeta.setdocmeta!(HolomorphicFun, :DocTestSetup, :(using HolomorphicFun); recursive=true)

makedocs(;
    modules=[HolomorphicFun],
    authors="Matthew Priddin and contributors",
    repo="https://github.com/mjp98/HolomorphicFun.jl/blob/{commit}{path}#{line}",
    sitename="HolomorphicFun.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://mjp98.github.io/HolomorphicFun.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/mjp98/HolomorphicFun.jl",
    devbranch="main",
)
