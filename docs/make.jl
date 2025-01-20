using Documenter
using DocumenterInterLinks
using CamiDiff

links = InterLinks(
    "CamiMath" => "https://walra356.github.io/CamiMath.jl/stable/"
)

makedocs(;
    modules=[CamiDiff],
    plugins=[links],
    authors="= <walra356@planet.nl> and contributors",
    sitename="CamiDiff.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        size_threshold_warn = 250000,
        size_threshold = 300000,
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Manual" => "manual.md"
    ],
)

deploydocs(;
    repo="github.com/walra356/CamiDiff.jl.git",
    devbranch = "main"
)