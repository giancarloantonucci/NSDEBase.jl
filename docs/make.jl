using Documenter
using NSDEBase

PAGES = ["Home" => "index.md", "Guided tour" => "guide.md", "API" => "api.md"]

makedocs(;
    root=@__DIR__,
    sitename="NSDEBase.jl",
    format=Documenter.HTML(),
    modules=[NSDEBase],
    pages=PAGES,
    checkdocs=:exports, # every export must carry a docstring, or the build fails
    authors="Giancarlo A. Antonucci <giancarlo.antonucci@icloud.com>",
)

deploydocs(; repo="github.com/giancarloantonucci/NSDEBase.jl.git")
