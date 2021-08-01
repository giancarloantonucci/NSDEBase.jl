using Documenter
using NSDEBase

PAGES = ["Home" => "index.md"]

makedocs(;
    sitename = "NSDEBase",
    format = Documenter.HTML(),
    modules = [NSDEBase],
    pages = PAGES,
    authors = "Giancarlo A. Antonucci <giancarlo.antonucci@icloud.com>"
)

deploydocs(;
    repo = "https://github.com/antonuccig/NSDEBase.jl"
)
