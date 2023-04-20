using Documenter, BiologicalOscillations

include("pages.jl")

makedocs(
    sitename="BiologicalOscillations.jl",
    authors="Franco Tavella",
    doctest=false,
    pages=pages,
    )

deploydocs(repo = "github.com/YangLab-um/BiologicalOscillations.jl.git",)