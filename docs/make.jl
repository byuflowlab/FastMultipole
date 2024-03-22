using FLOWFMM
using Documenter

DocMeta.setdocmeta!(FLOWFMM, :DocTestSetup, :(using FLOWFMM); recursive=true)

makedocs(;
    modules=[FLOWFMM],
    authors="Ryan Anderson <rymanderson@gmail.com> and contributors",
    sitename="FLOWFMM.jl",
    format=Documenter.HTML(;
        canonical="https://flow.byu.edu/FLOWFMM",
        edit_link="main",
        assets=String[],
    ),
    pages=[
		"Intro" => "index.md",
        "Quick Start" => "tutorial.md",
        "Guided Examples" => "howto.md",
        "API Reference" => "reference.md",
        "Theory" => "theory.md"
    ],
)

deploydocs(;
    repo="github.com/byuflowlab/FLOWFMM",
    devbranch="main",
)
