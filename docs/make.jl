using FastMultipole
using Documenter

DocMeta.setdocmeta!(FastMultipole, :DocTestSetup, :(using FastMultipole); recursive=true)

makedocs(;
    modules=[FastMultipole],
    authors="Ryan Anderson <rymanderson@gmail.com> and contributors",
    sitename="FastMultipole.jl",
    format=Documenter.HTML(;
        canonical="https://flow.byu.edu/FastMultipole",
        edit_link="main",
        assets=String[],
    ),
    pages=[
		"Introduction" => "index.md",
        "Quick Start" => "quickstart.md",
        "Guided Examples" => "guided_examples.md",
        "Reference" => "reference.md",
        "Theory" => "theory.md"
    ],
)

deploydocs(;
    repo="github.com/byuflowlab/FastMultipole",
    devbranch="main",
)
