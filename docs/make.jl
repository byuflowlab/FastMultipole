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
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/byuflowlab/FLOWFMM",
    devbranch="main",
)
