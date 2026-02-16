using Documenter
using ERGMCount

DocMeta.setdocmeta!(ERGMCount, :DocTestSetup, :(using ERGMCount); recursive=true)

makedocs(
    sitename = "ERGMCount.jl",
    modules = [ERGMCount],
    authors = "Statistical Network Analysis with Julia",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://Statistical-network-analysis-with-Julia.github.io/ERGMCount.jl",
        edit_link = "main",
    ),
    repo = "https://github.com/Statistical-network-analysis-with-Julia/ERGMCount.jl/blob/{commit}{path}#{line}",
    pages = [
        "Home" => "index.md",
        "Getting Started" => "getting_started.md",
        "User Guide" => [
            "Reference Measures" => "guide/references.md",
            "Count Terms" => "guide/terms.md",
            "Estimation" => "guide/estimation.md",
        ],
        "API Reference" => [
            "Types" => "api/types.md",
            "Terms" => "api/terms.md",
            "Estimation" => "api/estimation.md",
        ],
    ],
    warnonly = [:missing_docs, :docs_block],
)

deploydocs(
    repo = "github.com/Statistical-network-analysis-with-Julia/ERGMCount.jl.git",
    devbranch = "main",
    versions = [
        "stable" => "dev",
        "dev" => "dev",
    ],
    push_preview = true,
)
