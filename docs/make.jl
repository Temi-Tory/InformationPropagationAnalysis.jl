using Documenter
using InformationPropagationAnalysis

DocMeta.setdocmeta!(InformationPropagationAnalysis, :DocTestSetup,
                    :(using InformationPropagationAnalysis); recursive = true)

makedocs(
    sitename = "InformationPropagationAnalysis.jl",
    authors = "T. Ohiani",
    modules = [
        InformationPropagationAnalysis,
        InformationPropagationAnalysis.Input,
        InformationPropagationAnalysis.Diamonds,
        InformationPropagationAnalysis.Probability,
        InformationPropagationAnalysis.CriticalPath,
        InformationPropagationAnalysis.Flow,
    ],
    # the staging build has no git remote; name it explicitly so source links resolve
    repo = Remotes.GitHub("Temi-Tory", "InformationPropagationAnalysis.jl"),
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://Temi-Tory.github.io/InformationPropagationAnalysis.jl",
        edit_link = "main",
    ),
    pages = [
        "Home" => "index.md",
        "Toolkits" => [
            "Input" => "input.md",
            "Diamonds" => "diamonds.md",
            "Probability" => "probability.md",
            "CriticalPath" => "criticalpath.md",
            "Flow" => "flow.md",
        ],
        "Reproducibility" => "reproducibility.md",
        "API reference" => "api.md",
    ],
    # ProbabilityBoundsAnalysis cannot precompile on Julia 1.12, so doctests run
    # interpreted and are slow; keep them limited to the stable-output blocks.
    doctest = true,
    checkdocs = :none,
)

deploydocs(;
    repo = "github.com/Temi-Tory/InformationPropagationAnalysis.jl",
    devbranch = "main",
)
