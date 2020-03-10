using Documenter, SWAXS

makedocs(;
    modules=[SWAXS],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/jipq6175/SWAXS.jl/blob/{commit}{path}#L{line}",
    sitename="SWAXS.jl",
    authors="Yen-Lin Chen, Cornell University",
    assets=String[],
)

deploydocs(;
    repo="github.com/jipq6175/SWAXS.jl",
)
