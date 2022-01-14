using Documenter
using SubstitutionTilings

makedocs(
    sitename = "SubstitutionTilings",
    format = Documenter.HTML(prettyurls = false),
    modules = [SubstitutionTilings],
)


deploydocs(
    repo = "github.com/danrocag/SubstitutionTilings.jl.git"
)
