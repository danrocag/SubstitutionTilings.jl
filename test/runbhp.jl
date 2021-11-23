using SubstitutionTilings
using SubstitutionTilings.Nilpotent

using Test

@time substitute(bhp_subst(), [a()], 3)
bhp = @time substitute(bhp_subst(), [a()], 4)

isolated_a = [
    a(-2, -2, -2),
    a(0, -2, -2),
    a(2, -2, -2),
    a(-2, -2, 0),
    b(0, -2, 0),
    a(2, -2, 0),
    a(-2, -2, 2),
    a(0, -2, 2),
    a(2, -2, 2),
]

@time empirical_frequency(isolated_a, bhp)