using SubstitutionTilings
using SubstitutionTilings.Nilpotent

using Test
using LinearAlgebra

bhp = @time Dict(substitute(bhp_subst(), [a()], 3))

center_a = Dict([
    a(-2,-2,-2),
    a(0,-2,-2),
    a(2,-2,-2),
    a(-2,0,-2),
    a(0,0,-2),
    a(2,0,-2),
    a(-2,2,-2),
    a(0,2,-2),
    a(2,2,-2),
    b(-2,-2,0),
    b(0,-2,0),
    b(2,-2,0),
    b(-2,0,0),
    a(0,0,0),
    b(2,0,0),
    b(-2,2,0),
    b(0,2,0),
    b(2,2,0),
    a(-2,-2,2),
    a(0,-2,2),
    a(2,-2,2),
    a(-2,0,2),
    a(0,0,2),
    a(2,0,2),
    a(-2,2,2),
    a(0,2,2),
    a(2,2,2),
])
@time empirical_frequency(center_a, bhp)

patch = Dict([
    b(0,-2,0),
    a(0,0,0),
    b(0,2,0),
])
@time float(empirical_frequency(patch, bhp))


(collars, S) = total_collaring(bhp_subst(), center_a)
# 263 tiles

@time frequency(bhp_subst(), center_a, patch, 2)