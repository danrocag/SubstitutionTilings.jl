using SubstitutionTilings
using SubstitutionTilings.Nilpotent
using SubstitutionTilings.Collaring

using Test
using LinearAlgebra

t = substitute(bhp_subst(), [a()], 2)
bhp = @time substitute(bhp_subst(), [a()], 3)

center_a = [
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
]
@time empirical_frequency(center_a, bhp, 14)


(collars, S) = Collaring.accessible_subst(bhp_subst(), center_a)
collars[1] == center_a
length(collars)
function support(patch)
    return map(t -> t[1], patch)
end
supports = unique((support.(collars)))

A_bhp = transition_matrix(S, 1:263)

(e, V) = eigen(A_bhp)
Rational.(A_bhp)*[1//263 for i=1:263] == [81//263 for i=1:263]

function forget_collar(i)
    return Dict(collars[i])[He(0,0,0)]
end

count(i -> forget_collar(i))

function recognize_collar(tiling, g)
    try
        return findfirst(isequal(Collaring.collar_class(tiling, g)), collars)
    catch e
        return nothing
    end
end

function collared_frequency(patch, depth)
    freq = 0//1

    center_ptile = Dict(patch)[He(0,0,0)]

    for label in 1:264
        domain = substitute(S, [(He(0,0,0), label)], depth-1)
        forced_uncollared_domain = substitute(bhp_subst(), collars[label], depth)
        forced_domain = Dict(filter(!isnothing, map(t -> (t[1], recognize_collar(forced_uncollared_domain, t[1])), forced_uncollared_domain)))
        for tile in domain
            if tile[2] == center_ptile
                translated_patch = tile[1] * patch
                if all(t -> (t[1] => t[2]) in forced_domain, translated_patch)
                    freq += 1//264//81^(depth-1)
                end
            end
        end
    end
    return freq
end

collared_frequency([(He(0,0,0), 2)],2)