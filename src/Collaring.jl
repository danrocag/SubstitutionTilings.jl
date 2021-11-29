module Collaring

using ...CoreDefs
using StructEquality
using LinearAlgebra

export collar_in

function collar_in end

function collar_class(tiling, g)
    return inv(g)*collar_in(tiling, g)
end

function center_label(patch)
    return Dict(patch)[id(patch[1][1])]
end
function center_tile(patch)
    e = id(patch[1][1])
    return (e, Dict(patch)[e])
end

function accessible_subst(S :: SubSystem{G,D,L}, initial_collar) where {G,D,L}
    collars = [initial_collar]
    visited = 0
    sub = Dict{Int,Vector{Tuple{G, Int}}}([])
    while length(collars) > visited
        new_visited = length(collars)
        for i in (visited+1):length(collars)
            sub[i] = []
            patch = substitute(S, collars[i], 1)
            for t in substitute(S, [center_tile(collars[i])], 1)
                class = collar_class(patch, t[1])
                search = findfirst(c -> issetequal(class, c), collars)
                if !isnothing(search)
                    push!(sub[i], (t[1], search))
                else
                    push!(collars, class)
                    push!(sub[i], (t[1], length(collars)))
                end
            end
        end
        visited = new_visited
    end
    return (collars, SubSystem(sub, S.λ))
end

function frequency(S, initial_collar, patch, depth)
    (collars, Sc) = accessible_subst(S, initial_collar)
    n = length(collars)
    println(n)
    A_tr = transition_matrix(Sc, 1:n)
    (eigenvalues, eigenvectors) = eigen(A_tr)
    λ_PF = eigenvalues[n]
    v_PF = eigenvectors[n]

    

    e = id(patch[1][1])

    patches_c = Base.product(([(t[1], l) for l in filter(c -> center_label(c) == t[2], collars)] for t in patch)...)

    freq = 0//1
    for label in 1:n
        domain = substitute(Sc, [(e, label)], depth-1)
        forced_uncollared_domain = substitute(S, collars[label], depth-1)
        forced_domain = Dict(filter(!isnothing, map(t -> (t[1], recognize_collar(forced_uncollared_domain, t[1])), forced_uncollared_domain)))
        for tile in domain
            if center_tile(tile[2]) == center_ptile
                translated_patches = Ref(tile[1]) .* patches_c
                for patch_c in translated_patches_c
                    if Dict(patch_c) ⊆ forced_domain
                        freq += v_PF[label]/λ_PF^(depth-1)
                    end
                end
            end
        end
    end
    return freq
end

end