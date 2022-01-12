module Collaring

using ...CoreDefs
using StructEquality
using LinearAlgebra

export collar_in, is_interior

function collar_in end
function is_interior end

function collar_class(tiling, g)
    return inv(g)*collar_in(tiling, g)
end

function center_label(patch)
    e = id(typeof(collect(keys(patch))[1]))
    return patch[e]
end
function center_tile(patch)
    e = id(typeof(collect(keys(patch))[1] ))
    return e => patch[e]
end

function accessible_subst(S :: SubSystem{G,D,L}, initial_collar) where {G,D,L}
    collars = [initial_collar]
    visited = 0
    sub = Dict{Int,Vector{Pair{G, Int}}}([])
    while length(collars) > visited
        new_visited = length(collars)
        for i in (visited+1):length(collars)
            sub[i] = []
            patch = substitute(S, collars[i], 1)
            for t in substitute(S, [center_tile(collars[i])], 1)
                class = collar_class(patch, t[1])
                search = findfirst(c -> issetequal(class, c), collars)
                if !isnothing(search)
                    push!(sub[i], t[1] => search)
                else
                    push!(collars, class)
                    push!(sub[i], t[1] => length(collars))
                end
            end
        end
        visited = new_visited
    end
    return (collars, SubSystem(sub, S.λ))
end

function product_patch_subset(patch, tiling)
    for (g, collars) in patch
        for collar in collars
            if g*collar ⊊ tiling
                return false
            end 
        end
    end
    return true
end

function frequency(S :: SubSystem{G, D, L}, initial_collar, patch, depth) where {G, D, L}
    (collars, Sc) = accessible_subst(S, initial_collar)
    n = length(collars)
    A_tr = transition_matrix(Sc, 1:n)
    (eigenvalues, eigenvectors) = eigen(A_tr)
    λ_PF = eigenvalues[n]
    v_PF = eigenvectors[:,n]/sum(eigenvectors[:,n])

    ptiles = unique([collar[id(G)] for collar in collars])
    collars_of_ptile = Dict([ptile => Int[] for ptile in ptiles])
    for i = 1:n
        push!(collars_of_ptile[collars[i][id(G)]], i)
    end


    patch_c = Dict([])
    for (g,l) in patch
        if is_interior(patch, g)
            patch_c[g] = (:interior, l)
        else
            patch_c[g] = (:exterior, collars_of_ptile[l])
        end
    end


    freq = 0.0
    for label in 1:n
        domain = substitute(Sc, [id(G) => label], depth-1)
        forced_uncollared_domain = substitute(S, collars[label], depth-1)

        for tile in domain
            translated_patch_c = tile[1] * patch_c
            is_subset = true
            for (g, (kind, detect)) in translated_patch_c
                if kind == :interior
                    if (g => detect) ∉ forced_uncollared_domain
                        is_subset = false
                        break
                    end
                else
                    if !any(c -> c ⊆ forced_uncollared_domain, Ref(g).*collars[detect])
                        is_subset = false
                        break
                    end
                end
            end
            if is_subset
                freq += v_PF[label]/λ_PF^(depth-1)
            end
        end
    end
    return freq
end

end