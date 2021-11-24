module Collaring

using ...CoreDefs
using StructEquality
using Graphs

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
                if class ∉ collars
                    push!(collars, class)
                end
                push!(sub[i], (t[1], findfirst(isequal(class), collars)) )
            end
        end
        visited = new_visited
    end
    return (collars, SubSystem(sub, S.λ))
end

end