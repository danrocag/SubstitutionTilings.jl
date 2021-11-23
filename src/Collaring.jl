module Collaring

using ...CoreDefs
using StructEquality
using Graphs

@def_structequal struct Collar{G, L}
    label :: L
    patch :: Tiling{G, L} # assumed to be a canonical patch with center label "label"
end

function recognize(tiling, collars)
    tiling_dict = Dict(tiling)
    result = []
    for tile in tiling
        for collar in collars
            if all(t -> (t[1] => t[2]) in tiling, t[1]*collar.patch)
                push!(result, collar)
                break
            end
        end
    end
    return result
end

struct UnrecognizedCollar <: Exception
end
function recognize(tiling, collars, tile)
    tiling_dict = Dict(tiling)
    for collar in collars
        if all(t -> (t[1] => t[2]) in tiling_dict, t[1]*collar.patch)
            return collar
        end
    end
    throw(UnrecognizedCollar)
end

function collared_subst(S, collars)
    sub = Dict([])
    for collar in collars
        image = substitute(S, collar.patch, 1)
        for t in S.sub[collar.label]
            return (t[1], recognize(tiling, collars, t))
        end
    end
    return SubSystem(sub, S.λ)
end


function sub_graph(S, collars)
    n = length(collars)
    graph = SimpleDiGraph(n)
    for i=1:n
        for t in S.sub[collars[i]]
            j = first(collars, t[2])
            add_edge!(graph, i, j)
        end
    end
    return graph
end

end