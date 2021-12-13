module CoreDefs

export SubSystem, Tiling, rational_to_float, substitute, check_subset, empirical_frequency, dilate, id, draw, embed_aff, Tiling, SetTiling, embed_field_elem, GroupElem, EGroupElem, vertices, in_border
export transition_matrix

using Luxor
using Dictionaries

import Base: *



"""
    GroupElem

A subtype of `GroupElem` is the type of elements of some group.
It should implement the function `id`, `*` and `inv`.
"""
abstract type GroupElem end

"""
   EGroupElem

A subtype of `GroupElem is the type of elements of some group which can be embedded in `Aff(2)`.
Thus, it should additionally implement the function `embed_aff`
"""
abstract type EGroupElem <: GroupElem end


function dilate end
"""
    embed_aff(g)

`embed_aff` should return the embedding of `g` as an element of `Aff(2)`.
This should be a vector length `6` with the three columns of the augmented matrix of the transformation, ignoring its last row.
For example, if `g` should be embedded as a reflection along the x axis followed by a translation by `(1,0)`,
the output should be `[1., 0., 0., -1., 1., 0.]`.
This is the format Luxor uses to represent affine transformations.
"""
function embed_aff end

"""
    draw(t :: T, s :: Symbol)

`T` should be a type of prototiles and `s` should be an action Luxor recognizes when drawing shapes (such as `:fill` or `:stroke`).
`draw(t,s)` should then draw some geometric primitives in a Luxor context using the action `s`.
"""
function draw end

"""
    id(g)


Returns the identity element of the group `g` belongs to.
"""
function id end


function Base.:*(g :: G, x :: Pair) where {G<:GroupElem}
    return g*x[1] => x[2]
end
function Base.:*(g :: G, X :: AbstractVector) where {G<:GroupElem}
    return map(x -> g*x, X)
end
function Base.:*(g :: G, X :: AbstractSet) where {G<:GroupElem}
    return map(x -> g*x, X)
end

function dilate(λ, x :: Pair)
    return dilate(λ,x[1]) => x[2]
end
function dilate(λ, X :: AbstractVector)
    return map(x -> dilate(λ,x), X)
end
function dilate(λ, X :: AbstractSet)
    return map(x -> dilate(λ,x), X)
end

"""
    draw(t, sc, hue, action)


`t` should be a tile: that is, a Pair consisting of some euclidean group element and some prototile.
`draw` draws the tile with the set scale, hue and Luxor action.
"""
function draw(t :: Pair{<:EGroupElem, T}, sc, hue, action :: Symbol) where T
    origin()
    scale(sc)
    sethue(hue)
    transform(embed_aff(t[1]))
    draw(t[2], action)
end


function vertices(t :: Pair{EGroupElem, T}) where T
    return t[1]*vertices(t[2])
end
function in_border(x, t :: Pair{EGroupElem, T}) where T
    return in_border(inv(t[1])*x, t[2])
end


"""
    SubSystem

The description of an inflation system.
`G` should be a subset of `EGroupElem` and `dilate(λ' :: D, g :: G)` should be implemented.
"""
struct SubSystem{G, D, L}
    sub::Dict{L, Vector{Pair{G, L}}}
    λ :: D
end
VectorTiling{G,L} = AbstractVector{Pair{G, L}}
DictTiling{G,L} = AbstractDict{G, L}


function rational_to_float(x)
    return Int(numerator(x))/Int(denominator(x))
end
function embed_field_elem(map_coeffs, map_gen, x)
    deg = degree(parent(x))
    result = map_coeffs(coeff(x, 0))
    for i in 1:(deg-1)
        result += map_coeffs(coeff(x, i))*map_gen^i
    end
    return result
end

"""
    substitute(S, tiling, n, [in_bounds, window])

Applies `n` inflation steps of the tiling:
each step dilates the tiling and then substitutes each tile according to the substitution rule of S.
This version computes the substitution depth first: this is more efficient unless the substitution rule has lots of overlaps
(see `substitute_bf`)

Optionally, one can filter which tiles are computed by providing `in_bounds` and `window`.
If provided, `in_bounds` should be a test function `in_bounds(tile, n, window)`
which recieves a tile, a depth and the supplied window,
and outputs whether the tile should be substituted further or discarded.
This can be used to greatly reduce compute times when rendering images of tilings.
"""
function substitute(S, tiling, n, in_bounds = nothing, window=nothing)
    result = typeof(tiling)()
    for tile in tiling
        substitute_df_inner!(S, n, tile, result, in_bounds, window)
    end
    return result
end
function substitute_df_inner!(S, n, tile, result, in_bounds = nothing, window=nothing)
    if isnothing(in_bounds) || in_bounds(tile, n, window)
        if n == 0
            push!(result, tile)
            return
        else
            for new_tile in S.sub[tile[2]]
                substitute_df_inner!(S, n-1, (dilate(S.λ, tile[1]) * new_tile), result, in_bounds, window)
            end
        end
    end
end

# set version of the substitute function: use when the substitution rule has overlaps
function substitute_bf(S, tiling, n=1) where {G, L,D}
    if n == 0
        return tiling
    else
        step = mapreduce(tile -> (dilate(S.λ, tile[1]) * S.sub[tile[2]]), union, tiling)
        return substitute(S, Set(step), n-1)
    end
end

"""
    center_patch

Translates the patch so that the first tile is located at the identity.
"""
function center_patch(patch)
    centered = inv(patch[1][1])*patch
    @assert centered[1][1] == id(patch[1][1])
    return centered
end

"""
    empirical_frequency

Computes the empirical relative frequency of `patch` in the `tiling`:
that is, computes how many translates of `patch` are subsets of `tiling`
divided by the total amount of tiles in `tiling`.
"""

function empirical_frequency(patch :: AbstractArray, tiling :: Dictionary)
    return empirical_frequency(Dictionary(Dict(patch)), tiling)
end
function empirical_frequency(patch :: Dictionary, tiling :: Dictionary)
    freq = 0//1
    n = 0

    origin_ptile = patch[id(collect(keys(patch))[1])]
    for g in keys(tiling)
        n += 1
        if origin_ptile == tiling[g]
            translated_patch = g*patch
            println(translated_patch)
            if all(g -> isassigned(tiling,g) && tiling[g] == translated_patch[g], keys(translated_patch))
                freq += 1
            end
        end
    end
    freq = freq//n
    return freq
end


struct InconclusiveSubsetError <: Exception
end
Base.showerror(io::IO, e::InconclusiveSubsetError) = print(io, "The given patch is neither a subset nor incompatible")
function check_subset(patch, tiling, incompatible)
    subset = true
    collision = false
    for x in patch
        x_found = false
        for t in tiling
            if x == t
                println("Equal:",x,",",t)
                x_found = true
                break
            end
            # We check if x and t are incompatible
            if incompatible(t[1], x[1])
                println("Incompatible:",x,",",t)
                return false
            end
        end
        if !x_found
            subset = false
        end
    end

    if subset
        return true
    else
        throw(InconclusiveSubsetError)
    end
end

function transition_matrix(S, order)
    n = length(order)
    return [count(t -> t[2]==order[j],S.sub[order[i]]) for i=1:n, j=1:n ]
end


end