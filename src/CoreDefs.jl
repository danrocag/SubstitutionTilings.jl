module CoreDefs

export SubSystem, substitute, check_subset, empirical_frequency, dilate, id, draw, embed_aff, DGroupElem
export collar_in, is_interior, frequency, total_collaring, UnrecognizedCollar, transition_matrix
export vertices, @collar_in_from_vertices
export embed_center
export autocorrelation


using Luxor
using LinearAlgebra


"""
    DGroupElem

A subtype of `DGroupElem` is the type of elements of some group
with a dilation.
It should implement the functions `id`, `*`, `inv` and `dilate`.
"""
abstract type DGroupElem end

"""
    dilate(λ, g)

`embed_aff` dilates an group element `g` by the coefficient `λ`.
"""
function dilate end


"""
    id(g)


Returns the identity element of the group `g` belongs to.
"""
function id end

function Base.:*(g :: G, x :: Pair) where {G<:DGroupElem}
    return g*x[1] => x[2]
end
function Base.:*(g :: G, X :: AbstractVector) where {G<:DGroupElem}
    return map(x -> g*x, X)
end
function Base.:*(g :: G, X :: AbstractSet) where {G<:DGroupElem}
    return map(x -> g*x, X)
end
function Base.:*(g :: G, X :: AbstractDict) where {G<:DGroupElem}
    return typeof(X)(g*x for x in X)
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
    embed_aff(g)

`embed_aff` should return the embedding of `g` as an element of `Aff(2)`.
This should be a vector length `6` with the three columns of the augmented matrix of the transformation, ignoring its last row.
For example, if `g` should be embedded as a reflection along the x axis followed by a translation by `(1,0)`,
the output should be `[1., 0., 0., -1., 1., 0.]`.
This is the format Luxor uses to represent affine transformations.
"""
function embed_aff end

"""
    draw(t :: T, action :: Symbol)

`T` should be a type of prototiles and `s` should be an action Luxor recognizes when drawing shapes (such as `:fill` or `:stroke`).
`draw(t,s)` should then draw some geometric primitives in a Luxor context using the action `s`.
"""
function draw end


"""
    draw(t :: Pair{<:DGroupElem, T}, sc, hue, action)


If `t` is a tile, `draw(t, sc, hue, action)` draws it using
`draw(t[1], action)` after applying the transformation `embed_aff(t[2])`
with scale `sc`, hue `hue` and action p `action`.

"""
function draw(t :: Pair{<:DGroupElem, T}, sc, hue, action :: Symbol) where T
    scale(sc)
    sethue(hue)
    transform(embed_aff(t[1]))
    draw(t[2], action)
end


"""
    embed_center(t :: Pair{<:DGroupElem, T})

If `t` is a tile with coordinates that can be embedded in the euclidean plane,
outputs the center of the tile as a Luxor point
"""

function embed_center(g :: DGroupElem)
    return Point(embed_aff(g)[5:6]...)
end

"""
    SubSystem{G, D, L}

The description of a substitution system with coordinates in `G`,
labels of type `L`
and dilation coefficients of type `D`.
`G` should be a subset of `DGroupElem` implementing `dilate(λ :: D, g :: G)`.
"""
struct SubSystem{G, D, L}
    sub::Dict{L, Vector{Pair{G, L}}}
    λ :: D
end
VectorTiling{G,L} = AbstractVector{Pair{G, L}}
DictTiling{G,L} = AbstractDict{G, L}


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

"""
    substitute_bf(S, tiling,n=1) 

Set version of the substitute function: use when the substitution rule has overlaps
"""
function substitute_bf(S, tiling, n=1)
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
    empirical_frequency(patch, tiling)

Computes the empirical relative frequency of `patch` in the `tiling`:
that is, computes how many translates of `patch` are subsets of `tiling`
divided by the total amount of tiles in `tiling`.
"""
function empirical_frequency(patch :: Dict{G, L}, tiling :: Dict{G, L}) where {G<:DGroupElem, L}
    freq = 0//1
    n = 0

    origin_ptile = patch[id(G)]
    for g in keys(tiling)
        n += 1
        if origin_ptile == tiling[g]
            translated_patch = g*patch
            if translated_patch ⊆ tiling
                freq += 1
            end
        end
    end
    freq = freq//n
    return freq
end
function empirical_frequency(patch :: AbstractArray, tiling :: Dict)
    return empirical_frequency((Dict(patch)), tiling)
end
function empirical_frequency(patch :: AbstractArray, tiling :: AbstractArray)
    freq = 0//1
    n = 0

    origin_ptile = patch[1]
    for (g,l) in tiling
        n += 1
        if (g => l) in tiling
            translated_patch = g*patch
            if translated_patch ⊆ tiling
                freq += 1
            end
        end
    end
    freq = freq//n
    return freq
end

"""
    transition_matrix(S, order)

Given a substitution system `S` and a list of labels `order`,
computes the transition matrix of `S` in that order
"""
function transition_matrix(S, order)
    n = length(order)
    return [count(t -> t[2]==order[i],S.sub[order[j]]) for i=1:n, j=1:n ]
end


struct UnrecognizedCollar <: Exception end
function collar_in end
function is_interior end
function vertices end

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

macro collar_in_from_vertices(G, L, zero_angle, full_circle)
    return quote
        function precollar_in(tiling :: Dict{$(esc(G)), $(esc(L))}, g :: $(esc(G)))
            gvertices = vertices(g => tiling[g])
            precollar =  []
            for t in tiling
                if !isempty(keys(vertices(t)) ∩ keys(gvertices))
                    push!(precollar,t)
                end
            end
            return Dict(precollar)
        end
        function CoreDefs.is_interior(tiling :: Dict{$(esc(G)), $(esc(L))}, g :: $(esc(G)))
            angle_sums = Dict(k => $zero_angle for k in keys(vertices(g => tiling[g])))
            for tile in tiling
                for (vertex,angle) in vertices(tile)
                    if vertex ∈ keys(angle_sums)
                        angle_sums[vertex] = angle_sums[vertex] .+ angle
                    end
                end
            end
        
            for (vertex, angle) in angle_sums
                if (angle ≠ $full_circle)
                    return false
                end
            end
            return true
        end
        function CoreDefs.collar_in(tiling :: Dict{$(esc(G)), $(esc(L))}, g :: $(esc(G)); check=false)
            precollar =  precollar_in(tiling , g)
            if check && !is_interior(precollar, g)
                throw(UnrecognizedCollar)
            else
                return precollar
            end
        end
    end
end

"""
    total_collaring(S, initial_collar)

Given a (primitive) substitution system `S` and some collar `initial_collar` (centered at the origin),
construct its total collaring.
The result is given as a tuple `(collars, Sc)`
where `collars` is the list of all collars that happen in `S`-legal tilings
and `Sc` is the corresponding substitution system,
where the labels are integers constituting indexes in `collars`.

Requires `collar_in` to be implemented for the given substitution system.
"""
function total_collaring(S :: SubSystem{G,D,L}, initial_collar) where {G,D,L}
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

"""
Given a (primitive) substitution system `S`,
some collar `initial_collar` and some patch `patch` (both centered at the origin) and some depth `depth`
calculates the frequency of `patch` in `S`.
`initial_collar` is used to construct a total collaring and doesn't matter as long as it's a legal collar.
`depth` determines how the level at which the frequency is calculated:
the result is only exact if `depth` is high enough, but lower values result in faster computation.

Requires `collar_in` to be implemented for the given substitution system.
"""
function frequency(S :: SubSystem{G, D, L}, initial_collar, patch, depth) where {G, D, L}
    (collars, Sc) = total_collaring(S, initial_collar)
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

# Radius should be forced by the chosen prototile depth!
"""
Given a (primitive) substitution system `S`,
some collar `initial_collar` and some patch `patch` (both centered at the origin) and some depth `depth`,
calculates the autocorrelation measure of a tiling for all differences between supertiles of depth `depth`.
The values are accurate as long as `depth` is high enough to force collars.
`weights(i :: L,g :: G,j :: L)` gives the weight used to compute the autocorrelation: by default it is constantly `1`

`collar_in` should be implemented.
"""
function autocorrelation(S :: SubSystem{G, D, L}, initial_collar, depth; weights=nothing) where {G, D, L}
    (collars, Sc) = total_collaring(S, initial_collar)
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

    measure = Dict{G, Float64}()
    if isnothing(weights)
        for label in 1:n
            domain = substitute(S, [id(G) => center_label(collars[label])], depth-1)
            forced_domain = substitute(S, collars[label], depth-1)

            for tile in domain
                translated_f_domain = inv(tile[1]) * forced_domain
                for t in translated_f_domain
                    g = t[1]
                    if g in Set(keys(measure))
                        measure[g] += v_PF[label]/λ_PF^(depth-1)
                    else
                        measure[g] = v_PF[label]/λ_PF^(depth-1)
                    end
                end
            end
        end
    else
        for label in 1:n
            domain = substitute(S, [id(G) => center_label(collars[label])], depth-1)
            forced_domain = substitute(S, collars[label], depth-1)

            for tile in domain
                translated_f_domain = inv(tile[1]) * forced_domain
                for t in translated_f_domain
                    g = t[1]
                    if g in Set(keys(measure))
                        measure[g] += v_PF[label]/λ_PF^(depth-1)*weights(tile[2], t)
                    else
                        measure[g] = v_PF[label]/λ_PF^(depth-1)*weights(tile[2], t)
                    end
                end
            end
        end
    end
    return measure
end

end