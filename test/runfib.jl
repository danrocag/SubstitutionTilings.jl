# The Fibonacci tiling

using SubstitutionTilings
using SubstitutionTilings.NumFields
using StructEquality
using Luxor

@NumFields.simple_number_field_concrete Qτ [1,1] τ
Base.promote_rule(::Type{Qτ}, ::Type{<:Integer}) = Qτ

@def_structequal struct FibElem <: DGroupElem
    a :: Qτ
end

@enum FibTile begin
    A=1
    B=2
end

function Base.:*(x :: FibElem, y :: FibElem)
    return FibElem(x.a + y.a)
end

function SubstitutionTilings.dilate(λ :: Qτ, x :: FibElem)
    return FibElem(λ*x.a)
end

function Base.inv(x :: FibElem)
    return FibElem(-x.a)
end

function SubstitutionTilings.id(::Type{FibElem})
    return FibElem(0)
end

function SubstitutionTilings.embed_aff(g :: FibElem)
    fτ = (1+sqrt(5))/2
    return [1, 0, 0, 1, embed_field(float,fτ, g.a), 0]
end

function SubstitutionTilings.draw(ptile::FibTile, action)
    fτ = (1+sqrt(5))/2
    if ptile == A
        return Luxor.poly([
            Point(-fτ*0.9,-10),
            Point(-fτ*0.9,10),
            Point(fτ*0.9,10),
            Point(fτ*0.9,-10)
        ], close = true, action)
    else
        return Luxor.poly([
            Point(-0.9,-10),
            Point(-0.9,10),
            Point(0.9,10),
            Point(0.9,-10)
        ], close = true, action)
    end
end

fib = SubSystem(Dict(A => [FibElem(Qτ(-1)//2) => A, FibElem(τ//2) => B], B => [FibElem(0) => A]),τ)

fib_tiling = substitute(fib, Dict([FibElem(0) => A]), 3)

function color(tile :: Pair{FibElem,FibTile})
    return (tile[2] == A) ? 1 : 2
end

width = 800
height = 20
sc = 20
@png begin
    colors = ["#DD93FC", "#E7977A"]
    patch = [FibElem(-τ//2) => A, FibElem(τ//2) => A]
    tiling = substitute(fib, patch, 4)
    setline(1)

    for tile in tiling
        origin()
        draw(tile, sc, colors[color(tile)], :fill)
    end

    origin()
    scale(sc)
    sethue("black")
    Luxor.poly([
            Point(-0.05,-10),
            Point(-0.05,10),
            Point(0.05,10),
            Point(0.05,-10)
        ], close = true, :fill)
end width height

function SubstitutionTilings.is_interior(tiling :: Dict, t :: FibElem)
    return haskey(tiling, t) && (haskey(tiling, FibElem(t.a-τ)) || haskey(tiling, FibElem(t.a-(1+τ)//2))) && (haskey(tiling, FibElem(t.a+τ)) || haskey(tiling, FibElem(t.a+(1+τ)//2)))
end

function SubstitutionTilings.collar_in(tiling :: Dict, t :: FibElem)
    if !is_interior(tiling, t)
        throw(UnrecognizedCollar)
    end
    collar = Dict{FibElem, FibTile}([])
    shifts = ([0, -τ, -(1+τ)//2, τ, (1+τ)//2])
    for shift in shifts
        if haskey(tiling, FibElem(t.a+shift))
            collar[FibElem(t.a+shift)] = tiling[FibElem(t.a+shift)]
        end
    end
    return collar
end

initial_collar = collar_in(fib_tiling, FibElem(0))
(collars, fib_c) = total_collaring(fib, initial_collar)

collars

width = 800
height = 40
sc = 20
@png begin
    colors = ["#DD93FC", "#E7977A"]
    first_tile = [FibElem(0) => A]
    tiling = collars[4]
    setline(1)

    for tile in tiling
        origin()
        draw(tile, sc, colors[color(tile)], :fill)
    end
end width height

frequency(fib, initial_collar, collars[4], 4)
