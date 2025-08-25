using SubstitutionTilings
using SubstitutionTilings.CoreDefs
using StructEquality
using LinearAlgebra

using Luxor

@struct_hash_equal_isequal struct AT <: DGroupElem
    a::Int
end

@enum ATTile begin
    A=1
    B=2
end

function Base.:*(x :: AT, y :: AT)
    return AT(x.a + y.a)
end

function SubstitutionTilings.dilate(λ, x :: AT)
    return AT(λ*x.a)
end

function SubstitutionTilings.inv(x :: AT)
    return AT(-x.a)
end

function SubstitutionTilings.id(::Type{AT})
    return AT(0)
end

function SubstitutionTilings.is_interior(tiling :: Dict, t :: AT)
    return haskey(t.a, tiling) && haskey(t.a-1, tiling) && haskey(t.a+1, tiling)
end

function SubstitutionTilings.embed_aff(g :: AT)
    return [1, 0, 0, 1, g.a, 0]
end

function SubstitutionTilings.collar_in(tiling :: Dict, t :: AT)
    collar_shape = t*[AT(x) for x in [-1,0,1]]
    collar = []
    tiling_dict = tiling
    for s in collar_shape
        if !haskey(tiling_dict, s)
            throw(UnrecognizedCollar)
        end
        push!(collar, (s => tiling_dict[s]))
    end
    return Dict(collar)
end

at = SubSystem(Dict(
    A =>[
        (AT(0) => A),
        (AT(1) => A),
        (AT(2) => A),
        (AT(3) => A),
        (AT(4) => B)
    ], B => [
        (AT(0) => A),
        (AT(1) => B),
        (AT(2) => B),
        (AT(3) => B),
        (AT(4) => B)
    ]),5)

substitute(at, [(AT(0) => A)], 4)

(collars, S) = total_collaring(at, Dict([(AT(-1) => A), (AT(0) => A), (AT(1) => B)]))
length(collars)
collars

function SubstitutionTilings.draw(ptile::ATTile, action)
    fτ = 1
    if ptile == A
        return Luxor.poly([
            Point(-fτ*0.5,-100),
            Point(-fτ*0.5,100),
            Point(fτ*0.5,100),
            Point(fτ*0.5,-100)
        ], close = true, action)
    else
        return Luxor.poly([
            Point(-0.5,-100),
            Point(-0.5,100),
            Point(0.5,100),
            Point(0.5,-100)
        ], close = true, action)
    end
end

function fcolor(tile :: Pair{AT,ATTile})
    return (tile[2] == A) ? 1 : 2
end

n= 4
tiling = substitute(at, [AT(0) => A], n)
tiling

sc = 1
width = div(5^n*sc*2,2)
height = 30

@pdf begin
    colors = ["#5FD9A1", "#99706A"]
    setline(1)

    for tile in tiling
        origin()
        draw(AT(-div(5^(n),2))*tile, sc, colors[fcolor(tile)], :fillstroke)
        #sethue("black")
        #draw(tile[2], :stroke)
    end

    origin()
    scale(sc)
    sethue("black")
    """Luxor.poly([
            Point(-0.05,-10),
            Point(-0.05,10),
            Point(0.05,10),
            Point(0.05,-10)
        ], close = true, :fill)"""
end width height "anti-1"