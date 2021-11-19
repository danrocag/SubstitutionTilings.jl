module HeisenbergImproper
__precompile__(false)

export He, heis, hn, left, right

using Nemo
using Luxor
import Luxor.transform
using Colors
using StructEquality

QQ = Rational{Int}

using ...CoreDefs


@def_structequal struct He <: EGroupElem
    a::QQ
    b::QQ
    c::QQ
end
@enum HePTile begin
    Left=1
    Right=2
end

function Base.:*(x :: He, y :: He) where {T}
    return He( x.a + y.a, x.b + y.b + x.a*y.c, x.c + y.c)
end

function CoreDefs.dilate(λ, x :: He)
    return He(x.a, λ * x.b, λ * x.c )
end

function Base.inv(x :: He)
    return He(-x.a, -x.b, x.a*x.b - x.c)
end

function left(a,b,c)
    return (He(a,b,c), Left)
end
function left()
    return left(0,0,0)
end

function right(a,b,c)
    return (He(a,b,c), Right)
end
function right()
    return right(0,0,0)
end

function heis()
    heisenberg_subst = Dict([
        (Left, [
            left(0,-2,-1),
            right(0,2,-1),
            left(-1,1,-1),
            left(0,2,3),
        ]),
        (Right, [
            right(0,-2,3),
            left(0,2,3),
            right(1,-1,3),
            right(0,2,-1),
        ])
    ])

    return SubSystem(heisenberg_subst, 2//1) 
end


function CoreDefs.embed_aff(x :: He)
    return float([1, 0, x.a, 1, x.b, x.c])
end

function CoreDefs.draw(p::HePTile, action)
    if p == Left
        Luxor.poly([
            Point(-2,-1),
            Point(2,3),
            Point(2,-1)
        ], action, close=true)
    else
        Luxor.poly([
            Point(-2,3),
            Point(2,3),
            Point(2,-1)
        ], action, close=true)
    end

end


function hn(n)
    substitute(heis(), [left()], n)
end

function color(tile)
    return Int(tile[2])
end

end