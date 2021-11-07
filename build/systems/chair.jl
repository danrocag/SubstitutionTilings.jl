module Chair
__precompile__(false)

export chair, chair_system, ChairElem, ChairPTile

using Nemo
using Luxor
import Luxor.transform
using Colors
import AbstractAlgebra
using StructEquality

using ...CoreDefs


@def_structequal struct ChairElem
    angle :: Int # should be 0 ≤ i ≤ 3
    x :: Int
    y :: Int
end

struct ChairPTile end

function chair(angle, x, y)
    return (ChairElem(angle, x, y), ChairPTile())
end

function sin(i :: Integer)
    if i == 0 || i == 2
        return 0
    elseif i == 1
        return 1
    else
        return -1
    end
end


function cos(i :: Integer)
    return sin((i+1) % 4 )
end


function CoreDefs.:*(s :: ChairElem, t :: ChairElem)
    return ChairElem(
        (s.angle + t.angle) % 4,
        s.x + cos(s.angle)*t.x + sin(s.angle)*t.y,
        s.y - sin(s.angle)*t.x + cos(s.angle)*t.y
    )
end


function CoreDefs.dilate(λ :: Int, t :: ChairElem)
    return ChairElem(t.angle, λ*t.x, λ*t.y)
end


chair_subst = Dict([
    (ChairPTile(), [
        (ChairElem(0, -1, -1), ChairPTile()),
        (ChairElem(0, 1, 1), ChairPTile()),    
        (ChairElem(1, -1, 5), ChairPTile()),
        (ChairElem(3, 5, -1), ChairPTile()),
    ])
])

chair_system = SubSystem(chair_subst, 2)



function CoreDefs.embed_aff(t :: ChairElem)
    return float([
        cos(t.angle), -sin(t.angle),
        sin(t.angle), cos(t.angle),
        t.x, t.y])
end


function CoreDefs.draw(::ChairPTile, action)
    Luxor.poly([
        Point(-1, -1),
        Point(-1, 3),
        Point(1, 3),
        Point(1, 1),
        Point(3, 1),
        Point(3, -1),
    ], close = true, action)
end

origin = chair(0,0,0)

function in_bounds(t :: Tuple{ChairElem, ChairPTile}, n, window)
    return 2^(n+1)*(abs(t[1].x)) < window.w + 2^(n+1)*3 && 2^(n+1)*(abs(t[1].y)) < window.h + 2^(n+1)*3
end

end