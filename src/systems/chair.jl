module Chair

export chair, chair_system, ChairElem, ChairPTile

using Luxor
import Luxor.transform
using StructEquality

using ...CoreDefs


@struct_hash_equal_isequal struct ChairElem <: DGroupElem
    angle :: Int # should be 0 ≤ i ≤ 3
    x :: Int
    y :: Int
end

struct ChairPTile end

function chair(angle, x, y)
    return ChairElem(angle, x, y) => ChairPTile()
end

function sin(i :: Integer)
    ix = mod(i,4)
    if ix == 0 || ix == 2
        return 0
    elseif ix == 1
        return 1
    else
        return -1
    end
end


function cos(i :: Integer)
    return sin(i+1)
end


function Base.:*(s :: ChairElem, t :: ChairElem)
    return ChairElem(
        mod(s.angle + t.angle,4),
        s.x + cos(s.angle)*t.x + sin(s.angle)*t.y,
        s.y - sin(s.angle)*t.x + cos(s.angle)*t.y
    )
end
function Base.:*(s :: ChairElem, t :: Tuple{Integer, Integer})
    return [
        s.x + cos(s.angle)*t[1] + sin(s.angle)*t[2],
        s.y - sin(s.angle)*t[1] + cos(s.angle)*t[2],
    ]
end


function CoreDefs.dilate(λ :: Int, t :: ChairElem)
    return ChairElem(t.angle, λ*t.x, λ*t.y)
end
function CoreDefs.id(::Type{ChairElem})
    return ChairElem(0,0,0)
end
function Base.inv(g::ChairElem)
    return ChairElem(-g.angle,0,0)*ChairElem(0,-g.x, -g.y)
end

chair_subst = Dict([
    (ChairPTile(), [
        ChairElem(0, -1, -1) => ChairPTile(),
        ChairElem(0, 1, 1) => ChairPTile(),    
        ChairElem(1, -1, 5) => ChairPTile(),
        ChairElem(3, 5, -1) => ChairPTile(),
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

const origin = chair(0,0,0)

function in_bounds(t :: Pair{ChairElem, ChairPTile}, n, window)
    return 2^(n+1)*(abs(t[1].x)) < window.w + 2^(n+1)*3 && 2^(n+1)*(abs(t[1].y)) < window.h + 2^(n+1)*3
end


function squares(g)
    return g*[
        (0,0),
        (0,2),
        (2,0),
    ]
end
function adjacent(xs,ys)
    for x in xs
        for y in ys
            if max(abs.(x-y)...) ≤ 2
                return true
            end
        end
    end
    return false
end
function CoreDefs.collar_in(tiling :: Dict{ChairElem, ChairPTile}, g::ChairElem)
    result = Pair{ChairElem, ChairPTile}[]
    squares_g = tiling[g]
    for (h,l) in tiling
        if adjacent(squares(g), squares(h))
            push!(result, h => l)
        end
    end
    return Dict(result)
end

end