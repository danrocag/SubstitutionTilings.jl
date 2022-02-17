using SubstitutionTilings
using SubstitutionTilings.Chair
using SubstitutionTilings.Collaring
using Test

using Luxor


cross = [
    chair(0,0,0),
    chair(2,-2,-2),
    chair(1,0,-2),
    chair(3,-2,0),
]

tiling = substitute(chair_system, [chair(0,0,0)], 8, Chair.in_bounds, (w=80, h=80))
tiling_dict = Dict(tiling)
@time empirical_frequency(cross, tiling_dict)
@draw begin
    colors = ["#3CD0E6", "#CA7EE6", "#E6873C", "#B4E647"]
    sc = 10
    tiling = substitute(chair_system, [chair(0,0,0)], 8, Chair.in_bounds, (w=80, h=80))
    collar = substitute(chair_system, collar_class(tiling, ChairElem(0,1,1)), 0)

    setline(0.1)

    for tile in collar
        draw(tile, sc, colors[1 + tile[1].angle % 4], :fill)
        draw(tile, sc, "black", :stroke)
    end
end 800 800



frequency(chair_system, collar, cross, 4)
@time sizeof(substitute(chair_system, [chair(0,0,0)], 20, Chair.in_bounds, (w=800, h=800)))
 