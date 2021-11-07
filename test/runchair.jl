
using SubstitutionTilings
using SubstitutionTilings.Chair
using Test

using Luxor


@draw begin
    colors = ["#3CD0E6", "#CA7EE6", "#E6873C", "#B4E647"]
    sc = 10
    tiling = substitute(chair_system, [chair(0,0,0)], 10, Chair.in_bounds, (w=80, h=80))
    setline(0.1)

    for tile in tiling
        draw(tile, sc, colors[1 + tile[1].angle % 4], :fill)
        draw(tile, sc, "black", :stroke)
    end
end 800 800

@time sizeof(substitute(chair_system, [chair(0,0,0)], 20, Chair.in_bounds, (w=800, h=800)))
 