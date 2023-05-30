using SubstitutionTilings
using SubstitutionTilings.Chair
using Test

using Luxor


@pdf begin
    colors = ["#3CD0E6", "#CA7EE6", "#E6873C", "#B4E647"]
    sc = 10
    init = Chair.ChairElem(0,1,1)substitute(chair_system, [chair(0,0,0)], 2)
    tiling = substitute(chair_system, init, 4, Chair.in_bounds, (w=80, h=80))

    for tile in tiling
        origin()
        setopacity(1)
        draw(tile, sc, colors[tile[1].angle+1], :fill)
        origin()
        setline(0.1)
        setopacity(0.3)
        draw(tile, sc, "black", :stroke)
    end
end 800 800 "chair-tiling"

width = 400
height = 500
sc = 20
@pdf begin
    colors = ["#3CD0E6", "#CA7EE6", "#E6873C", "#B4E647"]
    quadrants = Tiler(width, height, 2, 1, margin=5)

    for (pos, n) in quadrants
        first_tile = chair(0,0,0)
        tiling = substitute(chair_system, [first_tile], n-1)
        for tile in tiling
            println(typeof(tile))
            origin()
            translate(pos)
            scale(sc)
            sethue(colors[tile[1].angle+1])
            setopacity(1)
            transform(embed_aff(tile[1]))
            draw(tile[2], :fill)
            origin()
            translate(pos)
            setline(0.1)
            setopacity(0.3)
            draw(tile, sc, "black", :stroke)
        end
    end
end width height "chair-rule"

frequency(chair_system, collar, cross, 4)
@time sizeof(substitute(chair_system, [chair(0,0,0)], 20, Chair.in_bounds, (w=800, h=800)))
 