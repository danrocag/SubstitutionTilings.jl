using SubstitutionTilings
using SubstitutionTilings.Pinwheel

using Luxor

using Profile

K = Pinwheel.K


w = 800
h = 800
sc = 30
@png begin
    colors = ["#DD93FC", "#E7977A",]

    first_tile = wheel()
    tiling = Pinwheel.PinwheelElem(0,0,0,-i*20-20)*substitute(pinwheel(), [first_tile], 5)
    setline(1)

    for tile in tiling
        origin()
        draw(tile, sc, colors[Pinwheel.color(tile)], :fill)
        origin()
        draw(tile, sc, "black", :stroke)
    end
    #setline(2)
    #draw(first_tile, sc, "black", :stroke)
end w h "pinwheel-tiling.png"

@png begin
    w = 800
    h = 800
    sc = 100
    colors = ["#DD93FC", "#E7977A",]
    quadrants = Tiler(width, height, 2, 1, margin=5)

    for (pos, n) in quadrants
        first_tile = Pinwheel.PinwheelElem(0,0,0,-1//2*i-1//2)*wheel()
        tiling = substitute(pinwheel(), [first_tile], n-1)
        for tile in tiling
            origin()
            translate(pos)
            scale(sc)
            sethue(colors[Pinwheel.color(tile)])
            transform(embed_aff(tile[1]))
            draw(tile[2], :fill)
            origin()
            translate(pos)
            scale(sc)
            setline(0.5)
            sethue("black")
            transform(embed_aff(tile[1]))
            draw(tile[2], :stroke)
        end
    end
end width height "pinwheel-subst.png"

vertex_star = Dict([wheel(), wheel(0,0,1,K(0)), wheel(0,1,0,-i)])
initial_collar = collar_in(Dict(Pinwheel.PinwheelElem(0,0,0,-i-1)*substitute(pinwheel(), [first_tile], 2)), Pinwheel.PinwheelElem(0,0,0,K(0)))

vertex_star_2 = Dict([wheel(), wheel(0,0,1,K(0)), wheel(0,2,0,K(0)), wheel(0,2,1,K(0))])
@time frequency(pinwheel(), initial_collar, vertex_star,2)
@time frequency(pinwheel(), initial_collar, vertex_star_2,2) # 4x bigger than in Baake-Grimm because every tile gets counted

@png begin
    colors = ["#DD93FC", "#E7977A",]

    tiling = vertex_star
    setline(1)

    for tile in tiling
        origin()
        draw(tile, sc, colors[Pinwheel.color(tile)], :fill)
        origin()
        draw(tile, sc, "black", :stroke)
    end
    setline(2)
    #draw(first_tile, sc, "black", :stroke)
end w h "vertex_star_1.png"

@png begin
    colors = ["#DD93FC", "#E7977A",]

    tiling = vertex_star_2
    setline(1)

    for tile in tiling
        origin()
        draw(tile, sc, colors[Pinwheel.color(tile)], :fill)
        origin()
        draw(tile, sc, "black", :stroke)
    end
    setline(2)
    #draw(first_tile, sc, "black", :stroke)
end w h "vertex_star_2.png"