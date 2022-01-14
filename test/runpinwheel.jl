using SubstitutionTilings
using SubstitutionTilings.Pinwheel
using SubstitutionTilings.Collaring
using Test

using Luxor
using Hecke

using Profile

K = Pinwheel.K


w = 800
h = 800
sc = 60
@draw begin
    colors = ["#DD93FC", "#E7977A",]

    tiling = [wheel(), wheel(0,0,1,K(0)), wheel(0,2,0,K(0)), wheel(0,2,1,K(0))]
    setline(1)

    for tile in tiling
        origin()
        draw(tile, sc, colors[Pinwheel.color(tile)], :fill)
        draw(tile, sc, "black", :stroke)
    end
    setline(2)
    draw(first_tile, sc, "black", :stroke)
end w h

vertex_star = Dict([wheel(), wheel(0,0,1,K(0)), wheel(0,1,0,-i)])
initial_collar = collar_in(Dict(Pinwheel.PinwheelElem(0,0,0,-i-1)*substitute(pinwheel(), [first_tile], 2)), Pinwheel.PinwheelElem(0,0,0,K(0)))

vertex_star_2 = Dict([wheel(), wheel(0,0,1,K(0)), wheel(0,2,0,K(0)), wheel(0,2,1,K(0))])
@time .frequency(pinwheel(), initial_collar, vertex_star,2)
@time .frequency(pinwheel(), initial_collar, vertex_star_2,2) # 4x bigger than in Baake-Grimm because every tile gets counted

@draw begin
    colors = ["#DD93FC", "#E7977A",]

    tiling = collars[24]
    setline(1)

    for tile in tiling
        origin()
        draw(tile, sc, colors[Pinwheel.color(tile)], :fill)
        draw(tile, sc, "black", :stroke)
    end
    setline(2)
    draw(first_tile, sc, "black", :stroke)
end w h