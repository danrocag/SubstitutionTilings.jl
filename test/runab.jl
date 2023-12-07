using SubstitutionTilings
using SubstitutionTilings.AmmannBeenker
using Luxor

L = AmmannBeenker.L

L(1)
width = 400
height = 400
sc = 40
@pdf begin
    colors =  ["#AFB1E4", "#E0CA6C", "#948238"]
    quadrants = Tiler(width, height, 2, 2, margin=5)

    for (pos, n) in quadrants
        println(n)
        first_tile = div(n-1,2) == 0 ? rhomb() : square()
        tiling = substitute(ammannbeenker(), [ABElem(0,0,-L(div(n-1,2)*sq2//2))*first_tile], 1-n%2)
        for tile in tiling
            origin()
            sethue(colors[color(tile)])
            translate(pos)
            scale(sc)
            transform(embed_aff(ABElem(0,0,L(0))*tile[1]))
            draw(tile[2], :fill)
        end
        
    end
end width height "ab-rule.pdf"

width = 800
height = 800
sc = 40
@pdf begin
    background("black")
    colors = ["#AFB1E4", "#E0CA6C", "#948238"]
    first_tile = rhomb()
    tiling = substitute(ammannbeenker(), [ABElem(0,0,ζ^5//2)*first_tile], 6)


    for tile in tiling
        origin()
        draw(tile, sc, colors[color(tile)], :fill)
        origin()
    end
end width height "ab-tiling.pdf"