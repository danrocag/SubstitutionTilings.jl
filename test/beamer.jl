using SubstitutionTilings
using SubstitutionTilings.Penrose
import SubstitutionTilings.Penrose: ψ
using Test

using Luxor

L = Penrose.Qζ


width = 1000*4
height = 1000*3
sc = 80
@draw begin
    colors = ["#DD93FC", "#E7977A", "#9B70AF", "#A0644F",]
    first_tile = hkite(0, false, L(0))
    tiling = substitute(penrose(), [first_tile], 12, Penrose.in_bounds, (w=width/sc, h=height/sc))
    #tiling  = pentagon
    println(typeof(tiling))
    setline(0.1)

    for tile in tiling
        origin()
        draw(tile, sc, colors[Penrose.color(tile)], :fill)
        draw(tile, sc, "black", :stroke)
    end
    origin()
    setline(10)
    draw(first_tile, sc, "black", :stroke)
end width height
@draw begin
    colors = ["#DD93FC", "#E7977A", "#9B70AF", "#A0644F",]
    first_tile = hkite(0, false, L(0))
    tiling = substitute(penrose(), [first_tile], 12, Penrose.in_bounds, (w=width/sc, h=height/sc))
    #tiling  = pentagon
    println(typeof(tiling))
    setline(0.1)

    for tile in tiling
        origin()
        draw(tile, sc, colors[Penrose.color(tile)], :fill)
        draw(tile, sc, "black", :stroke)
    end
    setline(10)
    origin()
    draw(first_tile, sc, "black", :stroke)
end width height# "canonical1.png"
@png begin
    colors = ["#DD93FC", "#E7977A", "#9B70AF", "#A0644F",]
    first_tile = hkite(0, false, L(0))
    tiling = substitute(penrose(), [first_tile], 12, Penrose.in_bounds, (w=width/sc, h=height/sc))
    #tiling  = pentagon
    println(typeof(tiling))
    setline(0.1)

    for tile in tiling
        origin()
        draw(tile, sc, colors[Penrose.color(tile)], :fill)
        draw(tile, sc, "black", :stroke)
    end
    first_tile = hkite(0, false, L(0))
    tiling = substitute(penrose(), [first_tile], 11, Penrose.in_bounds, (w=width/sc, h=height/sc))
    #tiling  = pentagon
    println(typeof(tiling))
    setline(2)

    for tile in tiling
        origin()
        draw(tile, sc*1.618, "black", :stroke)
    end
    setline(10)
    draw(first_tile, sc, "black", :stroke)
end width height "canonical2.png"
@png begin
    colors = ["#DD93FC", "#E7977A", "#9B70AF", "#A0644F",]
    first_tile = hkite(0, false, L(0))
    tiling = substitute(penrose(), [first_tile], 12, Penrose.in_bounds, (w=width/sc, h=height/sc))
    #tiling  = pentagon
    println(typeof(tiling))
    setline(0.1)

    for tile in tiling
        origin()
        draw(tile, sc, colors[Penrose.color(tile)], :fill)
        draw(tile, sc, "black", :stroke)
    end
    first_tile = hkite(0, false, L(0))
    tiling = substitute(penrose(), [first_tile], 10, Penrose.in_bounds, (w=width/sc, h=height/sc))
    #tiling  = pentagon
    println(typeof(tiling))
    setline(2)

    for tile in tiling
        origin()
        draw(tile, sc*1.618^2, "black", :stroke)
    end
    setline(10)
    draw(first_tile, sc, "black", :stroke)
end width height "canonical3.png"


width = 800
height = 600
@png begin
    sc = 80
    colors = ["#DD93FC", "#E7977A", "#9B70AF", "#A0644F",]
    quadrants = Tiler(width, height, 2, 2, margin=5)

    for (pos, n) in quadrants
        if n % 2 == 0
            first_tile = hkite(0, false, L(0))
        else
            first_tile = hdart(0, false, L(0))
        end
        tiling = substitute(penrose(), [first_tile], div(n-1,2),)
        for tile in tiling
            origin()
            translate(pos)
            scale(sc)
            sethue(colors[Penrose.color(tile)])
            transform(embed_aff(tile[1]))
            draw(tile[2], :fill)
        end
    end
end width height
@png begin
    sc = 120
    colors = ["#DD93FC", "#E7977A", "#9B70AF", "#A0644F",]
    first_tile = hkite(0, false, L(0))
    quadrants = Tiler(width, height, 2, 3, margin=5)

    for (pos, n) in quadrants
        tiling = substitute(penrose(), [first_tile], n,)
        for tile in tiling
            origin()
            translate(pos)
            scale(sc/1.618^n)
            sethue(colors[Penrose.color(tile)])
            transform(embed_aff(tile[1]))
            draw(tile[2], :fill)
        end
    end
end width height
@png begin
    sc = 120
    colors = ["#DD93FC", "#E7977A", "#9B70AF", "#A0644F",]
    first_tile = hdart(0, false, L(0))
    quadrants = Tiler(width, height, 2, 3, margin=5)

    for (pos, n) in quadrants
        tiling = substitute(penrose(), [first_tile], n,)
        for tile in tiling
            origin()
            translate(pos)
            scale(sc/1.618^n)
            sethue(colors[Penrose.color(tile)])
            transform(embed_aff(tile[1]))
            draw(tile[2], :fill)
        end
    end
end width height
@png begin
    sc = 80
    setline(1)
    colors = ["#DD93FC", "#E7977A", "#9B70AF", "#A0644F",]
    quadrants = Tiler(width, height, 1, 2, margin=5)

    for (pos, n) in quadrants
        if n % 2 == 0
            first_tile = hkite(0, false, L(0))
        else
            first_tile = hdart(0, false, L(0))
        end
        tiling = Penrose.force([first_tile])
        for tile in tiling
            origin()
            translate(pos)
            scale(sc)
            sethue(colors[Penrose.color(tile)])
            transform(embed_aff(tile[1]))
            draw(tile[2], :fill)
        end
        origin()
        translate(pos)
        scale(sc)
        sethue("black")
        transform(embed_aff(first_tile[1]))
        draw(first_tile[2], :stroke)
    end
end width height
@png begin
    sc = 18
    setline(1)
    colors = ["#DD93FC", "#E7977A", "#9B70AF", "#A0644F",]
    quadrants = Tiler(width, height, 1, 2, margin=80)

    for (pos, n) in quadrants
        if n % 2 == 0
            first_tile = hkite(0, false, L(0))
        else
            first_tile = hdart(0, false, L(0))
        end
        tiling = unique(substitute(forced_penrose(), ([first_tile]), 4))
        for tile in tiling
            origin()
            translate(pos)
            scale(sc)
            sethue(colors[Penrose.color(tile)])
            transform(embed_aff(tile[1]))
            draw(tile[2], :fill)
        end
        origin()
        translate(pos)
        scale(sc*1.618^4)
        sethue("black")
        transform(embed_aff(first_tile[1]))
        draw(first_tile[2], :stroke)
    end
end width height

width = 1000
height = 1000
sc = 80
@draw begin
    colors = ["#DD93FC", "#E7977A", "#9B70AF", "#A0644F",]
    first_tile = hkite(0, false, L(0))
    tiling = substitute(penrose(), [first_tile], 12, Penrose.in_bounds, (w=width/sc, h=height/sc))
    #tiling  = pentagon
    println(typeof(tiling))
    setline(1)

    for tile in tiling
        origin()
        draw(tile, sc, colors[Penrose.color(tile)], :fill)
        draw(tile, sc, "black", :stroke)
    end
    origin()
    setline(5)
    draw(first_tile, sc, "black", :stroke)
end width height
@draw begin
    colors = ["#DD93FC", "#E7977A", "#9B70AF", "#A0644F",]
    first_tile = hkite(0, false, L(0))
    tiling = substitute(penrose(), [first_tile], 11, Penrose.in_bounds, (w=width/sc, h=height/sc))
    #tiling  = pentagon
    println(typeof(tiling))
    setline(1)

    for tile in substitute(penrose(), tiling, 1)
        origin()
        draw(tile, sc, colors[Penrose.color(tile)], :fill)
        draw(tile, sc, "black", :stroke)
    end
    origin()
    setline(5)
    draw(collect(tiling)[80], sc*1.618, "black", :stroke)
    origin()
    setline(1)
    draw(first_tile, sc, "black", :stroke)
end width height
@draw begin
    colors = ["#DD93FC", "#E7977A", "#9B70AF", "#A0644F",]
    first_tile = hkite(0, false, L(0))
    preimg = substitute(penrose(), [first_tile], 7)
    tiling = substitute(penrose(), preimg, 1)
    g = inv(collect(tiling)[49][1])
    tiling = g*tiling
    preimg = dilate(ϕ, g)*preimg
    preimg_tile = collect(preimg)[600]
    setline(1)

    for tile in tiling
        origin()
        draw(tile, sc, colors[Penrose.color(tile)], :fill)
        draw(tile, sc, "black", :stroke)
    end
    origin()
    setline(1)
    draw(first_tile, sc, "black", :stroke)
    origin()
    setline(5)
    draw(preimg_tile, sc*1.618, "black", :stroke)
end width height