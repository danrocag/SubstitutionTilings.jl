using SubstitutionTilings
using SubstitutionTilings.Chair
using Test

using Luxor
import Luxor: translate

using LinearAlgebra
using Plots

# Code to generate chair tiling pictures for presentations and such
width = 1200
height = 800
sc = 10
@pdf begin
    colors = ["#3CD0E6", "#CA7EE6", "#E6873C", "#B4E647"]
    init = Chair.ChairElem(0,1,1)substitute(chair_system, [chair(0,0,0)], 2)
    tiling = substitute(chair_system, init, 6, Chair.in_bounds, (w=width/sc, h=height/sc))

    for tile in tiling
        
        setopacity(1)
        draw(tile, sc, colors[tile[1].angle+1], :fill)
        
        setline(0.1)
        setopacity(0.3)
        draw(tile, sc, "black", :stroke)
        
        scale(sc)
        translate(tile[1].x, tile[1].y)
        sethue("black")
        setopacity(1)
        if true
            circle(O, .2, :fillstroke)
        end
        # fontsize(1)
        # Luxor.text(string(tile[1].angle+1), halign=:center, valign=:middle)
    end
end width height "chair-tiling-beamer-total"

@pdf begin
    colors = ["#3CD0E6", "#CA7EE6", "#E6873C", "#B4E647"]
    sc = 10
    init = Chair.ChairElem(0,1,1)substitute(chair_system, [chair(0,0,0)], 2)
    tiling = substitute(chair_system, init, 4, Chair.in_bounds, (w=80, h=80))

    for tile in tiling
        
        setopacity(1)
        draw(tile, sc, colors[tile[1].angle+1], :fill)
        
        setline(0.1)
        setopacity(0.3)
        draw(tile, sc, "black", :stroke)
        
        scale(sc)
        translate(tile[1].x, tile[1].y)
        sethue("black")
        setopacity(1)
        circle(O, 0.1, :fill)
    end
end 800 800 "chair-tiling"

width = 900
height = 300
sc = 10
@pdf begin
    colors = ["#3CD0E6", "#CA7EE6", "#E6873C", "#B4E647"]
    
    t = Table([300], [50,150,250,100,300])

    for (pos, n) in t

        if t.currentcol % 2 == 0
            
            translate(pos)
            sethue("black")
            Luxor.arrow(Point(-0, 0), Point(75, 0))
            continue
        end

        first_tile = chair(0,0,0)
        tiling = substitute(chair_system, [first_tile], div(t.currentcol,2))
        for tile in tiling
            println(typeof(tile))
            
            translate(pos)
            scale(sc)
            sethue(colors[tile[1].angle+1])
            transform(embed_aff(tile[1]))
            draw(tile[2], :fill)
            
            translate(pos)
            setline(0.1)
            setopacity(0.3)
            draw(tile, sc, "black", :stroke)
            setopacity(1)
        end
    end
end width height "chair-rule-3"

width = 900
height = 300
sc = 5
@pdf begin
    colors = ["#3CD0E6", "#CA7EE6", "#E6873C", "#B4E647"]
    
    t = Table([300], [50,50,250,50,200, 50, 250])

    for (pos, n) in t

        if t.currentcol % 2 == 0
            
            translate(pos)
            sethue("black")
            Luxor.arrow(Point(-0, 0), Point(75, 0))
            continue
        end

        first_tile = chair(0,0,0)
        tiling = substitute(chair_system, [first_tile], div(t.currentcol,2))
        for tile in tiling
            println(typeof(tile))
            
            translate(pos)
            scale(sc)
            sethue(colors[tile[1].angle+1])
            transform(embed_aff(tile[1]))
            draw(tile[2], :fill)
            
            translate(pos)
            setline(0.1)
            setopacity(0.3)
            draw(tile, sc, "black", :stroke)
            setopacity(1)
            
            translate(pos)
            scale(sc)
            translate(tile[1].x, tile[1].y)
            sethue("black")
            circle(O, .5, :stroke)
            fontsize(1)
            Luxor.text(string(tile[1].angle+1), halign=:center, valign=:middle)
        end
    end
end width height "chair-rule-4-numbered"

width = 600
height = 300
sc = 20
@pdf begin
    colors = ["#3CD0E6", "#CA7EE6", "#E6873C", "#B4E647"]
    
    t = Table([300], [250,100,250])
    quadrants = Tiler(width, height, 2, 2, margin=5)

    for (pos, n) in t

        if t.currentcol == 2
            
            translate(pos)
            sethue("black")
            Luxor.arrow(Point(-50, 0), Point(50, 0))
            continue
        end

        first_tile = chair(0,0,0)
        tiling = substitute(chair_system, [first_tile], t.currentcol == 1 ? 0 : 1)
        for tile in tiling
            println(typeof(tile))
            
            translate(pos)
            scale(sc)
            sethue(colors[tile[1].angle+1])
            transform(embed_aff(tile[1]))
            draw(tile[2], :fill)
            
            translate(pos)
            setline(0.1)
            setopacity(0.3)
            draw(tile, sc, "black", :stroke)
            setopacity(1)
        end
    end
end width height "chair-rule"

width = 600
height = 300
sc = 20
@pdf begin
    colors = ["#3CD0E6", "#CA7EE6", "#E6873C", "#B4E647"]
    
    t = Table([300], [250,100,250])
    quadrants = Tiler(width, height, 2, 2, margin=5)

    for (pos, n) in t

        if t.currentcol == 2
            
            translate(pos)
            sethue("black")
            Luxor.arrow(Point(-50, 0), Point(50, 0))
            continue
        end

        first_tile = chair(3,0,0)
        tiling = substitute(chair_system, [first_tile], t.currentcol == 1 ? 0 : 1)
        for tile in tiling
            println(typeof(tile))
            
            translate(pos)
            scale(sc)
            sethue(colors[tile[1].angle+1])
            transform(embed_aff(tile[1]))
            draw(tile[2], :fill)
            
            translate(pos)
            setline(0.1)
            setopacity(0.3)
            draw(tile, sc, "black", :stroke)
            setopacity(1)
        end
    end
end width height "chair-rule-rot4"

width = 600
height = 300
sc = 20
@pdf begin
    colors = ["#3CD0E6", "#CA7EE6", "#E6873C", "#B4E647"]
    
    t = Table([300], [250,100,250])
    quadrants = Tiler(width, height, 2, 2, margin=5)

    for (pos, n) in t

        if t.currentcol == 2
            
            translate(pos)
            sethue("black")
            Luxor.arrow(Point(-50, 0), Point(50, 0))
            continue
        end

        first_tile = chair(3,0,0)
        tiling = substitute(chair_system, [first_tile], t.currentcol == 1 ? 0 : 1)
        for tile in tiling
            println(typeof(tile))
            
            translate(pos)
            scale(sc)
            sethue(colors[tile[1].angle+1])
            transform(embed_aff(tile[1]))
            draw(tile[2], :fill)
            
            translate(pos)
            setline(0.1)
            setopacity(0.3)
            draw(tile, sc, "black", :stroke)
            setopacity(1)
            if true
                
                translate(pos)
                scale(sc)
                translate(tile[1].x, tile[1].y)
                sethue("black")
                circle(O, .5, :stroke)
                fontsize(1)
                Luxor.text(string(tile[1].angle+1), halign=:center, valign=:middle)
            end
        end
    end
end width height "chair-rule-numbered-rot4"

width = 350
height = 4*100
sc = 6
@pdf begin
    colors = ["#3CD0E6", "#CA7EE6", "#E6873C", "#B4E647"]
    
    t = Table([100, 100, 100, 100], [150,50,150])

    for (pos, n) in t

        if t.currentcol == 2
            
            translate(pos)
            sethue("black")
            Luxor.arrow(Point(-50, 0), Point(50, 0))
            continue
        end

        first_tile = chair(t.currentrow-1,0,0)
        tiling = substitute(chair_system, [first_tile], t.currentcol == 1 ? 0 : 1)
        for tile in tiling
            println(typeof(tile))
            
            translate(pos)
            scale(sc)
            sethue(colors[tile[1].angle+1])
            transform(embed_aff(tile[1]))
            draw(tile[2], :fill)
            if t.currentcol == 3
                sethue("black")
                Luxor.box(O,0.5,0.5,:fill)
            end
            
            translate(pos)
            setline(0.1)
            setopacity(0.3)
            draw(tile, sc, "black", :stroke)
            setopacity(1)
            
            translate(pos)
            circle(O,1, :fill)
        end
    end
end width height "chair-rule-total"

width = 520
height = 130
sc = 20
@pdf begin
    colors = ["#3CD0E6", "#CA7EE6", "#E6873C", "#B4E647"]
    quadrants = Tiler(width, height, 1, 4, margin=5)

    for (pos, n) in quadrants
        first_tile =  ChairElem(n-1,0,0)*chair(0,0,0)
        tiling = [first_tile]
        for tile in tiling
            
            translate(pos)
            scale(sc)
            sethue(colors[tile[1].angle+1])
            setopacity(1)
            transform(embed_aff(tile[1]))
            draw(tile[2], :fill)

            
            translate(pos)
            setline(0.1)
            setopacity(0.3)
            draw(tile, sc, "black", :stroke)


            
            translate(pos)
            scale(sc)
            sethue("black")
            setopacity(1)
            transform(embed_aff(tile[1]))
            
            translate(pos)
            scale(sc)
            circle(O, .5, :stroke)
            fontsize(1)
            Luxor.text(string(tile[1].angle+1), halign=:center, valign=:middle)
        end
    end
end width height "chair-prototiles-numbered"

width = 400
height = 400
sc = 15
@pdf begin
    colors = ["#3CD0E6", "#CA7EE6", "#E6873C", "#B4E647"]
    init = substitute(chair_system, [chair(0,0,0)], 4)
    tiling = ChairElem(4,0,1)*substitute(chair_system, init, 6, Chair.in_bounds, (w=width/sc*2, h=height/sc*2))

    w = [1, 3]
    quadrants = Tiler(width, height, 1, 4, margin=5)
    for tile in tiling
        
        setopacity(1)
        draw(tile, sc, colors[tile[1].angle+1], :fill)
        
        setline(0.1)
        setopacity(0.3)
        draw(tile, sc, "black", :stroke)
        
        scale(sc)
        translate(tile[1].x, tile[1].y)
        sethue("black")
        setopacity(1)
        if tile[1].angle in w
            circle(O, 0.15, :fill)
        end
    end
end width height "chair-tiling-weight2"

@draw begin
    colors = ["#3CD0E6", "#CA7EE6", "#E6873C", "#B4E647"]
    sc = 10
    init = Chair.ChairElem(0,1,1)substitute(chair_system, [chair(0,0,0)], 2)
    #tiling = substitute(chair_system, init, 0, Chair.in_bounds, (w=80, h=80))
    tiling = initial_collar

    for tile in tiling
        
        setopacity(1)
        draw(tile, sc, colors[tile[1].angle+1], :fill)
        
        setline(0.1)
        setopacity(0.3)
        draw(tile, sc, "black", :stroke)
    end
end 800 800

init[2][1]
initial_collar = collar_in(Dict(init),init[2][1])

(collars, Sc) = total_collaring(chair_system, initial_collar)

@draw begin
    colors = ["#3CD0E6", "#CA7EE6", "#E6873C", "#B4E647"]
    sc = 10
    init = Chair.ChairElem(0,1,1)substitute(chair_system, [chair(0,0,0)], 2)
    #tiling = substitute(chair_system, init, 0, Chair.in_bounds, (w=80, h=80))
    tiling = collars[3]

    for tile in tiling
        
        setopacity(1)
        draw(tile, sc, colors[tile[1].angle+1], :fill)
        
        setline(0.1)
        setopacity(0.3)
        draw(tile, sc, "black", :stroke)
    end
end 800 800


N = 3

@time supertile = substitute(chair_system, [chair(0,0,0)], N)
sizeof(supertile)/2.0^30

xs = [t[1].x/2.0^N for t in supertile]
ys = [t[1].y/2.0^N for t in supertile]
#weights = [real((-1)^t[1].angle) for t in supertile]
K = 6
weights = [imag((im)^t[1].angle) for t in supertile].*2^K./2^N

histogram2d(xs, ys, bins=2^K, aspect_ratio=1, weights=weights,
    show_empty_bins=true,
    size=(2^10,2^10),
    color=:bam10)


N = 12
sums = zeros(ComplexF64, N)
for n=1:N
    supertile = substitute(chair_system, [chair(0,0,0)], n)
    xs = [t[1].x/2.0^n for t in supertile]
    ys = [t[1].y/2.0^n for t in supertile]
    weights = [im^(t[1].angle) for t in supertile]
    sums[n] = sum(exp.(-2*π*im*xs*1e-5).*weights)
end
sums
plot(abs.(sums)./(2.0.^(1:N)))


nu = @time autocorrelation(chair_system, initial_collar, 6, 1)
k = collect(keys(nu))
R = ChairElem(1,0,0)
k[200]
R*k[200]*inv(R)
nu[k[199]]
nu[R*R*k[199]*inv(R*R)]
frequency(chair_system, collar, cross, 4)
@time sizeof(substitute(chair_system, [chair(0,0,0)], 20, Chair.in_bounds, (w=800, h=800)))
 