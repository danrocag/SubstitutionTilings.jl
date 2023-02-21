using SubstitutionTilings
using SubstitutionTilings.Pinwheel

using Luxor

using Bessels
using LinearAlgebra
using Plots

Qθ = Pinwheel.Qθ

selfsimilar = wheel(0,0,0,(-1-i)//4)
selfsimilar in substitute(pinwheel(), [selfsimilar], 2)
w = 800
h = 800
sc = 30
@pdf begin
    colors = ["#DD93FC", "#E7977A",]

    tiling = substitute(pinwheel(), [selfsimilar], 6)
    setline(1)

    for tile in tiling
        origin()
        draw(tile, sc, colors[Pinwheel.color(tile)], :fill)
        origin()
        #draw(tile, sc, "black", :stroke)
    end
    #setline(2)
    #draw(selfsimilar, sc, "black", :stroke)
end w h "pinwheel-tiling"

width = 800
height = 800
sc = 100
@pdf begin
    colors = ["#DD93FC", "#E7977A",]
    quadrants = Tiler(width, height, 2, 1, margin=5)

    for (pos, n) in quadrants
        first_tile = Pinwheel.PinwheelElem(0,0,0,-1//2*i-Qθ(1//2))*wheel()
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
            sethue("black")
            transform(embed_aff(tile[1]))
            setline(0.0000001)
            #draw(tile[2], :stroke)
        end
    end
end width height "pinwheel-rule"

vertex_star = Dict([wheel(), wheel(0,0,1,Qθ(0)), wheel(0,1,0,-i)])
substitute(pinwheel(), [wheel()], 3)
initial_collar = inv(selfsimilar[1])*collar_in(Dict(substitute(pinwheel(), [selfsimilar], 2)), selfsimilar[1])

vertex_star_2 = Dict([wheel(), wheel(0,0,1,Qθ(0)), wheel(0,2,0,Qθ(0)), wheel(0,2,1,Qθ(0))])
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


function balanced(_ :: Pinwheel.PinwheelPTile, t :: Pair{Pinwheel.PinwheelElem, Pinwheel.PinwheelPTile})
    t[1].refl ? 1 : -1
end
# 3rd iteration in 45 seconds with nf_field
@time nu = autocorrelation(pinwheel(), initial_collar, 3, weights=balanced)
maximum(norm.(Pinwheel.embed_float.(keys(nu))))

xs = 0.0001:0.001:5
ys = zeros(length(xs))

R = 20
count = 0
for (i,j) in nu
    r = norm(Pinwheel.embed_float(i))
    if true
        count += 1
        ys += j*besselj0.(2*pi*r*xs)
    end
end
plot(xs,ys,xticks=0:0.5:10)

Is = zeros(length(xs))
for i in 1:length(xs) 
    Is[i] = sum(ys[1:i].*xs[1:i])
end
hs = log.(abs.(Is))./log.(xs)
plot(xs[1:50],hs[1:50],xticks=0:0.5:10)