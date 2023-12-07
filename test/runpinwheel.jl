using SubstitutionTilings
using SubstitutionTilings.Pinwheel

using Luxor


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
        setline(1)
        setopacity(0.5)
        draw(tile, sc, "black", :stroke)
    end
end w h "pinwheel-tiling"

width = 600
height = 200
sc = 40
@pdf begin
    colors = ["#DD93FC", "#E7977A",]
    quadrants = Tiler(width, height, 1, 2, margin=5)

    for (pos, n) in quadrants
        first_tile = Pinwheel.PinwheelElem(0,0,0,Qθ(-1)-i//2)*wheel()
        tiling = substitute(pinwheel(), [first_tile], n-1)
        for tile in tiling
            origin()
            translate(pos)
            scale(sc)
            sethue(colors[Pinwheel.color(tile)])
            transform(embed_aff(tile[1]))
            draw(tile[2], :fill)
            origin()
            sethue("black")
            translate(pos)
            scale(sc)
            transform(embed_aff(tile[1]))
            setline(1)
            setopacity(0.4)
            draw(tile[2], :stroke)
            setopacity(1)
        end
    end
end width height "pinwheel-rule"

initial_collar = inv(selfsimilar[1])*collar_in(Dict(substitute(pinwheel(), [selfsimilar], 2)), selfsimilar[1])

vertex_star = Dict([wheel(), wheel(0,0,1,Qθ(0)), wheel(0,1,0,-i)])
vertex_star_2 = Dict([wheel(), wheel(0,0,1,Qθ(0)), wheel(0,2,0,Qθ(0)), wheel(0,2,1,Qθ(0))])
@time frequency(pinwheel(), initial_collar, vertex_star, 3)
@time frequency(pinwheel(), initial_collar, vertex_star_2, 2) # 4x bigger than in Baake-Grimm because every tile gets counted
@time Pinwheel.CoreDefs.frequencies(pinwheel(), initial_collar, [vertex_star, vertex_star, vertex_star, vertex_star, vertex_star], 2)


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

collars, Sc = total_collaring(pinwheel(), initial_collar)

@draw begin
    colors = ["#DD93FC", "#E7977A",]

    tiling = collars[2]
    setline(1)

    for tile in tiling
        origin()
        draw(tile, sc, colors[Pinwheel.color(tile)], :fill)
        origin()
        draw(tile, sc, "black", :stroke)
    end
    setline(2)
end w h



function balanced(_ :: Pinwheel.PinwheelPTile, t :: Pair{Pinwheel.PinwheelElem, Pinwheel.PinwheelPTile})
    t[1].refl ? 1 : -1
end
using Serialization
@time nu = autocorrelation(pinwheel(), initial_collar, 4);
# serialize("pinwheel_nu_5", nu)


nu = deserialize("pinwheel_nu_5")
nu_arr = collect(nu)
using LinearAlgebra
maximum(norm.(Pinwheel.embed_float.(keys(nu))))

function window_bump(r :: Float64)
    if 1 <= r
        return 0
    else
        return exp(-1/(1-r^2))
    end
    # integral = 0.46651254959873767
end
I_bump = 0.46651239317832216

Rs = 1:0.1:20
vars = zeros(length(Rs))
@time for (g, freq) in nu_arr
    r = norm(Pinwheel.embed_float(g))
    bump_weight = window_bump.(r./Rs)
    vars += freq*bump_weight
end
vars -= I_bump.*(Rs.^2)

using Plots

plot(Rs,log.(abs.(vars)))
plot!(Rs,log.(abs.(vars)))
vars_5 = vars