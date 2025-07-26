using SubstitutionTilings
using SubstitutionTilings.Penrose
import SubstitutionTilings.Penrose: ψ, ϕ, ζ
using Plots
using Luxor
import Luxor: translate
using Bessels

const L = Penrose.Qζ

Penrose.embed_nf(ζ)
#collars = SubstitutionTilings.CoreDefs.collars(penrose(), 6);
#length(collars)

@draw juliacircles()

width = 1200
height = 1200
sc = 20
@draw begin
    #setantialias(6)
    colors = ["#DD93FC", "#E7977A", "#9B70AF", "#A0644F",]
    first_tile = hkite()
    tiling = substitute(penrose(), [first_tile], 10, Penrose.in_bounds, (w=width/sc, h=height/sc))
    #tiling = initial_collar

    for tile in tiling
        origin()
        draw(tile, sc,  "grey", :stroke)
        #rect(Point(-0.5,-0.01), 0.4,0.02, action=:fill)
        #rect(Point(-0.5,-0.1), 0.02,0.1, action=:fill)
    end
    for tile in tiling
        origin()
        draw(tile, sc, colors[Penrose.color(tile)], :fill)
        origin()
        setline(0.5)
        setopacity(0.4)
        draw(tile, sc, "#DD93FC", :stroke)
        setopacity(1)
    end
    origin()
    setline(1)
    #draw(first_tile, sc/0.618, "black", :stroke)
end width height



width = 800
height = 800
sc = 60
@draw begin
    colors = ["#DD93FC", "#E7977A", "#9B70AF", "#A0644F",]
    #tiling = initial_collar
    println(typeof(tiling))
    quadrants = Tiler(width, height, 2, 2, margin=5)

    for (pos, n) in quadrants
        println(n)
        first_tile = n % 2 == 0 ? hkite() : hdart()
        tiling = substitute(penrose(), [first_tile], div(n-1,2))
        
        for tile in tiling
            origin()
            Luxor.translate(pos)
            scale(sc)
            sethue(colors[Penrose.color(tile)])
            transform(embed_aff(tile[1]))
            draw(tile[2], :fill)
            origin()
            sethue("black")
            Luxor.translate(pos)
            scale(sc)
            transform(embed_aff(tile[1]))
            setline(1)
            setopacity(0.5)
            draw(tile[2], :stroke)
            setopacity(1)
        end
        origin()
        circle(pos, 1, :fill)
    end
end width height #"penrose-rule"

width = 1280*2
height = 800*2
sc = 400
@draw begin
    colors = ["#DD93FC", "#E7977A", "#9B70AF", "#A0644F",]

    quadrants = Tiler(width, height, 1, 2, margin=5)

    begin
        pos = quadrants[1][1]
        first_tile = hkite()
        tiling = substitute(penrose(), [first_tile], 0, Penrose.in_bounds, (w=width/sc, h=height/sc))
        println(typeof(tiling))
        setline(1)
    
        for tile in tiling
            origin()
            Luxor.translate(pos)
            draw(tile, sc, colors[Penrose.color(tile)], :fill)
        end
        origin()
        Luxor.translate(pos)
        scale(sc)
        setline(2)
        setcolor("black")
        Luxor.line(Penrose.embed_nf_p(ζ^4),Penrose.embed_nf_p(L(1)), :stroke)
        Luxor.arrow(Penrose.embed_nf_p(ζ^4),Penrose.embed_nf_p(ζ^4+(L(1)-ζ^4)*(7//16)), arrowheadlength=0.1)
        Luxor.arrow(Penrose.embed_nf_p(ζ^4),Penrose.embed_nf_p(ζ^4+(L(1)-ζ^4)*(9//16)), arrowheadlength=0.1)
        Luxor.line(Penrose.embed_nf_p(ζ^6),Penrose.embed_nf_p(L(1)), :stroke)
        Luxor.arrow(Penrose.embed_nf_p(L(1)),Penrose.embed_nf_p(ζ^6+(L(1)-ζ^6)//2), arrowheadlength=0.1)
        Luxor.arrow(Penrose.embed_nf_p(L(1)),Penrose.embed_nf_p(ζ^6+(L(1)-ζ^6)*(3//8)), arrowheadlength=0.1)
        Luxor.arrow(Penrose.embed_nf_p(L(1)),Penrose.embed_nf_p(ζ^6+(L(1)-ζ^6)*(5//8)), arrowheadlength=0.1)
        Luxor.line(Penrose.embed_nf_p(ζ^4),Penrose.embed_nf_p(ζ^6), :stroke)
        Luxor.arrow(Penrose.embed_nf_p(ζ^6),Penrose.embed_nf_p(ζ^4+(ζ^6-ζ^4)*(1//2)), arrowheadlength=0.1)
    end
    
    begin
        pos = quadrants[2][1]
        first_tile = hdart()
        tiling = substitute(penrose(), [first_tile], 0, Penrose.in_bounds, (w=width/sc, h=height/sc))
        println(typeof(tiling))
        setline(1)
    
        for tile in tiling
            origin()
            Luxor.translate(pos)
            draw(tile, sc, colors[Penrose.color(tile)], :fill)
        end
        origin()
        Luxor.translate(pos)
        scale(sc)
        setline(2)
        setcolor("black")
        Luxor.line(Penrose.embed_nf_p(ζ^2-ζ^0//ϕ),Penrose.embed_nf_p(L(1-ζ^0//ϕ)), :stroke)
        Luxor.arrow(Penrose.embed_nf_p(ζ^2-ζ^0//ϕ),Penrose.embed_nf_p(ζ^2-ζ^0//ϕ+(L(1-ζ^0//ϕ)-(ζ^2-ζ^0//ϕ))//2), arrowheadlength=0.1)
        Luxor.line(Penrose.embed_nf_p(L(1-ζ^0//ϕ)),Penrose.embed_nf_p(ζ^8-ζ^0//ϕ), :stroke)
        Luxor.arrow(Penrose.embed_nf_p(L(1-ζ^0//ϕ)),Penrose.embed_nf_p(1-ζ^0//ϕ+((ζ^8-ζ^0//ϕ)-(1-ζ^0//ϕ))*(3//16)), arrowheadlength=0.1)
        Luxor.arrow(Penrose.embed_nf_p(L(1-ζ^0//ϕ)),Penrose.embed_nf_p(1-ζ^0//ϕ+((ζ^8-ζ^0//ϕ)-(1-ζ^0//ϕ))*(5//16)), arrowheadlength=0.1)
        Luxor.arrow(Penrose.embed_nf_p(L(1-ζ^0//ϕ)),Penrose.embed_nf_p(1-ζ^0//ϕ+((ζ^8-ζ^0//ϕ)-(1-ζ^0//ϕ))*(7//16)), arrowheadlength=0.1)
        Luxor.arrow(Penrose.embed_nf_p(L(1-ζ^0//ϕ)),Penrose.embed_nf_p(1-ζ^0//ϕ+((ζ^8-ζ^0//ϕ)-(1-ζ^0//ϕ))*(9//16)), arrowheadlength=0.1)
        Luxor.line(Penrose.embed_nf_p(ζ^8-ζ^0//ϕ),Penrose.embed_nf_p(ζ^2-ζ^0//ϕ), :stroke)
        Luxor.arrow(Penrose.embed_nf_p(ζ^8-ζ^0//ϕ),Penrose.embed_nf_p(ζ^8-ζ^0//ϕ-(ζ^8-ζ^2)*(7//16)), arrowheadlength=0.1)
        Luxor.arrow(Penrose.embed_nf_p(ζ^8-ζ^0//ϕ),Penrose.embed_nf_p(ζ^8-ζ^0//ϕ-(ζ^8-ζ^2)*(9//16)), arrowheadlength=0.1)
   end
end width height #"penrose.png"

pentagon = (([
    PenroseElem(mod(k+s,10),s,L(1))*PenroseElem(0,0,L(-1)) => Penrose.Hkite
    for k=0:2:10 for s=0:1
]));

width = 300
height = 300
sc = 50
@png begin
    colors = ["#DD93FC", "#E7977A", "#9B70AF", "#A0644F",]
    first_tile = hkite(0, false, L(0))
    tiling = pentagon

    for tile in tiling
        origin()
        draw(tile, sc, colors[Penrose.color(tile)], :fill)
    end
    origin()
    setline(1)
    draw(first_tile, sc, "black", :stroke)
end width height "pentagon.png"

trapezoid = Dict([
    hkite(0,0,L(0)),
    hdart(2,1, ζ^8//ϕ )])
width = 300
height = 300
sc = 50
@png begin
    colors = ["#DD93FC", "#E7977A", "#9B70AF", "#A0644F",]
    first_tile = hkite(0, false, L(0))
    tiling = trapezoid

    for tile in tiling
        origin()
        draw(tile, sc, colors[Penrose.color(tile)], :fill)
    end
    origin()
    setline(1)
    draw(first_tile, sc, "black", :stroke)
end width height "trapezoid.png"

# Comparison of Array and Dict operations
@assert substitute(penrose(), Dict([hkite(0,0,L(0))]), 20, Penrose.in_bounds, (w=50, h=50)) == Dict(substitute(penrose(), ([hkite(0,0,L(0))]), 20, Penrose.in_bounds, (w=50, h=50)))

t_dict = zeros(35)
n_dict = zeros(35)
for i in 1:35
    p = substitute(penrose(), Dict([hkite(0,0,L(0))]), 20, Penrose.in_bounds, (w=i*2, h=i*2));
    time = @timed empirical_frequency(pentagon, Dict(p))
    t_dict[i] = time.time
    n_dict[i] = length(p)
end # essentially linear in the amont of tiles: pretty good!

t_arr = zeros(35)
n_arr = zeros(35)
for i in 1:35
    p = substitute(penrose(), ([hkite(0,0,L(0))]), 20, Penrose.in_bounds, (w=i*2, h=i*2));
    time = @timed empirical_frequency(pentagon, (p))
    t_arr[i] = time.time
    n_arr[i] = length(p)
end
plot(n_dict,t_dict)
plot!(n_arr, t_arr)
plot((n_dict), (t_dict))
plot(log.(n_arr), log.(t_arr))
# Frequency counting with arrays is exponential in the amount of tiles
# with Dicts it is more like linear

tiling = substitute(penrose(), [hkite(0,0,L(0))], 20, Penrose.in_bounds, (w=50, h=50))
i = argmin(abs.(Penrose.embed_nf.([t[1].z for t in tiling])))
@time initial_collar = collar_in(Dict(tiling), tiling[i][1], check=true)
frequency(penrose(), initial_collar, Dict(pentagon), 4)

function balanced(_ :: Penrose.PenrosePTile, t :: Pair{Penrose.PenroseElem, Penrose.PenrosePTile})
    t[1].refl ? 1 : -1
end
# 3rd iteration in 45 seconds with nf_field
#@time nu = autocorrelation(penrose(), initial_collar, 9);
using Serialization
nu = deserialize("penrose_nu_9")
nu_arr = collect(nu)y
rs = zeros(length(nu_arr))
freqs = zeros(length(nu_arr))
for j in eachindex(nu_arr)
    g, freq = nu_arr[j]
    rs[j] = abs(Penrose.embed_float(g))
    freqs[j] = freq
end

using AccurateArithmetic
maximum(abs.(Penrose.embed_float.(keys(nu))))
xs = 0.2:0.005:1
ys = zeros(length(xs))
@time for i=eachindex(xs)
    ys[i] = sum_kbn(freqs.*besselj0.(2*pi*rs*xs[i]))
end
plot(xs,ys)
plot!(xs,xs*0)
zs = zeros(length(xs))
@time for i=eachindex(xs)
    zs[i] = sum_kbn(ys[1:i])
end

plot(log.(zs)./log.((xs)))
empirical_frequency(pentagon, tiling)
Penrose.frequency(pentagon, 4)
@testset "Frequencies" begin
    @test Penrose.frequency([hkite()], 4) == ψ
    @test Penrose.frequency([hdart()], 4) == 1-ψ
    @test Penrose.frequency(pentagon, 4) == 7*ψ - 4
    @test Penrose.frequency(pentagon, 5) == 7*ψ - 4
    @test Penrose.frequency(pentagon, 6) == 7*ψ - 4
end

## Calculating window

function τ(x :: Qζ)
    galois_gen(x)
end

function τ(g :: PenroseElem)
    PenroseElem(mod(g.rot*3, 10), g.refl, τ(g.z))
end

function τ(t :: Pair{PenroseElem, SubstitutionTilings.Penrose.PenrosePTile})
    τ(t[1]) => t[2]
end

first_tile ∈ substitute(penrose(), [first_tile], 8)

window_approximant = τ.(substitute(penrose(), [first_tile], 12))
center_slice = [Penrose.embed_float(t[1]) for t=window_approximant if (t[1].refl && t[2]==SubstitutionTilings.Penrose.Hdart)]

scatter(center_slice, xlim=(-3,3), ylim=(-3,3), aspect_ratio=:equal)


## Empirical diffraction
centers = [Penrose.embed_float(t[1]) for t=substitute(penrose(), [hkite()], 9)]

function erf(t)
    exp(-t^2)
end

using AccurateArithmetic
using LinearAlgebra
xs = 0:0.0001:1
@time ys = 2*π*sum(besselj0.(π*norm(t1-t2)*xs) for t1=centers for t2=centers)
cutoff = 15
plot(xs[cutoff:end], ys[cutoff:end])