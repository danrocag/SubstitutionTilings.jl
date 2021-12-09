using SubstitutionTilings
using SubstitutionTilings.Penrose
import SubstitutionTilings.Penrose: ψ
using Test

using Luxor

const L = Penrose.Qζ




#collars = SubstitutionTilings.CoreDefs.collars(penrose(), 6);
#length(collars)

width = 2560*4
height = 1600*4
sc = 80
@draw begin
    colors = ["#DD93FC", "#E7977A", "#9B70AF", "#A0644F",]
    first_tile = hkite(0, false, L(0))
    tiling = substitute(penrose(), [first_tile], 12, Penrose.in_bounds, (w=width/sc, h=height/sc))
    println(typeof(tiling))
    setline(1)

    for tile in tiling
        origin()
        draw(tile, sc, colors[Penrose.color(tile)], :fill)
    end
    draw(first_tile, sc, "black", :stroke)
end width height

pentagon = [
    PenroseElem(k+s,s,L(1))*PenroseElem(0,0,L(-1)) => Penrose.Hkite
    for k=0:2:10 for s=0:1
];
p = @time substitute(penrose(), [hkite(0,0,L(0))], 20, Penrose.in_bounds, (w=5, h=5));
p = @time substitute(penrose(), Dict([hkite(0,0,L(0))]), 20, Penrose.in_bounds, (w=500, h=500));
@time SubstitutionTilings.CoreDefs.empirical_frequency2([hkite()], p)
@time empirical_frequency([hkite()], Dict(p))
for i in 1:10
    p = substitute(penrose(), Dict([hkite(0,0,L(0))]), 20, Penrose.in_bounds, (w=i*50, h=i*50));
    @time empirical_frequency(pentagon, p)
end
length(p)

@testset "Frequencies" begin
    @test Penrose.frequency([hkite()], 4) == ψ
    @test Penrose.frequency([hdart()], 4) == 1-ψ
    @test Penrose.frequency(pentagon, 4) == 7*ψ - 4
    @test Penrose.frequency(pentagon, 5) == 7*ψ - 4
    @test Penrose.frequency(pentagon, 6) == 7*ψ - 4
end

