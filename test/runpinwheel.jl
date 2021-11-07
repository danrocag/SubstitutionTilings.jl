using SubstitutionTilings
using SubstitutionTilings.Pinwheel
using Test

using Luxor
using Hecke

using Profile

K = Pinwheel.K



@time collars = SubstitutionTilings.CoreDefs.partial_collars(pinwheel(), 0, 6);
length(collars)


w = 800
h = 800
sc = 60
@draw begin
    colors = ["#DD93FC", "#E7977A",]
    """
    first_tile = wheel(0,0,0,-i//2 - 1//4)
    tiling = substitute(pinwheel(), [first_tile],3)
    filter!(
        tile -> any(v -> in_border(v//sq5^3, first_tile), vertices(tile)),
        tiling)
    """
    #"""
    first_tile = wheel()
    tiling  = collars[4]
    #"""
    #display(map(vertices, tiling))
    setline(1)

    for tile in tiling
        origin()
        draw(tile, sc, colors[Pinwheel.color(tile)], :fill)
        draw(tile, sc, "black", :stroke)
    end
    setline(2)
    draw(first_tile, sc*sqrt(5)^0, "black", :stroke)
end w h




pentagon = [
    (PenroseElem(k+s,s,L(1))*PenroseElem(0,0,L(-1)),Penrose.Hkite)
    for k=0:2:10 for s=0:1
];
p = @time substitute(penrose(), [hkite(0,0,L(0))], 20, Penrose.in_bounds, (w=100, h=100));
@time empirical_frequency(pentagon, p)


@testset "Frequencies" begin
    @test Penrose.frequency([hkite()], 4) == ψ
    @test Penrose.frequency([hdart()], 4) == 1-ψ
    @test Penrose.frequency(pentagon, 4) == 7*ψ - 4
    @test Penrose.frequency(pentagon, 5) == 7*ψ - 4
    @test Penrose.frequency(pentagon, 6) == 7*ψ - 4
end

