using SubstitutionTilings
using SubstitutionTilings.Heisenberg
using Luxor


width = 1280
height = 1280
sc = 10
@draw begin
    colors = ["#21ABC2", "#43848F", "#F6504F", "#C2216C"]
    tiling = He(0,0,0)*substitute(heis(), dilate(-1//7 ,Heisenberg.He(0,10,1))*hn(0),9);
    println(sizeof(tiling))
    
    setline(0.5)
    for tile in tiling
        draw(tile, sc, colors[Heisenberg.color(tile)], :fill)
        
        draw(tile, sc, "black", :stroke)
    end
end width height 