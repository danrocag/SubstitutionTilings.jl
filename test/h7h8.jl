using SubstitutionTilings
using SubstitutionTilings.NumFields
using Luxor

@NumFields.simple_number_field_concrete Qθ  [-4,0,16,0] θ
promote_rule(::Type{Qθ}, ::Type{Rational{Int64}}) = Qθ

function embed_float(x::Qθ)
    embed_field(Float64, sqrt(3) + sqrt(5), x)
end

function embed_float(v::Vector{Qθ})
    return embed_float.(v)
end
function Luxor.translate(v::Vector{Qθ})
    translate(Point(Tuple(embed_float(v))))
end

sq3 = (θ^3-14*θ)//4
sq5 = θ-sq3
ϕ = (1 + sq5)//2

rot = Qθ[1//2 sq3//2; -sq3//2 1//2]

hex1 = Qθ[1; 0]
hex2 = Qθ[1//2; sq3//2]

kite_outline = [
    Qθ[0,0],
    Qθ[0;sq3//2],
    Qθ[1//2; sq3//2],
    Qθ[3//4;sq3//4],
]

hat_kites = [
    [0,0,5],[0,0,0],[0,0,1],[0,0,2],
    [1,1,2], [1,1,3],
    [2,-1,4], [2,-1,5]
]

function draw_kite(;action = :fillstroke)
    poly(Point.(Tuple.(embed_float.(kite_outline))), action=action, close=true)    
end

function draw_hat(;action=:fillstroke)
    for (x,y,k) in hat_kites
        @layer begin
            translate(x*hex1+y*hex2)
            rotate(-k*2π/6)
            draw_kite(action=action)
        end
    end
end

function hat_centroid(x,y,r)
    x*hex1 + y*hex2 + rot^r*Qθ[3//4; 0]
end

h7 = [
    (0,0,1,1)
    (2,-1,2,0)
    (-1,2,0,0)
    (-2,1,1,0)
    (-1,-1,2,0)
    (1,-2,4,0)
    (0,-3,3,0)
]
h8 = [
    h7
    (-2,-2,2,0)
]

meta_list = Dict([:h7 => h7, :h8 => h8])

function draw_meta(t; action=:fillstroke)
    for (x,y,k,s) = meta_list[t]
        @layer begin
            translate(x*hex1+y*hex2)
            rotate(k*2π/6)
            if s == 1
                transform([-1 0 0 1 0 0])
            end
            draw_hat(action=action)
        end
    end
end

h77 = [
    (0,0,0)
]
h78 = [
    (-4,-1,0)
    (-1,5,2)
    (-5,4,1)
    (0,-6,4)
    (-4,-7,5)
]
h87=h77
h88 = [
    h78
    (-9,-3,0)
]

h = Dict([
    (:h7,:h7) => h77,
    (:h7,:h8) => h78,
    (:h8,:h7) => h87,
    (:h8,:h8) => h88,
])

shifts = Dict([
    :h7 => Point(20,-1.7),
    :h8 => Point(1.8,-0.2)
    ])
sc=0.37196
function inflate(i, f=draw_meta)
    for j=[:h7,:h8]
        for (x,y,k)=h[(i,j)]
            @layer begin
                translate(x*hex1+y*hex2)
                translate(shifts[i])
                rotate(-k*2π/6)
                f(j)
            end
        end
    end
end

@draw begin
    setstrokescale(false)
    circle(O, 5, :fill)
    scale(60)
    sethue("blueviolet")
    inflate(:h8, j -> @layer begin
        scale(sc)
        inflate(j, k -> draw_meta(k, action=:stroke))
    end)
end 1200 1200
phi^(-2)