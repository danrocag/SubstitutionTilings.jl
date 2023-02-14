module Pinwheel

export i, sq5, wheel, pinwheel

using Luxor
using StructEquality

using ...CoreDefs
using ...NumFields

@NumFields.simple_number_field_concrete Qθ [-36,0,8,0] θ
Base.promote_rule(::Type{Qθ}, ::Type{<:Integer}) = Qθ

const i = (θ^3-2θ)/12
const sq5 = θ-i

function conj(x :: Qθ)
    return embed_field(a -> a, sq5-i, x)
end

function rconj(x :: Qθ)
    return embed_field(a -> a, -sq5+i, x)
end

function Base.://(x :: Qθ, y :: Qθ)
    y1 = conj(y)
    y2 = rconj(y)
    y3 = rconj(y1)
    alg_conj = y1*y2*y3
    norm = y*alg_conj
    return x*alg_conj*norm.denom//norm.coeffs[1]
end
function Base.:^(x :: Qθ, k :: Integer)
    if k == 0
        Qθ(1)
    elseif k < 0
        Qθ(1)//(x^(-k))
    elseif k%2 == 0
        y = x^(div(k,2))
        y*y
    else
        x*x^(k-1)
    end
end

ζ = (i-2)//sq5

@def_structequal struct PinwheelElem <: DGroupElem
    rot_ζ :: Int
    rot_i :: Int
    refl :: Bool
    z :: Qθ
end

function embed_float(x  :: Qθ)
    embed_field(Complex{Float64}, im + sqrt(5), x)
end

function Base.:*(g :: PinwheelElem, h :: PinwheelElem)
    if !g.refl
        return PinwheelElem(
            g.rot_ζ + h.rot_ζ,
            mod(g.rot_i + h.rot_i, 4),
            h.refl,
            g.z + i^g.rot_i * (-ζ)^-g.rot_ζ * h.z
        )
    else
        return PinwheelElem(
            g.rot_ζ - h.rot_ζ,
            mod(g.rot_i - h.rot_i, 4),
            !h.refl,
            g.z + i^g.rot_i * (-ζ)^-g.rot_ζ * conj(h.z)
        )
    end
end
function CoreDefs.id(::Type{PinwheelElem})
    return PinwheelElem(0,0,0,Qθ(0))
end
function Base.:*(g :: PinwheelElem, x :: Qθ)
    if !g.refl
        return g.z + i^g.rot_i * (-ζ)^-g.rot_ζ * x
    else
        return g.z + i^g.rot_i * (-ζ)^-g.rot_ζ * conj(x)
    end
end
function CoreDefs.dilate(λ :: Qθ, g :: PinwheelElem)
    return PinwheelElem(g.rot_ζ, g.rot_i, g.refl, λ*g.z)
end
function CoreDefs.id(::PinwheelElem)
    return PinwheelElem(0,0,0,Qθ(0))
end
function Base.inv(p::PinwheelElem)
    return PinwheelElem(0,0,p.refl,Qθ(0))*PinwheelElem(-p.rot_ζ,-p.rot_i,0,Qθ(0))*PinwheelElem(0,0,0,-p.z)
end
function CoreDefs.embed_aff(g :: PinwheelElem)
    t = -atan(1/2)
    m = [0 -1; 1 0]^(g.rot_i)*[cos(t) sin(t); -sin(t) cos(t)]^g.rot_ζ

    z = embed_float(g.z)
    return [m[1,1], m[2,1], (-1)^g.refl*m[1,2], (-1)^g.refl*m[2,2], real(z), imag(z)]
end

struct PinwheelPTile end

function CoreDefs.draw(::PinwheelPTile, action)
    Luxor.poly([
        Point(0,0),
        Point(2,0),
        Point(0,1),
    ], close = true, action)
end

function wheel(rot_ζ, rot_i, refl, z)
    return PinwheelElem(rot_ζ, rot_i, refl, z) => PinwheelPTile()
end
function wheel()
    return wheel(0,0,0,Qθ(0))
end

function color(p)
    return 1+Int(p[1].refl)
end


function pinwheel()
    sub  = Dict([
        (PinwheelPTile(), [
            wheel(-1,3,1,i*sq5-ζ)
            wheel(-1,0,1,sq5+2*ζ)
            wheel(-1,0,0,sq5+2*ζ)
            wheel(-1,0,1,2*sq5+2*ζ)
            wheel(-1,2,0,sq5*i-3*ζ)
        ])
    ])
    return SubSystem(sub, sq5)
end

const right = (1,0)
const down = (0,1)


function embed_float(x  :: Pair{PinwheelElem, PinwheelPTile})
    embed_float(x[1])
end
function vertices(t :: Pair{PinwheelElem, PinwheelPTile})
    return Dict([t[1]*Qθ(0) => right, t[1]*Qθ(1) => 2 .* right, t[1]*Qθ(2) => down, t[1]*(i) => right.-down])
end

function precollar_in(tiling :: Dict, g :: PinwheelElem)
    gvertices = vertices(g => tiling[g])
    precollar =  []
    for t in tiling
        if !isempty(keys(vertices(t)) ∩ keys(gvertices))
            push!(precollar,t)
        end
    end
    return Dict(precollar)
end

function CoreDefs.is_interior(tiling :: Dict, g :: PinwheelElem)
    angle_sums = Dict(g => (0,0) for k in keys(vertices(g => tiling[g])))
    for tile in tiling
        for (vertex,angle) in vertices(tile)
            if vertex ∈ keys(angle_sums)
                angle_sums[vertex] = angle_sums[vertex] .+ angle
            end
        end
    end

    for (vertex, angle) in angle_sums
        if (angle[1] % 4 ≠ 0) || angle[2] ≠ 0
            return false
        end
    end
    return true
end

function CoreDefs.collar_in(tiling :: Dict, g :: PinwheelElem)
    precollar =  precollar_in(tiling , g)
    if !is_interior(precollar, g)
        throw(UnrecognizedCollar)
    else
        return precollar
    end
end

end