module Pinwheel
__precompile__(false)

export i, sq5, wheel, pinwheel

using Hecke
using Luxor
using StructEquality

using ...CoreDefs


Qx, x = QQ["x"]
K_ns, a = NumberField([x^2 - 5, x^2 + 1], ["√5", "i"])
K, nsmap = simple_extension(K_ns)
sq5 = inv(Pinwheel.nsmap)(a[1])
i = inv(Pinwheel.nsmap)(a[2])

ζ = (i - 2)//sq5
conj = complex_conjugation(K)

@def_structequal struct PinwheelElem <: DGroupElem
    rot_ζ :: Int
    rot_i :: Int
    refl :: Bool
    z :: nf_elem
end

# Hack because the nf_elem hash doesn't work for some reason
function Base.hash(x :: PinwheelElem, h::UInt)
    return hash((x.rot_ζ, x.rot_i, x.refl, string(x.z)), h)
end

function embed_float(x  :: PinwheelElem)
    cpx = Hecke.complex_embeddings(K)[2](x.z, 64)
    return Float64.([real(cpx); imag(cpx)])
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
    return PinwheelElem(0,0,0,K(0))
end
function Base.:*(g :: PinwheelElem, x :: nf_elem)
    if !g.refl
        return g.z + i^g.rot_i * (-ζ)^-g.rot_ζ * x
    else
        return g.z + i^g.rot_i * (-ζ)^-g.rot_ζ * conj(x)
    end
end
function CoreDefs.dilate(λ :: nf_elem, g :: PinwheelElem)
    return PinwheelElem(g.rot_ζ, g.rot_i, g.refl, λ*g.z)
end
function CoreDefs.id(::PinwheelElem)
    return PinwheelElem(0,0,0,K(0))
end
function Base.inv(p::PinwheelElem)
    return PinwheelElem(0,0,p.refl,K(0))*PinwheelElem(-p.rot_ζ,-p.rot_i,0,K(0))*PinwheelElem(0,0,0,-p.z)
end
function CoreDefs.embed_aff(g :: PinwheelElem)
    t = -atan(1/2)
    m = [0 -1; 1 0]^(g.rot_i)*[cos(t) sin(t); -sin(t) cos(t)]^g.rot_ζ

    z = convert(Complex{Float64}, evaluate(g.z, infinite_places(K)[2]))
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
    return wheel(0,0,0,K(0))
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
const ϑ = (0,1)


function embed_float(x  :: Pair{PinwheelElem, PinwheelPTile})
    embed_float(x[1])
end
function vertices(t :: Pair{PinwheelElem, PinwheelPTile})
    return Dict([t[1]*K(0) => right, t[1]*K(1) => 2 .* right, t[1]*K(2) => ϑ, t[1]*(i) => right.-ϑ])
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