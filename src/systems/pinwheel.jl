module Pinwheel

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

ζ = (2 - i)//sq5

conj = complex_conjugation(K)

@def_structequal struct PinwheelElem <: EGroupElem
    rot :: Int
    refl :: Bool
    z :: nf_elem
end


function Base.:*(g :: PinwheelElem, h :: PinwheelElem)
    if !g.refl
        return PinwheelElem(
            g.rot + h.rot,
            h.refl,
            g.z + ζ^-g.rot * h.z
        )
    else
        return PinwheelElem(
            g.rot - h.rot,
            !h.refl,
            g.z + ζ^-g.rot * conj(h.z)
        )
    end
end
function Base.:*(g :: PinwheelElem, x :: nf_elem)
    if !g.refl
        return g.z + ζ^-g.rot * x
    else
        return g.z + ζ^-g.rot * conj(x)
    end
end
function CoreDefs.dilate(λ :: nf_elem, g :: PinwheelElem)
    return PinwheelElem(g.rot, g.refl, λ*g.z)
end
function CoreDefs.id(::PinwheelElem)
    return PinwheelElem(0,0,K(0))
end
function Base.inv(p::PinwheelElem)
    return PinwheelElem(0,p.refl,K(0))*PinwheelElem(-p.rot,0,K(0))*PinwheelElem(0,0,-p.z)
end
function CoreDefs.embed_aff(g :: PinwheelElem)
    t = -atan(1/2)
    m = [cos(t) sin(t); -sin(t) cos(t)]^g.rot

    z = convert(Complex{Float64}, evaluate(g.z, infinite_places(K)[1]))
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

function wheel(rot, refl, z)
    return (PinwheelElem(rot, refl, z), PinwheelPTile())
end
function wheel()
    return wheel(0,0,K(0))
end

function color(p)
    return 1+Int(p[1].refl)
end


function pinwheel()
    sub  = Dict([
        (PinwheelPTile(), [
            wheel(4,0,K(i+1))
        ])
    ])
    return SubSystem(sub, sq5)
end
end