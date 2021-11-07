module Pinwheel

export i, sq5, wheel, pinwheel

using Hecke
using Luxor
using StructEquality

using ...CoreDefs

Qx, x = QQ["x"]
K_ns, a = NumberField([x^2 - 5, x^2 + 1], ["√5", "i"])
K, nsmap = simple_extension(K_ns)
sq5 = -inv(Pinwheel.nsmap)(a[1])
i = inv(Pinwheel.nsmap)(a[2])

ζ = (2 - i)//sq5

conj = complex_conjugation(K)

function eval(z)
    return convert(Complex{Float64}, evaluate(z, infinite_places(K)[1], 1000))
end

@def_structequal struct PinwheelElem <: EGroupElem
    rotI :: Int
    rot4 :: Int
    refl :: Bool
    z :: nf_elem
end


function Base.:*(g :: PinwheelElem, h :: PinwheelElem)
    if !g.refl
        return PinwheelElem(
            g.rotI + h.rotI,
            mod(g.rot4 + h.rot4,4),
            h.refl,
            g.z + i^g.rot4 * ζ^g.rotI * h.z
        )
    else
        return PinwheelElem(
            g.rotI - h.rotI,
            mod(g.rot4 - h.rot4,4),
            !h.refl,
            g.z + i^g.rot4 * ζ^g.rotI * conj(h.z)
        )
    end
end
function Base.:*(g :: PinwheelElem, x :: nf_elem)
    if !g.refl
        return g.z + i^g.rot4 * ζ^g.rotI * x
    else
        return g.z + i^g.rot4 * ζ^g.rotI * conj(x)
    end
end
function CoreDefs.dilate(λ :: nf_elem, g :: PinwheelElem)
    return PinwheelElem(g.rotI, g.rot4, g.refl, λ*g.z)
end
function CoreDefs.id(::PinwheelElem)
    return PinwheelElem(0,0,0,K(0))
end
function Base.inv(p::PinwheelElem)
    return PinwheelElem(0,0,p.refl,K(0))*PinwheelElem(-p.rotI, -p.rot4, 0, K(0))*PinwheelElem(0,0,0,-p.z)
end
function CoreDefs.embed_aff(g :: PinwheelElem)
    t = -atan(1/2)
    m = [cos(t) -sin(t); sin(t) cos(t)]^g.rotI * [0 -1; 1 0]^g.rot4 

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

function wheel(rotI, rot4, refl, z)
    return (PinwheelElem(rotI, rot4, refl, z), PinwheelPTile())
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
            wheel(1,3,1,ζ^-1+i*3//sq5),
            wheel(1,0,1,(2*i+1)//sq5),
            wheel(1,0,0,(2*i+1)//sq5),
            wheel(1,2,0,sq5+(2*i+1)//sq5),
            wheel(1,0,1,sq5+(2*i+1)//sq5),
        ])
    ])
    return SubSystem(sub, sq5)
end


function CoreDefs.vertices(::PinwheelPTile)
    return [K(0), K(1), K(2), i]
end

function in_interval(x :: nf_elem)
    x_float = eval(x)
    return 0 ≤ real(x_float) && real(x_float) ≤ 1 && abs(imag(x_float)) < eps()
end
function CoreDefs.in_border(x, :: PinwheelPTile)
    return any([
            in_interval( x//2),
            in_interval( (x - 2)//(i - 2)),
            in_interval( (x - i)//(0 - i)),
        ])
end


end