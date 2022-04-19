module Penrose

export ζ, ϕ, Qζ, hkite, hdart,  penrose, PenroseElem, forced_penrose, conj, galois_gen

using Luxor
using StructEquality

using ...CoreDefs
using ...NumFields

@NumFields.simple_number_field_concrete Qζ [-1, 1, -1, 1] ζ
promote_rule(::Type{Qζ}, ::Type{Rational{Int64}}) = Qζ

const ϕ = ζ + ζ^9
const L = Qζ

function conj(x :: Qζ)
    return embed_field(a -> a, ζ^9, x)
end

function galois_gen(x :: Qζ)
    return embed_field(a -> a, ζ^3, x)
end

function Base.://(x :: Qζ, y :: Qζ)
    y1 = galois_gen(y)
    y2 = galois_gen(y1)
    y3 = galois_gen(y2)
    alg_conj = y1*y2*y3
    norm = y*alg_conj
    return x*alg_conj*norm.denom//norm.coeffs[1]
end

@def_structequal struct PenroseElem <: DGroupElem
    rot :: Int
    refl :: Bool
    z :: Qζ
end

function Base.:*(g :: PenroseElem, h :: PenroseElem)
    if !g.refl
        return PenroseElem(
            mod(g.rot + h.rot, 10),
            h.refl,
            g.z + ζ^(10-g.rot) * h.z
        )
    else
        return PenroseElem(
            mod(g.rot - h.rot, 10),
            !h.refl,
            g.z + ζ^(10-g.rot) * conj(h.z)
        )
    end
end
function Base.:*(g :: PenroseElem, x :: Qζ)
    if !g.refl
        return g.z + ζ^(10-g.rot) * x
    else
        return g.z + ζ^(10-g.rot) * conj(x)
    end
end
function CoreDefs.dilate(λ :: Qζ, g :: PenroseElem)
    return PenroseElem(g.rot, g.refl, λ*g.z)
end
function CoreDefs.id(::Type{PenroseElem})
    return PenroseElem(0,0,L(0))
end
function Base.inv(p::PenroseElem)
    return PenroseElem(0,p.refl,L(0))*PenroseElem(10-p.rot,0,L(0))*PenroseElem(0,0,-p.z)
end


@enum PenrosePTile begin
    Hkite=1
    Hdart=2
end
function hkite(r,s,z)
    return PenroseElem(r,s,z) => Hkite
end
function hkite()
    return hkite(0,0,L(0))
end
function hdart(r,s,z)
    return PenroseElem(r,s,z) => Hdart
end
function hdart()
    return hdart(0,0,L(0))
end



function embed_nf(x)
    ζ_float = complex(cos(2*pi/10), sin(2*pi/10))
    return embed_field(Complex{Float64}, ζ_float, x)
end
function embed_nf_p(x)
    z = embed_nf(x)
    return Point(real(z), imag(z))
end
function CoreDefs.embed_aff(g :: PenroseElem)
    m = [cos(2*pi/10) sin(2*pi/10); -sin(2*pi/10) cos(2*pi/10)]^g.rot

    z = embed_nf(g.z)
    return [m[1,1], m[2,1], (-1)^g.refl*m[1,2], (-1)^g.refl*m[2,2], real(z), imag(z)]
end



function CoreDefs.draw(ptile::PenrosePTile, action)
    if ptile == Hkite
        Luxor.poly([
            embed_nf_p(L(1)),
            embed_nf_p(L(ζ^4)),
            embed_nf_p(L(ζ^6)),
        ], close = true, action)
    else
        Luxor.poly([
            embed_nf_p(L(1-ζ^0//ϕ)),
            embed_nf_p(L(ζ^2-ζ^0//ϕ)),
            embed_nf_p(L(ζ^8-ζ^0//ϕ))
        ], close = true, action)
    end
end

function in_interval(x :: Qζ)
    x_float = embed_nf(x)
    return 0 ≤ real(x_float) && real(x_float) ≤ 1 && imag(x_float) ≈ 0
end
function in_border(x, ptile :: PenrosePTile)
    if ptile == Hkite
        return any([
            in_interval( (x - 1)//(ζ^4 - 1)),
            in_interval( (x - ζ^4)//(ζ^6 - ζ^4)),
            in_interval( (x - ζ^6)//(1 - ζ^6)),
        ])
    else
        return any([ # TODO fix!
            in_interval( (x - 1+ζ^0//ϕ)//(ζ^2 - 1)),
            in_interval( (x - ζ^2+ζ^0//ϕ)//(ζ^8 - ζ^2)),
            in_interval( (x - ζ^8+ζ^0//ϕ)//(1 - ζ^8)),
        ])
    end
end


function penrose()
    # This is a function only so that Revise works
    pen_subst = Dict([
        (Hkite, [
            hkite(7, 0, ζ^6+ζ^3//ϕ),
            hkite(6, 1, ζ^4+ζ^9//ϕ^2),
            hdart(7, 0, ζ//ϕ+ζ^8//ϕ^2),
            ]),
        (Hdart, [
            hdart(4,0, ζ^3//ϕ),
            hkite(3,1,ζ^6-ζ^5-ζ^0//ϕ^2)
        ])
    ])
    return SubSystem(pen_subst, ϕ) :: CoreDefs.SubSystem{PenroseElem, Qζ, PenrosePTile}
end

function color(x)
    l = x[2] == Hdart ? 0 : 1
    return 1 + l + 2*x[1].refl
end

function in_bounds(tile, n, window)
    phi = (1 + sqrt(5))/2
    center = embed_nf(tile[1].z)
    return abs(real(center))*phi^(n+1) < window.w + 3/2*phi^(n+1) && abs(imag(center))*phi^(n+1) < window.h + 3/2*phi^(n+1)
end



function force(tiling)
    fkite = [
        hkite(0,0,L(0))
        hkite(9, 1, ζ^8//ϕ)
    ]

    fdart = unique([
        hdart(0,0,L(0));
        hdart(7,1,L(0));
        PenroseElem(9, 1, ζ^2)*fkite;
        PenroseElem(8, 0, ζ)*fkite
    ])

    result = convert(Vector{Pair{PenroseElem, PenrosePTile}}, [])

    for tile in tiling
        if tile[2] == Hdart
            append!(result, tile[1]*fdart)
        else
            append!(result, tile[1]*fkite)
        end
    end
    return result
end

function forced_penrose()
    forced_subst = Dict([
        (label, force(penrose().sub[label]) ) for label in [Hkite, Hdart]
    ])

    return SubSystem(forced_subst, ϕ) :: CoreDefs.SubSystem{PenroseElem, Qζ, PenrosePTile}
end

@NumFields.simple_number_field_concrete Qψ [1,-1] ψ

function embed_psi(x)
    psi_float = (sqrt(5) - 1)/2
    return embed_field(rational_to_float, psi_float, x)
end
# ψ = 1/ϕ, except that ψ is an element of the field Q(ψ) and ϕ is an element of the field Q(ζ_10) 

function frequency(patch :: AbstractVector, depth)
    return frequency(Dict(patch), depth)
end
function frequency(patch :: AbstractDict, depth)
    freq = Qψ(0)
    harmonic = Dict([
        (Hkite, (ψ^(2*depth - 1))),
        (Hdart, (ψ^(2*depth)))
    ])

    center_ptile = patch[id(PenroseElem)]

    for label in [Hkite, Hdart]
        domain = substitute(penrose(), [(PenroseElem(0,0,L(0)) => label)], depth-1)
        forced_domain = Dict(substitute(forced_penrose(), ([PenroseElem(0,0,L(0)) => label]), depth-1))
        for tile in domain
            if tile[2] == center_ptile
                translated_patch = tile[1] * patch
                if translated_patch ⊆ forced_domain
                    freq += harmonic[label]
                end
            end
        end
    end
    return freq
end

end