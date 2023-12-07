module AmmannBeenker

export ζ, δ, sq2, Qζ, rhomb, square, ammannbeenker, ABElem, color

using Luxor
using StructEquality

using ...CoreDefs
using ...NumFields

@NumFields.simple_number_field_concrete Qζ [-1, 0, 0, 0] ζ
promote_rule(::Type{Qζ}, ::Type{Rational{Int64}}) = Qζ

const sq2 = ζ + ζ^7
const δ = 1 + sq2
const L = Qζ

const ζ_powers = [ζ^i for i=0:8]
function conj(x :: Qζ)
    return embed_field(a -> a, ζ_powers[8], x)
end

function galois_gen(x :: Qζ)
    return embed_field(a -> a, ζ_powers[4], x)
end

function Base.://(x :: Qζ, y :: Qζ)
    y1 = galois_gen(y)
    y2 = galois_gen(y1)
    y3 = galois_gen(y2)
    alg_conj = y1*y2*y3
    norm = y*alg_conj
    return x*alg_conj*norm.denom//norm.coeffs[1]
end

@def_structequal struct ABElem <: DGroupElem
    rot :: Int
    refl :: Bool
    z :: Qζ
end

function Base.:*(g :: ABElem, h :: ABElem)
    if !g.refl
        return ABElem(
            mod(g.rot + h.rot, 8),
            h.refl,
            g.z + ζ_powers[g.rot+1] * h.z
        )
    else
        return ABElem(
            mod(g.rot - h.rot, 8),
            !h.refl,
            g.z + ζ_powers[g.rot+1] * conj(h.z)
        )
    end
end
function Base.:*(g :: ABElem, x :: Qζ)
    if !g.refl
        return g.z + ζ_powers[g.rot+1] * x
    else
        return g.z + ζ_powers[g.rot+1] * conj(x)
    end
end
function CoreDefs.dilate(λ :: Qζ, g :: ABElem)
    return ABElem(g.rot, g.refl, λ*g.z)
end
function CoreDefs.id(::Type{ABElem})
    return ABElem(0,false,L(0))
end
function Base.inv(p::ABElem)
    return ABElem(0,p.refl,L(0))*ABElem(8-p.rot,L(0))*ABElem(0,-p.z)
end


@enum ABPTile begin
    Rhomb=1
    Square=2
end
function rhomb(r,s,z)
    return ABElem(r,s,z) => Rhomb
end
function rhomb()
    return rhomb(0,false, L(0))
end
function square(r,s, z)
    return ABElem(r,s,z) => Square
end
function square()
    return square(0,false,L(0))
end


function embed_nf(x :: Qζ)
    ζ_float = complex(cos(2*pi/8), -sin(2*pi/8))
    return embed_field(Complex{Float64}, ζ_float, x)
end
function embed_nf_p(x :: Qζ)
    z = embed_nf(x)
    return Point(real(z), imag(z))
end
function embed_float(x :: ABElem)
    embed_nf(x.z)
end

function CoreDefs.embed_aff(g :: ABElem)
    m = [cos(2*pi/8) sin(2*pi/8); -sin(2*pi/8) cos(2*pi/8)]^g.rot

    z = embed_nf(g.z)
    return [m[1,1], m[2,1], (-1)^g.refl*m[1,2], (-1)^g.refl*m[2,2], real(z), imag(z)]
end



function CoreDefs.draw(ptile::ABPTile, action)
    if ptile == Square
        Luxor.poly([
            embed_nf_p(sq2),
            embed_nf_p(ζ^2*sq2//2+sq2//2),
            embed_nf_p(L(0)),
        ], close = true, action)
    else
        Luxor.poly([
            embed_nf_p(L(0)),
            embed_nf_p(L(1)),
            embed_nf_p(L(1+ζ^3)),
            embed_nf_p(L(ζ^3)),
        ], close = true, action)
    end
end

"""
function in_interval(x :: Qζ)
    x_float = embed_nf(x)
    return 0 ≤ real(x_float) && real(x_float) ≤ 1 && imag(x_float) ≈ 0
end
function in_border(x, ptile :: ABPTile)
    if ptile == Square
        return any([
            in_interval( (x - 1)//(ζ^2 - 1)),
            in_interval( (x - ζ^2)//(ζ^4 - ζ^2)),
            in_interval( (x - ζ^6)//(1 - ζ^6)),
            in_interval( (x - ζ^6)//(1 - ζ^6)),
        ])
    else
        return any([ # TODO fix!
        in_interval( (x - 0)//(ζ^2 - 1)),
        in_interval( (x - ζ^2 + 1)//(ζ^8 - ζ^2)),
        in_interval( (x - ζ^8 + 1)//(1 - ζ^8)),
        ])
    end
end
function in_border(x, tile :: Pair{ABElem, ABPTile})
    in_border(inv(t[1])*x, t[2])
end

function vertices(ptile :: ABPTile)
    Dict(ptile == Hkite ? Pair{Qζ, Int}[
        1 => 36,
        ζ^4 => 72,
        ζ^6 => 72] : [
        0 => 108,
        ζ^2-1 => 36,
        ζ^8-1 => 36])
end
function CoreDefs.vertices(tile :: Pair{ABElem, ABPTile})
    tile[1]*vertices(tile[2])
end
@collar_in_from_vertices ABElem ABPTile 0 360
"""


function ammannbeenker()
    subst = Dict([
        (Rhomb, [
            square(4,1,sq2),
            square(7,0,ζ^3*sq2),
            rhomb(2,0,ζ),
            rhomb(0,0,sq2),
            rhomb(0,0,ζ^3*sq2),
            square(0,1,ζ^2+ζ^3),
            square(3,0,L(1)+ζ),
            ]),
        (Square, [
            square(5,0, ζ*sq2),
            rhomb(6,0,ζ*sq2),
            square(3,0, ζ*(1+sq2)+L(1)-ζ^2),
            rhomb(0,0,δ),
            square(4,1,δ),
            ])])
    return SubSystem(subst, δ) :: CoreDefs.SubSystem{ABElem, Qζ, ABPTile}
end

function color(x)
    if x[2] == Rhomb
        1
    else
        x[1].refl ? 2 : 3
    end
end
"""
function in_bounds(tile, n, window)
    phi = (1 + sqrt(5))/2
    center = embed_nf(tile[1].z)
    return abs(real(center))*phi^(n+1) < window.w + 3/2*phi^(n+1) && abs(imag(center))*phi^(n+1) < window.h + 3/2*phi^(n+1)
end
function in_bounds_strict(tile, window)
    phi = (1 + sqrt(5))/2
    center = embed_nf(tile[1].z)
    return abs(real(center)) <= window.w/2 && abs(imag(center)) <= window.h/2
end



function force(tiling)
    fkite = [
        hkite(0,0,L(0))
        hkite(9, 1, ζ^8//ϕ)
    ]

    fdart = unique([
        hdart(0,0,L(0));
        hdart(7,1,L(0));
        ABElem(9, 1, ζ^2)*fkite;
        ABElem(8, 0, ζ)*fkite
    ])

    result = convert(Vector{Pair{ABElem, ABPTile}}, [])

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

    return SubSystem(forced_subst, ϕ) :: CoreDefs.SubSystem{ABElem, Qζ, ABPTile}
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

    center_ptile = patch[id(ABElem)]

    for label in [Hkite, Hdart]
        domain = substitute(penrose(), [(ABElem(0,0,L(0)) => label)], depth-1)
        forced_domain = Dict(substitute(forced_penrose(), ([ABElem(0,0,L(0)) => label]), depth-1))
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
"""
end