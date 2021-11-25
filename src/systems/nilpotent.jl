module Nilpotent
__precompile__(false)

export He, BHPTile, a, b, bhp_subst

using Nemo
using StructEquality
using Base.Iterators

QQ = Rational{Int}

using ...CoreDefs
using ...Collaring


@def_structequal struct He <: EGroupElem
    a::QQ
    b::QQ
    c::QQ
end

@enum BHPTile begin
    A=1
    B=2
end



function Base.:*(x :: He, y :: He)
    return He( x.a + y.a, x.b + y.b, x.c + y.c + (x.a*y.b-y.a*x.b)//2)
end

function CoreDefs.dilate(λ, x :: He)
    return He(λ*x.a, λ * x.b, λ^2 * x.c )
end

function Base.inv(x :: He)
    return He(-x.a, -x.b, -x.c)
end

function CoreDefs.id(::He)
    return He(0,0,0)
end

function a(a,b,c)
    return (He(a,b,c), A)
end
function a()
    return a(0,0,0)
end
function b(a,b,c)
    return (He(a,b,c), B)
end
function b()
    return b(0,0,0)
end

function bhp_subst()
    heisenberg_subst = Dict([
        (A, [
            a(-2,-2,-8),
            a(0,-2,-8),
            a(2,-2,-8),
            a(-2,0,-8),
            a(0,0,-8),
            a(2,0,-8),
            a(-2,2,-8),
            a(0,2,-8),
            a(2,2,-8),
            a(-2,-2,-6),
            b(0,-2,-6),
            b(2,-2,-6),
            a(-2,0,-6),
            a(0,0,-6),
            a(2,0,-6),
            b(-2,2,-6),
            a(0,2,-6),
            b(2,2,-6),
            b(-2,-2,-4),
            a(0,-2,-4),
            a(2,-2,-4),
            b(-2,0,-4),
            a(0,0,-4),
            b(2,0,-4),
            a(-2,2,-4),
            a(0,2,-4),
            a(2,2,-4),
            a(-2,-2,-2),
            a(0,-2,-2),
            a(2,-2,-2),
            a(-2,0,-2),
            a(0,0,-2),
            a(2,0,-2),
            a(-2,2,-2),
            a(0,2,-2),
            a(2,2,-2),
            b(-2,-2,0),
            b(0,-2,0),
            b(2,-2,0),
            b(-2,0,0),
            a(0,0,0),
            b(2,0,0),
            b(-2,2,0),
            b(0,2,0),
            b(2,2,0),
            a(-2,-2,2),
            a(0,-2,2),
            a(2,-2,2),
            a(-2,0,2),
            a(0,0,2),
            a(2,0,2),
            a(-2,2,2),
            a(0,2,2),
            a(2,2,2),
            a(-2,-2,4),
            a(0,-2,4),
            b(2,-2,4),
            a(-2,0,4),
            a(0,0,4),
            a(2,0,4),
            a(-2,2,4),
            a(0,2,4),
            b(2,2,4),
            b(-2,-2,6),
            b(0,-2,6),
            a(2,-2,6),
            b(-2,0,6),
            a(0,0,6),
            b(2,0,6),
            b(-2,2,6),
            a(0,2,6),
            a(2,2,6),
            a(-2,-2,8),
            a(0,-2,8),
            a(2,-2,8),
            a(-2,0,8),
            a(0,0,8),
            a(2,0,8),
            a(-2,2,8),
            a(0,2,8),
            a(2,2,8),
        ]),
        (B, [
            a(-2,-2,-8),
            a(0,-2,-8),
            a(2,-2,-8),
            a(-2,0,-8),
            a(0,0,-8),
            a(2,0,-8),
            a(-2,2,-8),
            a(0,2,-8),
            a(2,2,-8),
            a(-2,-2,-6),
            b(0,-2,-6),
            b(2,-2,-6),
            b(-2,0,-6),
            a(0,0,-6),
            a(2,0,-6),
            b(-2,2,-6),
            a(0,2,-6),
            b(2,2,-6),
            b(-2,-2,-4),
            a(0,-2,-4),
            a(2,-2,-4),
            a(-2,0,-4),
            a(0,0,-4),
            b(2,0,-4),
            a(-2,2,-4),
            a(0,2,-4),
            a(2,2,-4),
            a(-2,-2,-2),
            a(0,-2,-2),
            a(2,-2,-2),
            a(-2,0,-2),
            a(0,0,-2),
            a(2,0,-2),
            a(-2,2,-2),
            a(0,2,-2),
            a(2,2,-2),
            b(-2,-2,0),
            b(0,-2,0),
            b(2,-2,0),
            b(-2,0,0),
            a(0,0,0),
            b(2,0,0),
            b(-2,2,0),
            b(0,2,0),
            b(2,2,0),
            a(-2,-2,2),
            a(0,-2,2),
            a(2,-2,2),
            a(-2,0,2),
            a(0,0,2),
            a(2,0,2),
            a(-2,2,2),
            a(0,2,2),
            a(2,2,2),
            a(-2,-2,4),
            b(0,-2,4),
            a(2,-2,4),
            a(-2,0,4),
            a(0,0,4),
            a(2,0,4),
            b(-2,2,4),
            a(0,2,4),
            b(2,2,4),
            b(-2,-2,6),
            a(0,-2,6),
            b(2,-2,6),
            b(-2,0,6),
            a(0,0,6),
            b(2,0,6),
            a(-2,2,6),
            a(0,2,6),
            a(2,2,6),
            a(-2,-2,8),
            a(0,-2,8),
            a(2,-2,8),
            a(-2,0,8),
            a(0,0,8),
            a(2,0,8),
            a(-2,2,8),
            a(0,2,8),
            a(2,2,8),
        ]),
    ])

    return SubSystem(heisenberg_subst, 3//1) 
end

function Collaring.collar_in(tiling, t ::He)
    collar_shape = t*[He(x,y,z) for x=[-2,0,2] for y=[-2,0,2] for z=[-2,0,2]]
    collar = []
    tiling_dict = Dict(tiling)
    for s in collar_shape
        if !haskey(tiling_dict, s)
            throw(Collaring.UnrecognizedCollar)
        end
        push!(collar, (s,tiling_dict[s]))
    end
    return collar
end

end