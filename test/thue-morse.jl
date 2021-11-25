using SubstitutionTilings
using SubstitutionTilings.CoreDefs
using SubstitutionTilings.Collaring
using StructEquality
using LinearAlgebra

@def_structequal struct TM <: EGroupElem
    a::Int
end

@enum TMTile begin
    A=1
    B=2
end

function Base.:*(x :: TM, y :: TM)
    return TM(x.a + y.a)
end

function CoreDefs.dilate(λ, x :: TM)
    return TM(λ*x.a)
end

function Base.inv(x :: TM)
    return TM(-x.a)
end

function CoreDefs.id(::TM)
    return TM(0)
end

function Collaring.collar_in(tiling, t :: TM)
    collar_shape = t*[TM(x) for x in [-2,0,2]]
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

thue_morse = SubSystem(Dict(A => [(TM(-1), A), (TM(1),B)], B => [(TM(-1), B), (TM(1),A)]),2)

substitute(thue_morse, [(TM(0), A)], 4)

(collars, S) = Collaring.accessible_subst(thue_morse, [(TM(-2), A), (TM(0), A), (TM(2), B)])
length(collars)
A_thue = transition_matrix(S, 1:6)
(e, V) = eigen(A_thue)

function 