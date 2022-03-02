using SubstitutionTilings
using SubstitutionTilings.CoreDefs
using StructEquality
using LinearAlgebra

@def_structequal struct TM <: DGroupElem
    a::Int
end

@enum TMTile begin
    A=1
    B=2
end

function Base.:*(x :: TM, y :: TM)
    return TM(x.a + y.a)
end

function SubstitutionTilings.dilate(λ, x :: TM)
    return TM(λ*x.a)
end

function SubstitutionTilings.inv(x :: TM)
    return TM(-x.a)
end

function SubstitutionTilings.id(::Type{TM})
    return TM(0)
end

function SubstitutionTilings.is_interior(tiling :: Dict, t :: TM)
    return haskey(t.a, tiling) && haskey(t.a-1, tiling) && haskey(t.a+1, tiling)
end

function SubstitutionTilings.collar_in(tiling :: Dict, t :: TM)
    collar_shape = t*[TM(x) for x in [-2,0,2]]
    collar = []
    tiling_dict = tiling
    for s in collar_shape
        if !haskey(tiling_dict, s)
            throw(UnrecognizedCollar)
        end
        push!(collar, (s => tiling_dict[s]))
    end
    return Dict(collar)
end

thue_morse = SubSystem(Dict(A => [(TM(-1) => A), (TM(1) => B)], B => [(TM(-1) => B), (TM(1) => A)]),2)

substitute(thue_morse, [(TM(0) => A)], 4)

(collars, S) = total_collaring(thue_morse, Dict([(TM(-2) => A), (TM(0) => A), (TM(2) => B)]))
length(collars)