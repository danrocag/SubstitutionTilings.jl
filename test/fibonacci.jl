using SubstitutionTilings
using SubstitutionTilings.NumFields
using SubstitutionTilings.CoreDefs
using StructEquality


@NumFields.simple_number_field_concrete Qτ [1,1] τ
Base.promote_rule(::Type{Qτ}, ::Type{<:Integer}) = Qτ

@def_structequal struct FibElem <: DGroupElem
    a :: Qτ
end

@enum FibTile begin
    A=1
    B=2
end

function Base.:*(x :: FibElem, y :: FibElem)
    return FibElem(x.a + y.a)
end

function CoreDefs.dilate(λ :: Qτ, x :: FibElem)
    return FibElem(λ*x.a)
end

function Base.inv(x :: FibElem)
    return FibElem(-x.a)
end

function CoreDefs.id(::Type{FibElem})
    return FibElem(0)
end

function CoreDefs.is_interior(tiling :: Dict, t :: FibElem)
    return haskey(tiling, t) && (haskey(tiling, FibElem(t.a-τ)) || haskey(tiling, FibElem(t.a-(1+τ)//2))) && (haskey(tiling, FibElem(t.a+τ)) || haskey(tiling, FibElem(t.a+(1+τ)//2)))
end

function CoreDefs.collar_in(tiling :: Dict, t :: FibElem)
    if !is_interior(tiling, t)
        throw(UnrecognizedCollar)
    end
    collar = Dict{FibElem, FibTile}([])
    shifts = ([0, -τ, -(1+τ)//2, τ, (1+τ)//2])
    for shift in shifts
        if haskey(tiling, FibElem(t.a+shift))
            collar[FibElem(t.a+shift)] = tiling[FibElem(t.a+shift)]
        end
    end
    return collar
end


fib = SubSystem(Dict(A => [FibElem(Qτ(-1)//2) => A, FibElem(τ//2) => B], B => [FibElem(0) => A]),τ)

fib_tiling = substitute(fib, Dict([FibElem(0) => A]), 3)
initial_collar = collar_in(fib_tiling, FibElem(0))

(collars, fib_c) = total_collaring(fib, initial_collar)


frequency(fib, initial_collar, collars[4], 4)