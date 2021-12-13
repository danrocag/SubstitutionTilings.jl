using SubstitutionTilings
using SubstitutionTilings.CoreDefs
using SubstitutionTilings.NumFields
using StructEquality

@NumFields.simple_number_field Qτ [1,1] τ
Base.promote_rule(::Type{Qτ}, ::Type{<:Integer}) = Qτ

@def_structequal struct FibElem <: EGroupElem
    a :: Qτ
end

@enum FibTile begin
    A=1
    B=2
end

function Base.:*(x :: FibElem, y :: FibElem)
    return FibElem(x.a + y.a)
end

function CoreDefs.dilate(λ, x :: FibElem)
    return FibElem(λ*x.a)
end

function Base.inv(x :: FibElem)
    return FibElem(-x.a)
end

function CoreDefs.id(::FibElem)
    return FibElem(0)
end

fib = SubSystem(Dict(A => [FibElem(Qτ(-1)/2) => A, FibElem(τ/2) => B], B => [FibElem(0) => A]),τ)

fib_tiling = substitute(fib, Dict([FibElem(0) => A]),5)

for tile in (fib_tiling)
    (NumFields.reduce!(tile[1].a))
end