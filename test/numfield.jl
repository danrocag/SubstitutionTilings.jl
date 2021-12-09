using SubstitutionTilings.NumFields
using BenchmarkTools
using StaticArrays

@NumFields.simple_number_field Qζ [-1, 1, -1, 1] ζ
Base.promote_rule(::Type{Qζ}, ::Type{<:Integer}) = Qζ
@macroexpand @NumFields.simple_number_field Qζ ZZ[-1, 1, -1, 1] ζ

function func_macro(x::Qζ)
    return simplify!(x*x*x*x*x - x*x + x)
end
@benchmark func_macro(data) setup=(data=Qζ(rand(1:10,4),2)) #1.3 microseconds mutable

S, t = PolynomialRing(QQ, "t")
Qzeta, zeta = NumberField(t^4-t^3+t^2-t+1, "ζ")
function func_nemo(x::nf_elem)
    return x*x*x*x*x - x*x + x
end
@benchmark func_nemo(data) setup=(data=rand(1:10) + zeta*rand(1:10)+ zeta^2*rand(1:10) + zeta^3*rand(1:10))
#9 microsecond

function arith(x)
    return div(x[1] * x[2] - x[4], x[3]^5)
end
