using SubstitutionTilings.NumFields

@NumFields.simple_number_field_concrete Qρ  [-36, 0, 8, 0] ρ

@NumFields.simple_number_field_lazy Qζ [-1, 1, -1, 1] ζ
Base.promote_rule(::Type{Qζ}, ::Type{<:Integer}) = Qζ

function func_macro(x::Qζ)
    return reduce!(x*x*x*x*x - x*x + x)
end
@benchmark func_macro(data) setup=(data=Qζ(rand(1:1000,4)))
# median 749 ns microseconds mutable

@NumFields.simple_number_field_concrete Qζ2 [-1, 1, -1, 1] ζ2
function Base.show(io::IO, x::Qζ2)
    result = "("
    some_shown = false
    for (i, coeff) in enumerate(x.coeffs)
        if (coeff != 0)
            if some_shown
                result = result*" + "
            else
                some_shown = true
            end
            result = result*string(coeff)
            result = result*"ζ2^"
            result = result*string(i-1)
        end
    end
    result = result*")/"
    result = result*string(x.denom)
    print(io, result)
end


Base.promote_rule(::Type{Qζ2}, ::Type{<:Integer}) = Qζ2
function func_concrete(x)
    return x*x*x*x*x - x*x + x
end
@benchmark func_concrete(data) setup=(data=Qζ2(rand(1:1000,4)))
# median 1.7 micros

function arith(x)
    return div(x[1] * x[2] - x[4], x[3]^5)
end
