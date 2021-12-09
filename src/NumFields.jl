module NumFields

export embed_field, simplify!

using StaticArrays
using StructEquality


function embed_field end
function reduce! end
macro simple_number_field(name, polynomial, generator)
    p = eval(polynomial)
    T = eltype(p)
    N = length(p)

    A = zeros(T, N, N)
    for i = 1:N-1
        A[i+1,i] = 1
    end
    A[:, N] = p
    A_gen = SMatrix{N, N, T}(A)
    powers = [(A_gen^i)[:, 1] for i=0:2*N-2]
    println(powers)

    return quote
        mutable struct $name <: Integer
            coeffs :: MVector{$N, $(esc(T))}
            denom :: $(esc(T))
            reduced :: Bool
        end

        Base.promote_rule(::Type{$(esc(name))}, ::Type{$(esc(T))}) = $(esc(name))

        function ($(esc(name)))(n :: Integer)
            n_converted = convert($(esc(T)), n)
            return $(esc(name))([i==1 ? n_converted : 0 for i=1:$N], 1, true)
        end


        function ($(esc(name)))(n :: Rational{Integer})
            x = ($(esc(name)))(numerator(n))
            x.denom = denominator(n)
            return x
        end

        function NumFields.reduce!(x :: $(esc(name)))
            if !x.reduced
                divisor = gcd(x.coeffs..., x.denom)
                x.coeffs = [div(c,divisor) for c in x.coeffs]
                x.denom = div(x.denom, divisor)
                x.reduced = false
            end
            return x
        end


        const powers = $powers
        const $(esc(generator)) = $(esc(name))([i==2 ? 1 : 0 for i=1:$N],1,true)



        function Base.:+(x :: $(esc(name)), y :: $(esc(name)))
            return $(esc(name))(x.coeffs*y.denom + y.coeffs*x.denom, x.denom*y.denom, false)
        end


        function Base.:-(x :: $(esc(name)), y :: $(esc(name)))
            return $(esc(name))(x.coeffs*y.denom - y.coeffs*x.denom, x.denom*y.denom, false)
        end

        function Base.:*(x :: $(esc(name)), y :: $(esc(name)))
            coeffs = zeros(MVector{$N, $(esc(T))})

            @inbounds for i = 0:(2*$N-2)
                factor = 0
                for j = max(0,i-$N+1):min($N-1,i)
                    factor += x.coeffs[j+1] * y.coeffs[i-j+1]
                end
                coeffs += powers[i+1] * factor
            end
            
            return $(esc(name))(coeffs, x.denom*y.denom, false)
        end

        function Base.:*(x :: $(esc(name)), y :: Integer)
            $(esc(name))(x.coeffs*y,x.denom, x.reduced)
        end

        function Base.:*(x :: Integer, y :: $(esc(name)))
            $(esc(name))(y.coeffs*x,y.denom, y.reduced)
        end

        function Base.:(==)(x :: $(esc(name)), y :: $(esc(name)))
            reduce!(x)
            reduce!(y)
            return x.coeffs == y.coeffs && x.denom == y.denom
        end



        function Base.:/(x :: $(esc(name)), y :: Integer)
            return $(esc(name))(x.coeffs,x.denom*y, false)
        end


        function NumFields.embed_field(map_coeff, map_gen, x :: $(esc(name)))
            result = 0 :: typeof(map_gen)
            for i=1:$N
                result += map_coeff(x.coeffs[i])*map_gen^i
            end
            return result
        end
    end
end

end