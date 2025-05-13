# The Fibonacci tiling

using SubstitutionTilings
using SubstitutionTilings.NumFields
using StructEquality
using Luxor
using Statistics
using Plots
using LinearAlgebra

@NumFields.simple_number_field_concrete Qτ [1,1] τ
Base.promote_rule(::Type{Qτ}, ::Type{<:Integer}) = Qτ

@struct_hash_equal_isequal struct FibElem <: DGroupElem
    a :: Qτ
end

@enum FibTile begin
    A=1
    B=2
end

function Base.:*(x :: FibElem, y :: FibElem)
    return FibElem(x.a + y.a)
end

function SubstitutionTilings.dilate(λ :: Qτ, x :: FibElem)
    return FibElem(λ*x.a)
end

function Base.inv(x :: FibElem)
    return FibElem(-x.a)
end

function SubstitutionTilings.id(::Type{FibElem})
    return FibElem(0)
end

function embed_float(g :: FibElem)
    fτ = (1+sqrt(5))/2
    return embed_field(float,fτ, g.a)
end

function SubstitutionTilings.embed_aff(g :: FibElem)
    return [1, 0, 0, 1, embed_float(g), 0]
end

function LinearAlgebra.norm(g :: FibElem)
    abs(embed_float(g))
end

function SubstitutionTilings.draw(ptile::FibTile, action)
    fτ = (1+sqrt(5))/2
    if ptile == A
        return Luxor.poly([
            Point(-fτ*0.49,-10),
            Point(-fτ*0.49,10),
            Point(fτ*0.49,10),
            Point(fτ*0.49,-10)
        ], close = true, action)
    else
        return Luxor.poly([
            Point(-0.49,-10),
            Point(-0.49,10),
            Point(0.49,10),
            Point(0.49,-10)
        ], close = true, action)
    end
end

fib = SubSystem(Dict(A => [FibElem(Qτ(-1)//2) => A, FibElem(τ//2) => B], B => [FibElem(0) => A]),τ)

fib_tiling = substitute(fib, Dict([FibElem(0) => A]), 3)

function color(tile :: Pair{FibElem,FibTile})
    return (tile[2] == A) ? 1 : 2
end

width = 800
height = 30
sc = 30
@pdf begin
    colors = ["#DD93FC", "#E7977A"]
    patch = [FibElem(-τ//2) => A, FibElem(τ//2) => A]
    tiling = substitute(fib, patch, 5)
    setline(1)

    for tile in tiling
        origin()
        draw(tile, sc, colors[color(tile)], :fill)
    end

    origin()
    scale(sc)
    sethue("black")
    """Luxor.poly([
            Point(-0.05,-10),
            Point(-0.05,10),
            Point(0.05,10),
            Point(0.05,-10)
        ], close = true, :fill)"""
end width height "fibonacci-tiling"

function SubstitutionTilings.is_interior(tiling :: Dict, t :: FibElem)
    return haskey(tiling, t) && (haskey(tiling, FibElem(t.a-τ)) || haskey(tiling, FibElem(t.a-(1+τ)//2))) && (haskey(tiling, FibElem(t.a+τ)) || haskey(tiling, FibElem(t.a+(1+τ)//2)))
end

function SubstitutionTilings.collar_in(tiling :: Dict, t :: FibElem)
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

initial_collar = collar_in(fib_tiling, FibElem(0))
(collars, fib_c) = total_collaring(fib, initial_collar)

collars

width = 800
height = 200
sc = 110
@pdf begin
    colors = ["#DD93FC", "#E7977A"]
    quadrants = Tiler(width, height, 2, 2, margin=0)

    for (pos, n) in quadrants
        println(n)
        first_tile = div(n-1,2) == 0 ? (FibElem(0) => A) : (FibElem(0) => B)
        tiling = substitute(fib, [first_tile], (n-1)%2)
        for tile in tiling
            origin()
            translate(pos)
            scale(sc)
            box(O, 1000, 0.5, :clip)
            sethue(colors[color(tile)])
            transform(embed_aff(tile[1]))
            draw(tile[2], :fill)
            clipreset()
        end
    end
end width height "fibonacci-rule"

function empirical_autocorrelation(T)
    measure = Dict{FibElem, Float64}()
    N = length(T)
    for t in keys(T)
        for s in keys(T)
            g = t*inv(s)
            if haskey(measure, g)
                measure[g] += 1/N
            else
                measure[g] = 1/N
            end
        end
    end
    measure
end
nu_emp_small = empirical_autocorrelation(fib_tiling)

function small_autocorrelation(T, n)
    measure = Dict{FibElem, Float64}()
    measure[FibElem(0)] = 1
    for t in keys(T)
        for s in keys(T)
            g = t*inv(s)
            if !haskey(measure, g)
                measure[g] = frequency(fib, initial_collar,
                    Dict([FibElem(0) => A, g => A]), n)
                    + frequency(fib, initial_collar,
                    Dict([FibElem(0) => A, g => B]), n)
                    + frequency(fib, initial_collar,
                    Dict([FibElem(0) => B, g => A]), n)
                    + frequency(fib, initial_collar,
                    Dict([FibElem(0) => B, g => B]), n)
            end
        end
    end
    measure
end

ix = 0.01:0.001:1
plot(ix, abs.(empirical_diffraction(big_tiling, ix)))

nu = autocorrelation(fib, initial_collar, 15, 1)

function total_freq(g, depth)
    result = 0.
    result += frequency(fib, initial_collar, Dict(FibElem(0) => A, FibElem(g) => A), depth)
    result += frequency(fib, initial_collar, Dict(FibElem(0) => B, FibElem(g) => A), depth)
    result += frequency(fib, initial_collar, Dict(FibElem(0) => A, FibElem(g) => B), depth)
    result += frequency(fib, initial_collar, Dict(FibElem(0) => B, FibElem(g) => B), depth)
end
function check(nu, depth)
    for t in nu
        if t[1] != FibElem(0)
            @assert total_freq(t[1].a, depth) ≈ t[2]
        end
    end
end
check(nu, 18)

using AccurateArithmetic
function finite_diffraction_distr(Nu, x)
    sum([t[1] == FibElem(0) ? t[2]*x : t[2]*sin(norm(t[1])*x)/norm(t[1]) for t in Nu])
end
ix=0:0.01:1
plot(ix, finite_diffraction_distr.(Ref(nu), ix))
small_tiling = substitute(fib, Dict([FibElem(0) => A]), 2)
C = moments[end]/Rs[end]
variances = moments -C.*Rs
vs = variances
#cs = cumsum(variances)
#vs = (cs + circshift(cs,-5))[1:end-5]/5
plot(vs)
plot(log.(abs.(vs))./log.(Rs))