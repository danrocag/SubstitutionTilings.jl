# The Fibonacci tiling

using SubstitutionTilings
using SubstitutionTilings.NumFields
using StructEquality
using Luxor
using Plots
using Statistics

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

function SubstitutionTilings.draw(ptile::FibTile, action)
    fτ = (1+sqrt(5))/2
    if ptile == A
        return Luxor.poly([
            Point(-fτ*0.9,-10),
            Point(-fτ*0.9,10),
            Point(fτ*0.9,10),
            Point(fτ*0.9,-10)
        ], close = true, action)
    else
        return Luxor.poly([
            Point(-0.9,-10),
            Point(-0.9,10),
            Point(0.9,10),
            Point(0.9,-10)
        ], close = true, action)
    end
end

fib = SubSystem(Dict(A => [FibElem(Qτ(-1)//2) => A, FibElem(τ//2) => B], B => [FibElem(0) => A]),τ)

fib_tiling = substitute(fib, Dict([FibElem(0) => A]), 3)

function color(tile :: Pair{FibElem,FibTile})
    return (tile[2] == A) ? 1 : 2
end

width = 800
height = 20
sc = 20
@png begin
    colors = ["#DD93FC", "#E7977A"]
    patch = [FibElem(-τ//2) => A, FibElem(τ//2) => A]
    tiling = substitute(fib, patch, 4)
    setline(1)

    for tile in tiling
        origin()
        draw(tile, sc, colors[color(tile)], :fill)
    end

    origin()
    scale(sc)
    sethue("black")
    Luxor.poly([
            Point(-0.05,-10),
            Point(-0.05,10),
            Point(0.05,10),
            Point(0.05,-10)
        ], close = true, :fill)
end width height

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
height = 40
sc = 20
@png begin
    colors = ["#DD93FC", "#E7977A"]
    first_tile = [FibElem(0) => A]
    tiling = collars[4]
    setline(1)

    for tile in tiling
        origin()
        draw(tile, sc, colors[color(tile)], :fill)
    end
end width height

norm = sqrt(1 + phi^2)
frequency(fib, initial_collar, collars[4], 4)
nu = autocorrelation(fib, initial_collar, 19)
xs = Vector{Float64}()
ys = Vector{Float64}()
for delta in nu
    if abs(embed_float(delta[1])) < 6000
        push!(xs, embed_float(delta[1]))
        push!(ys, delta[2])
    end
end
plot(xs, ys, seriestype = :scatter)

moments = Vector{Float64}()
Rs = 1:6000
for R=Rs
    R = R
    k = R
    val = 0.0
    for delta in nu
        if abs(embed_float(delta[1])) < R
            val += delta[2]
        end
    end
    push!(moments, val)
end

plot(moments)
C = mean((moments./Rs)[5000:6000])
phi = (1 + sqrt(5))/2
covol = (1/phi+1/phi^3)
variances = moments -C.*Rs
plot(variances)
plot(variances./Rs)
plot(log.(abs.(variances))./log.(Rs))