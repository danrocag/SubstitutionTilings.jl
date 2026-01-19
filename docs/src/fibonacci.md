# Examples

## The Fibonacci substitution rule

In this page, we explain how to use this library to define the Fibonacci substitution system and associated tilings.


```@example 1
using SubstitutionTilings
```

In order to define a substitution system, first we need to define the field where the coordinates and scaling factor lie.
In order to be able to do exact frequency computation,
we need to work in a group with exact equality.
We use our own implementation of ``\mathbb{Q}(\tau)``,
where ``\tau = \frac{1 + \sqrt{5}}{2}``.
The libraries `Nemo.jl` and `Hecke.jl` provide implementations of number fields.

```@example 1
using SubstitutionTilings.NumFields
using StructEquality

@NumFields.simple_number_field_concrete Qτ [1,1] τ
Base.promote_rule(::Type{Qτ}, ::Type{<:Integer}) = Qτ
```

Now we define the group of coordinates of the tiles.
It must implement the operations `*`, `inv` and `dilate`.
In this case it is just a wrapper around `Qτ`.

```@example 1
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
```

Now we can define the substitution system as follows.

Then the substitution is defined by
```@example 1
fib = SubSystem(Dict(A => [FibElem(Qτ(-1)//2) => A, FibElem(τ//2) => B], B => [FibElem(0) => A]),τ)
```

```@example 1
fib_tiling = substitute(fib, Dict([FibElem(0) => A]), 3)
```

In order to be able to draw a tiling we need to define the functions `embed_aff`,
which tells us the affine coordinate matrix of a `FibElem`,
and `draw`, which returns the shape of a prototile.

```@example 1
function SubstitutionTilings.embed_aff(g :: FibElem)
    fτ = (1+sqrt(5))/2
    return [1, 0, 0, 1, embed_field(float,fτ, g.a), 0]
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

function color(tile :: Pair{FibElem,FibTile})
    return (tile[2] == A) ? 1 : 2
end

width = 800
height = 20
sc = 5
@draw begin
    colors = ["#DD93FC", "#E7977A"]
    first_tile = [FibElem(0) => A]
    tiling = substitute(fib, Dict([FibElem(0) => A]), 16)
    setline(1)

    for tile in tiling
        origin()
        draw(tile, sc, colors[color(tile)], :fill)
    end
end width height
```

In order to calculate frequencies, we need to implement `is_interior` and `collar_in`.

```@example 1
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
```

`total_collaring` returns a list of all possible collarings of every prototile.
For the Fibonacci sytem there are four of them.

The total collaring `SubstitutionTilings.jl` computes has 4 collars:
```@example 1
initial_collar = collar_in(fib_tiling, FibElem(0))
(collars, fib_c) = total_collaring(fib, initial_collar)

collars
```


Now we can compute the frequency of the following patch.

```@example 1
width = 800
height = 40
sc = 20
@draw begin
    colors = ["#DD93FC", "#E7977A"]
    first_tile = [FibElem(0) => A]
    tiling = collars[4]
    setline(1)

    for tile in tiling
        origin()
        draw(tile, sc, colors[color(tile)], :fill)
    end
end width height
```

```@example 1
frequency(fib, initial_collar, collars[4], 4)
```