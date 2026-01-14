# SubstitutionTilings

![splash image](penrose.png)

[![][docs-development-img]][docs-development-url]

**SubstitutionTilings.jl** is a Julia package designed for the exploration, generation, and analysis of substitution tilings. It provides tools for generating large patches of tilings through substitution rules, visualizing them, and computing statistical properties such as patch frequencies.
The package currently includes implementations for several well-known substitution systems:

- **Penrose P2**: The famous aperiodic tiling by Roger Penrose.
- **Ammann-Beenker**: An octagonal aperiodic tiling.
- **Chair**: Also known as the Tiling of the Chair.
- **Pinwheel**: A tiling where tiles appear in infinitely many orientations.
- **Nilpotent**: A substitution system illustrating nilpotent properties.

## Installation

This package is not yet registered in the General registry. You can install it directly from GitHub:

```julia
using Pkg
Pkg.add(url="https://github.com/DanRocag/SubstitutionTilings.jl")
```

## Basic Usage

Here is a simple example of how to generate and draw a Penrose tiling.

```julia
using SubstitutionTilings
using SubstitutionTilings.Penrose
using Luxor

# 1. Define the substitution system
S = penrose()

# 2. Define the initial set of tiles (a single half-kite)
initial_tiles = [hkite()]

# 3. Generate the tiling by applying substitution 5 times
# substitute(System, initial_tiles, iterations)
tiling = substitute(S, initial_tiles, 5)

# 4. Draw the result using Luxor
@png begin
    # Iterate over the generated tiles and draw them
    for tile in tiling
        origin()
        # draw(tile, scale, color, action)
        draw(tile, 20, "black", :stroke)
    end
end 800 800 "penrose_example.png"
```

For more detailed examples, such as the [Fibonacci tiling](https://danrocag.github.io/SubstitutionTilings.jl/dev/fibonacci.html), please refer to the documentation.

[docs-development-img]: https://img.shields.io/badge/docs-development-blue
[docs-development-url]: https://danrocag.github.io/SubstitutionTilings.jl/dev/