module SubstitutionTilings
__precompile__(false)


export substitute, dilate, empirical_frequency, draw, embed_aff, Tiling, SetTiling, id, vertices, in_border

include("CoreDefs.jl")
using .CoreDefs
#include("systems/heisenberg.jl")
#using .Heisenberg
include("systems/penrose.jl")
using .Penrose
include("systems/pinwheel.jl")
using .Pinwheel
#include("systems/heisenberg-improper.jl")
#using .HeisenbergImproper
include("systems/nilpotent.jl")
using .Nilpotent
include("systems/chair.jl")
using .Chair

end