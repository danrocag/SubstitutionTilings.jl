module SubstitutionTilings


export substitute, dilate, empirical_frequency, draw, embed_aff, Tiling, SubSystem, id, vertices, in_border, transition_matrix, DGroupElem
export transition_matrix
export collar_in, is_interior, frequency, total_collaring, total_collaring, UnrecognizedCollar
export vertices, @collar_in_from_vertices
export embed_center
export autocorrelation

include("CoreDefs.jl")
using .CoreDefs
include("NumFields.jl")
using .NumFields
#include("Collaring.jl")
#using .Collaring
#include("systems/heisenberg.jl")
#using .Heisenberg
include("systems/penrose.jl")
using .Penrose
include("systems/ammann-beenker.jl")
using .AmmannBeenker
include("systems/pinwheel.jl")
using .Pinwheel
#include("systems/heisenberg-improper.jl")
#using .HeisenbergImproper
include("systems/nilpotent.jl")
using .Nilpotent
include("systems/chair.jl")
using .Chair

end