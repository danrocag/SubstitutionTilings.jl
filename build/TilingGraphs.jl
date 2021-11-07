module TilingGraphs

struct TilingGraph{P, V}
    tiles :: Vector{P}
    vertices :: Vector{V}
    shared :: Dict{Int,Int}
    sta_v :: Vector{T}
end

