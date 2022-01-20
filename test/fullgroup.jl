@enum PFibTile begin
    α₁=1
    α₂=2
    β=3
end
L = [α₁, α₂, β]
Square = Tuple{PFibTile, PFibTile, PFibTile, PFibTile}
transitions = Dict{PFibTile, Vector{PFibTile}}([α₁ => [α₂,β], α₂ => [α₁,β], β => [α₁]])
H₊ = [[α₁, α₂], [α₁, β], [α₂, β], [β, α₁]]
H₀ = [[ℓ, ℓ] for ℓ in L]
S₊₊ = [ [transitions[h[1]][end], transitions[h[2]][1], h[1], h[2] ] for h in H₊ ]
S₊₀ = [[α₂, β, α₁, α₁], [α₁, β, α₂, α₂]]
S₊ = S₊₊ ∪ S₊₀
S₀ = [ [ℓ, ℓ, κ, κ] for κ in L for ℓ in transitions[κ]]
P = [ [ℓ, κ] for κ in L for ℓ in transitions[κ]]

function complete_up(h)
    return [square for square in S₀ ∪ S₊ if square[[1,2]] == h]
end

Path = Vector{PFibTile}
BBisection = Tuple{Path, Path}
Bisection = Vector{BBisection}

S = [ ([ℓ₁, κ₁], [ℓ₂, κ₂]) for κ₁ in L for ℓ₁ in transitions[κ₁] for κ₂ in L for ℓ₂ in transitions[κ₂]
    if (([κ₁, κ₂] ∈ H₊) || (ℓ₁ < ℓ₂ && κ₁ == κ₂))]
function promote_bisection(B :: Bisection, n)
    m = length(B[1][1])
    for m = m:n-1
        new_B = Bisection()
        for O in B
            h = [O[1][end], O[2][end]]
            for square in complete_up(h)
                push!(new_B, ([O[1];square[3]], [O[2];square[4]]))
            end
        end
        B = new_B
    end
    return B
end

function commutes_from(O :: BBisection, n)
    m = length(O[1])
    for i=n:m-1
        if [O[1][i], O[2][i], O[1][i+1], O[2][i+1]] ∉ S₊ ∪ S₀
            return false
        end
    end
    return true
end

function Base.:∘(B₁ :: Bisection, B₂ :: Bisection)
    if isempty(B₁) || isempty(B₂)
        return Bisection()
    end
    m = length(B₁[1][1])
    n = length(B₂[1][1])
    N = max(m,n)
    B₁ = promote_bisection(B₁, N+2)
    B₂ = promote_bisection(B₂, N+2)
    B = Bisection()
    for O₁ in B₁
        for O₂ in B₂
            if O₁[2] == O₂[1]
                push!(B, (O₁[1], O₂[2]))
            end
        end
    end
    return B
end

part = []
T = filter(!isempty, [[O₁] ∘ [O₂] ∘ [O₃] for O₁ in S for O₂ in S for O₃ in S])
promote_bisection([S[1]] ∘ [S[5]],5)