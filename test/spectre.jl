using SubstitutionTilings.NumFields
using Luxor
using StructEquality


@NumFields.simple_number_field_concrete Qθ [3,0] θ
promote_rule(::Type{Qθ}, ::Type{Rational{Int64}}) = Qθ
K = Qθ
sq3 = θ
rotmatrix = K[sq3//2 1//2; -1//2 sq3//2]

function embed_float(x::Qθ)
    embed_field(Float64, sqrt(3), x)
end
function embed_float(v::Vector{Qθ})
    return embed_float.(v)
end
function Luxor.translate(v::Vector{Qθ})
    translate(Point(Tuple(embed_float(v))))
end
function Luxor.Point(v :: Vector{Qθ})
    Point(embed_float(v)...)
end

@def_structequal struct G
    t :: Vector{Qθ}
    k :: Int
end

function Base.:*(g :: G, v :: Vector{Qθ})
    g.t + rotmatrix^mod(g.k,12)*v
end
function Base.:*(g, h :: G)
    G(g*h.t, mod(g.k+h.k,12))
end
function Base.inv(g :: G)
    G(-rotmatrix^mod(-g.k,12)*g.t, mod(-g.k,12))
end
function rot(k::Int)
    G(K[0;0],k)
end

function matrix(g :: G)
    [rotmatrix^mod(g.k,12) g.t]
end
function Luxor.Point(g :: G)
    Point(g.t)
end
function Luxor.transform(g :: G)
    m = embed_float.(matrix(g))
    m_arr = [m[1:2,1];m[1:2,2];m[1:2,3]]
    println(m_arr)
    transform(m_arr)
end

@def_structequal struct Figure
    angles :: Vector{Int}
    marks :: Vector{Int}
end


spectre_angles = [2, -3, 2, 0, 2, 3, -2, 3, -2, 3, 2, -3, 2, 3]
spectre_indices = [2,4,10,14]
spectre = Figure(spectre_angles, spectre_indices)
mystic_angles = [spectre_angles[1:5];-3; spectre_angles[5:end-1]; 0; spectre_angles[11:end]]
mystic_indices = [8,14]
mystic = Figure(mystic_angles, mystic_indices)

function displacements(fig :: Figure)
    displacements = G[G(K[0;0], 0)]
    current_displacement = G(K[1;0], 0)
    for a in fig.angles
        push!(displacements, current_displacement)
        current_displacement = current_displacement*rot(a)*G(Qθ[1;0],0)
    end
    @assert current_displacement == G(K[1;0], 0)
    return [displacements[end];displacements[2:end]]
end

function match(f :: Figure, i :: Int, g :: Figure, j::Int)
    df = dirs(f)
    pf = points(df)
    dg = dirs(g)
    pg = points(dg)
    ix = f.marks[i]
    jx = g.marks[j]
    
end

@draw begin
    f1 = mystic
    p1 = displacements(f1)
    checks1 = p1[f1.marks]
    #f2 = spectre
    #p2 = points(f2)
    #checks2 = p2[f2.marks]
    scale(40)
    poly(Point.(p1), action=:stroke)
    for check in checks1
        circle(Point(check), 0.1, action=:fillstroke)
        arrow(Point(check), Point(check*G(K[1;0],0));arrowheadlength=0.3)
    end
    #g = match(f1,1,f2,1)
    #transform(g)
    #poly(Point.(p2), action=:stroke)
    #for check in checks2
    #    circle(Point(check), 0.1, action=:fillstroke)
    #end

end 800 800
