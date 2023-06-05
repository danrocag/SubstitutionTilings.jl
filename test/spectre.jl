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
spectre_indices = [2,6,12,14,11,5]
spectre = Figure(spectre_angles, spectre_indices)
mystic_angles = [spectre_angles[1:5];-3; spectre_angles[5:end-1]; 0; spectre_angles[11:end]]
mystic_indices = [8,14]
mystic = Figure(mystic_angles, mystic_indices)

function displacements(angles :: Vector{Int})
    displacements = G[G(K[0;0], 0)]
    current_displacement = G(K[1;0], 0)
    for a in angles
        push!(displacements, current_displacement)
        current_displacement = current_displacement*rot(a)*G(Qθ[1;0],0)
    end
    #@assert current_displacement == G(K[1;0], 0)
    return [displacements[end];displacements[2:end]]
end
function displacements(fig :: Figure)
    displacements(fig.angles)
end

function points(displacements :: Vector{G})
    result = Vector{Qθ}[]
    for g in displacements
        push!(result, g.t)
    end
    result
end
function points(fig :: Figure)
    points(displacements(fig))
end
function directions(displacements :: Vector{G})
    result = Int[]
    for g in displacements
        push!(result, g.k)
    end
    result
end
function directions(figure :: Figure)
    directions(displacements(figure))
end


function match(f :: Figure, i :: Int, g :: Figure, j::Int)
    pf = displacements(f)
    pg = displacements(g)
    ix = f.marks[i]
    jx = g.marks[j]
    arrow_f = G(pf[ix].t, pf[ix+1].k)
    arrow_g = G(pg[jx].t, mod(6+pg[jx].k,12))
    arrow_f*inv(arrow_g)
end

function build_supertile(init_figure :: Figure, next :: Vector{Tuple{Int, Int, Figure}})
    #@assert init_figure == next[end][3]
    current_figure = init_figure
    result = Tuple{G, Figure}[]
    g = G(K[0;0],0)
    for (i,j,next_figure, _) in next
        g = g*match(current_figure, i, next_figure, j)
        push!(result, (g, next_figure))
        current_figure = next_figure
    end
    #@assert g == G(K[0;0],0)
    return result
end

function supertile_instructions(spectre, mystic)
    spectre_supertile_instructions = [
        (2,3,spectre, 4),
        (2,4,spectre, nothing),
        (1,3,spectre, nothing),
        (2,4,spectre, 1),
        (2,4,spectre, nothing),
        (1,3,spectre, 4),
        (2,4,spectre, 1),
        (2,1,mystic, nothing)
    ]
    mystic_supertile_instructions = [
        (2,3,spectre, 4),
        (2,4,spectre, nothing),
        (6,5,spectre, 1),
        (2,4,spectre, nothing),
        (1,3,spectre, 4),
        (2,4,spectre, 1),
        (2,1,mystic, nothing)
    ]
    return (spectre_supertile_instructions, mystic_supertile_instructions)
end

function spectre_supertile_instructions(spectre, mystic)
    supertile_instructions(spectre, mystic)[1]
end
function mystic_supertile_instructions(spectre, mystic)
    supertile_instructions(spectre, mystic)[2]
end

@draw begin
    scale(20)
    #poly(Point.(displacements(mystic)), action=:stroke)
    for (g,figure) in build_supertile(mystic, spectre_supertile_instructions)
        ps = Ref(g).*displacements(figure)
        poly(Point.(ps), action=:stroke)
        setopacity(0.5)
        for mark in figure.marks
            circle(Point(ps[mark]), 0.2, action=:fill)
        end
        setopacity(1)
    end
end 800 800
    
function outline(instructions)
    angles = Int[]
    marks = []
    last_connecting_mark = 1
    (i,j, figure, new_mark) = instructions[1]
    println(figure)
    jx = figure.marks[j]
    current_angles = circshift(figure.angles, 1-jx)
    last_connecting_angle = current_angles[end]
    last_figure = figure
    if !isnothing(new_mark)
        println("Marks:")
        println(last_connecting_mark)
        println(mod(figure.marks[new_mark]-jx, length(figure.angles)+1))
        push!(marks, last_connecting_mark+mod(figure.marks[new_mark]-jx, length(figure.angles)))
    end
    for (i,j,figure, new_mark) in instructions[2:end]
        ix = last_figure.marks[i]
        slice = 1:mod(ix-jx,length(current_angles))-1
        append!(angles, current_angles[slice])
        last_connecting_mark += length(slice)+1
        connecting_angle_i = current_angles[mod(ix-jx,length(current_angles))]
        jx = figure.marks[j]
        current_angles = circshift(figure.angles, 1-jx)
        connecting_angle_j = current_angles[end]
        push!(angles, connecting_angle_i+connecting_angle_j-6)
        last_figure = figure
        if !isnothing(new_mark)
            println("Marks:")
            println(last_connecting_mark)
            println(mod(figure.marks[new_mark]-jx, length(figure.angles)+1))
            push!(marks, last_connecting_mark+mod(figure.marks[new_mark]-jx, length(figure.angles)))
        end
    end
    ix = last_figure.marks[i]
    slice = 1:mod(ix-jx,length(current_angles))-1
    append!(angles, current_angles[slice])
    connecting_angle_i = current_angles[mod(ix-jx,length(current_angles))]
    push!(angles, connecting_angle_i+last_connecting_angle-6)
    return Figure(angles, marks)
end

function iterate(figures)
    (spectre, mystic) = figures
    superspectre_am = outline(spectre)
end

@draw begin
    instrs = supertile_instructions(spectre, mystic)
    (spectre2, mystic2) = outline.(instrs)
    figure = outline(spectre_supertile_instructions(spectre, mystic))
    println(figure)

    scale(20)
    poly(Point.(displacements(figure)), action=:stroke)
    ds = displacements(figure)
    for mark in figure.marks
        circle(Point(ds[mark]), 0.15, action=:fill)
    end
    circle(Point(0,0),0.2,action=:fill)
end 800 800
