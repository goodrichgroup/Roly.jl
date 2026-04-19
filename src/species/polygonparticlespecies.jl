mutable struct PolygonParticleSpecies{F,B<:BindingSite} <: ParticleSpecies{2,F,Pose{2,F,Angle2d{F}}}
    g::NautyDiGraph
    sites::Vector{B}
    corners::Matrix{F}
    rmin::F
    rmax::F
end
function PolygonParticleSpecies(n::Integer, a::F=1.0; colors=1:n, labels=colors) where {F<:Real}
    # TODO: this should have a constructor that takes polygon vertices and automatically creates sites on each polygon side
    r_in = convert(F, 0.5a * cot(π / n))
    r_out = convert(F, 0.5a * csc(π / n))

    sites = BindingSite{Pose{2,F,Angle2d{F}}}[]
    for (i, (c, l)) in enumerate(zip(colors, labels))
        ψ = Angle2d{F}(-F(π) * (1/2 + 2/n * (i-1)))
        x = SVector{2,F}(pol2cart(r_in, rotation_angle(ψ))) 
        push!(sites, BindingSite(Pose(x, ψ), c, i:i))
    end

    g = NautyDiGraph(cycle_digraph(n); vertex_labels=labels)
    corners = zeros(F, n, 2)
    return PolygonParticleSpecies{F,eltype(sites)}(g, sites, corners, r_in, r_out)
end

Base.show(io::Core.IO, ps::PolygonParticleSpecies) = print(io, "$(dimension(ps))d PolygonParticleSpecies with $(nsites(ps)) sites")

Base.copy(pps::PolygonParticleSpecies) = PolygonParticleSpecies(copy(pps.g), copy(pps.sites), copy(pps.corners), pps.rmin, pps.rmax)
dimension(::PolygonParticleSpecies) = 2

graphrep(p::PolygonParticleSpecies) = p.g
symmetrynumber(p::PolygonParticleSpecies) = 1
nsites(p::PolygonParticleSpecies) = length(p.sites)
bindingsites(p::PolygonParticleSpecies, i::Integer) = p.sites[i]
function setcolors!(p::PolygonParticleSpecies, colors::AbstractVector{<:Integer})
    length(colors) != nsites(p) && throw(ArgumentError("incorrect number of colors"))
    for k in eachindex(p.sites)
        s = p.sites[k]
        p.sites[k] = BindingSite(s.pose, colors[k], s.vertices)
    end
    return
end

function could_contact(p1::SpeciesAndPose{<:PolygonParticleSpecies}, p2::SpeciesAndPose{<:PolygonParticleSpecies}; kwargs...)
    spcs1, pose1 = p1
    spcs2, pose2 = p2
    return norm(pose1.x - pose2.x) < spcs1.rmax + spcs2.rmax
end
function overlap(p1::SpeciesAndPose{<:PolygonParticleSpecies}, p2::SpeciesAndPose{<:PolygonParticleSpecies}; buffer=0.1, kwargs...)
    spcs1, pose1 = p1
    spcs2, pose2 = p2
    #TODO: this only works if the buildingblocks tile
    return norm(pose1.x - pose2.x) < (spcs1.rmin + spcs2.rmin) * (1 - buffer)
end

const UnitTriangle = PolygonParticleSpecies(3)
const UnitSquare = PolygonParticleSpecies(4)
const UnitHexagon = PolygonParticleSpecies(6)


const INERT_COLOR = colorant"#E7E7E7"
const DEFAULT_COLORS = [
    [colorant"#EAF3FB", colorant"#C9E1F3", colorant"#B2D0EA", colorant"#54A3E4", colorant"#1A78C6", colorant"#00549A"],
    [colorant"#FDCFCE", colorant"#F8B8B7", colorant"#F39F9D", colorant"#ED7E7C", colorant"#E75451", colorant"#BE2E2C"],
    [colorant"#FFF3D9", colorant"#FEECCC", colorant"#FEE3B5", colorant"#FFC563", colorant"#FFB12F", colorant"#FFA000"],
    [colorant"#EDD5F3", colorant"#E2C3EA", colorant"#DBB9E4", colorant"#D99DE8", colorant"#B571C4", colorant"#8E4E9C"],
    [colorant"#D5F5FB", colorant"#BCEEF8", colorant"#AAF5FF", colorant"#86E4F1", colorant"#43CFE2", colorant"#20BCD2"],
    [colorant"#D9F5EC", colorant"#BFEEDD", colorant"#BCF3E1", colorant"#3DD4A4", colorant"#19B684", colorant"#008D60"],
    [colorant"#FFE9D5", colorant"#FFD6BA", colorant"#FFD4B6", colorant"#F7AF7D", colorant"#FF8835", colorant"#ED6304"],
    [colorant"#E0E2FF", colorant"#CED1FF", colorant"#C5C9FE", colorant"#A9AFFF", colorant"#7C85FF", colorant"#4854FD"],
    [colorant"#F9FFE8", colorant"#EDFAC8", colorant"#F4FFC4", colorant"#E0F290", colorant"#B9D63E", colorant"#859F13"],
    [colorant"#FFD5DE", colorant"#FFBCCA", colorant"#FFB2C5", colorant"#FF8EA9", colorant"#FF5C83", colorant"#E32654"],
    [colorant"#D2FDE4", colorant"#B5F8CC", colorant"#95FABA", colorant"#4DD980", colorant"#21AB53", colorant"#008832"]
]

function render!(ax, spcs::PolygonParticleSpecies{F}, pose::Pose=Pose{2,F}();
                 r=nothing, fillcolor=nothing, species_index::Int=1, strokewidth=3, kwargs...) where {F}
    n = nsites(spcs)
    a = 2spcs.rmin * tan(π / n)
    if isnothing(r)
        r = 0.15 * a
    end
    if isnothing(fillcolor)
        palette = DEFAULT_COLORS[mod1(species_index, length(DEFAULT_COLORS))]
        fillcolor = map(enumerate(bindingsites(spcs))) do (i, b)
            b.color == 0 ? INERT_COLOR : palette[mod1(i, length(palette))]
        end
    end
    return _draw_ngon!(ax, pose.x[1], pose.x[2], rotation_angle(pose.psi); n, a, r, color=fillcolor, strokewidth, kwargs...)
end
function render(spcs::PolygonParticleSpecies, pose::Pose=Pose{2}(); r=nothing, kwargs...)
    f = Figure()
    ax = Axis(f[1, 1])
    hidedecorations!(ax)
    return render!(ax, spcs, pose; r, kwargs...)
end

function _ngon_segment(n::Integer, a::Real, r::Real)
    ψ = π / n 
    R = a/2 * cot(ψ)
    x = r * tan(ψ)
    y = a - 2x

    start = Point(r * (tan(ψ) - sin(ψ)),  -r * (1 - cos(ψ)))
    com = Point(a / 2, -R)

    path = BezierPath([
        MoveTo(start - com),
        EllipticalArc(Point(x, -r) - com, r, r, 0, π/2 + ψ, π/2),
        LineTo(Point(x+y, 0) - com),
        EllipticalArc(Point(x+y, -r) - com, r, r, 0, π/2, π/2 - ψ),
        LineTo(Point(x+y, -r) - com),
        LineTo(Point(x, -r) - com),
        ClosePath()])
    return path
end

function _draw_ngon!(ax, x, y, ψ; n, a, r, kwargs...)
    rotation = [ψ + π - 2π/n * i for i in 0:n-1]
    marker = _ngon_segment(n, a, r)
    return scatter!(ax, fill(x, n), fill(y, n); 
                    marker, rotation, markersize=1, markerspace=:data, kwargs...)
end