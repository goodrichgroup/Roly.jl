mutable struct PolygonParticleSpecies{F,B<:BindingSite} <: ParticleSpecies
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
        push!(sites, BindingSite(Pose(x, ψ), c, l:l))
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
bindingsite(p::PolygonParticleSpecies, i::Integer) = p.sites[i]
function setcolor!(p::PolygonParticleSpecies, c::Integer)
    map!(p.sites) do s
        shift_color(s, c - 1)
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
    [colorant"#B2D0EA", colorant"#54A3E4", colorant"#1A78C6", colorant"#00549A"],
    [colorant"#F39F9D", colorant"#ED7E7C", colorant"#E75451", colorant"#BE2E2C"],
    [colorant"#FEE3B5", colorant"#FFC563", colorant"#FFB12F", colorant"#FFA000"],
    [colorant"#DBB9E4", colorant"#D99DE8", colorant"#B571C4", colorant"#8E4E9C"],
    [colorant"#AAF5FF", colorant"#86E4F1", colorant"#43CFE2", colorant"#20BCD2"],
    [colorant"#BCF3E1", colorant"#3DD4A4", colorant"#19B684", colorant"#008D60"],
    [colorant"#FFD4B6", colorant"#F7AF7D", colorant"#FF8835", colorant"#ED6304"],
    [colorant"#C5C9FE", colorant"#A9AFFF", colorant"#7C85FF", colorant"#4854FD"],
    [colorant"#F4FFC4", colorant"#E0F290", colorant"#B9D63E", colorant"#859F13"],
    [colorant"#FFB2C5", colorant"#FF8EA9", colorant"#FF5C83", colorant"#E32654"],
    [colorant"#95FABA", colorant"#4DD980", colorant"#21AB53", colorant"#008832"]
]

function render!(ax, spcs::PolygonParticleSpecies, pose::Pose=Pose{2}(); 
                 r=nothing, color=nothing, strokewidth=3, kwargs...)
    n = nsites(spcs)
    a = 2spcs.rmin * tan(π / n)
    if isnothing(r)
        r = 0.15 * a
    end
    if isnothing(color)
        sitecolors = [b.color for b in bindingsites(spcs)]
        spcs_color = cld(first(sitecolors), n)
        color = DEFAULT_COLORS[spcs_color][1:n]
        # grey out inert sides (set their color to zero to indicate that)
    end
    return _draw_ngon!(ax, pose.x[1], pose.x[2], rotation_angle(pose.ψ); n, a, r, color, strokewidth, kwargs...)
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