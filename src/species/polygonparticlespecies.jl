struct PolygonParticleSpecies{F,B<:BindingSite} <: ParticleSpecies{2,F,Pose{2,F,Angle2d{F}}}
    g::NautyDiGraph
    sites::Vector{B}
    corners::Vector{SVector{2,F}}
    rmin::F
    rmax::F
    isregular::Bool
    skin::F
end
function PolygonParticleSpecies(n::Integer, a::F=1.0; colors=1:n, labels=colors) where {F<:Real}
    # TODO: this should have a constructor that takes polygon vertices and automatically creates sites on each polygon side
    r_in = convert(F, 0.5a * cot(π / n))
    r_out = convert(F, 0.5a * csc(π / n))

    tol = sqrt(eps(F)) * r_out
    sites = BindingSite{Pose{2,F,Angle2d{F}},F}[]
    for (i, (c, l)) in enumerate(zip(colors, labels))
        ψ = Angle2d{F}(-F(π) * (1/2 + 2/n * (i-1)))
        x = SVector{2,F}(pol2cart(r_in, rotation_angle(ψ)))
        push!(sites, BindingSite(Pose(x, ψ), c, i:i, tol, tol / r_in))
    end

    g = NautyDiGraph(cycle_digraph(n); vertex_labels=labels)
    corners = [SVector{2,F}(r_out * cos(-π/2 - (2k - 1) * π / n),
                            r_out * sin(-π/2 - (2k - 1) * π / n)) for k in 1:n]
    return PolygonParticleSpecies{F,eltype(sites)}(g, sites, corners, r_in, r_out, true, tol)
end

Base.show(io::Core.IO, ps::PolygonParticleSpecies) = print(io, "$(dimension(ps))d PolygonParticleSpecies with $(nsites(ps)) sites")

Base.copy(pps::PolygonParticleSpecies) = PolygonParticleSpecies(copy(pps.g), copy(pps.sites), copy(pps.corners), pps.rmin, pps.rmax, pps.isregular, pps.skin)

dimension(::PolygonParticleSpecies) = 2

graphrep(p::PolygonParticleSpecies) = p.g
symmetrynumber(p::PolygonParticleSpecies) = 1
nsites(p::PolygonParticleSpecies) = length(p.sites)
bindingsites(p::PolygonParticleSpecies, i::Integer) = p.sites[i]
function setcolors!(p::PolygonParticleSpecies, colors::AbstractVector{<:Integer})
    length(colors) != nsites(p) && throw(ArgumentError("incorrect number of colors"))
    for k in eachindex(p.sites)
        s = p.sites[k]
        p.sites[k] = BindingSite(s.pose, colors[k], s.vertices, s.touching_tolerance, s.alignment_tolerance)
    end
    return
end

function can_skip_overlap_check(p1::PolygonParticleSpecies, p2::PolygonParticleSpecies)
    p1.isregular && p2.isregular || return false
    n = nsites(p1)
    n in (3, 4, 6) || return false
    return nsites(p2) == n
end

isconvex(::PolygonParticleSpecies) = true

function _SAT_overlap(p1::SpeciesAndPose{<:PolygonParticleSpecies},
                      p2::SpeciesAndPose{<:PolygonParticleSpecies})
    spcs1, pose1 = p1
    spcs2, pose2 = p2
    isconvex(spcs1) && isconvex(spcs2) || error("SAT overlap check requires convex polygons")

    # The SAT algorithm checks for intersections of convex polygons by exploiting a simple theorem:
    # For any nonintersecting convex shapes, there exists an axis in space on which the projections
    # of the shapes are nonoverlapping. Therefore, we project all polygon corners on axes given by
    # the surface normals of the polygons. As soon as we identify an axis where projections do not overlap
    # we know that the polygons are not intersecting.
    skin_sum = spcs1.skin + spcs2.skin
    corners1, corners2 = spcs1.corners, spcs2.corners
    for (corners, pose) in ((corners1, pose1), (corners2, pose2))
        n = length(corners)
        for i in 1:n
            edge = pose.psi * (corners[mod1(i + 1, n)] - corners[i])
            normal = SVector(-edge[2], edge[1])
            lo1, hi1 = extrema(dot(normal, pose1 * c) for c in corners1)
            lo2, hi2 = extrema(dot(normal, pose2 * c) for c in corners2)
            (hi2 < lo1 - skin_sum || hi1 < lo2 - skin_sum) && return false
        end
    end
    return true
end

function could_contact(p1::SpeciesAndPose{<:PolygonParticleSpecies}, p2::SpeciesAndPose{<:PolygonParticleSpecies}; kwargs...)
    spcs1, pose1 = p1
    spcs2, pose2 = p2
    return norm(pose1.x - pose2.x) < spcs1.rmax + spcs2.rmax
end
function overlap(p1::SpeciesAndPose{<:PolygonParticleSpecies}, p2::SpeciesAndPose{<:PolygonParticleSpecies}; kwargs...)
    spcs1, pose1 = p1
    spcs2, pose2 = p2
    skin_sum = spcs1.skin + spcs2.skin
    if can_skip_overlap_check(spcs1, spcs2)
        return norm(pose1.x - pose2.x) < (spcs1.rmin + spcs2.rmin) - skin_sum
    else
        return _SAT_overlap(p1, p2)
    end
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