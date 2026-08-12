"""
    PolygonParticleSpecies{F}

A 2D convex regular polygon with `n` edges and one binding site per edge.
"""
struct PolygonParticleSpecies{F,B<:BindingSite} <: ParticleSpecies{2,B}
    g::NautyDiGraph
    sites::Vector{B}
    corners::Vector{SVector{2,F}}
    rmin::F
    rmax::F
    isregular::Bool
    skin::F
end

"""
    PolygonParticleSpecies(n, a=1.0; colors=1:n)

Construct a regular `n`-gon with edge length `a`. Each edge carries one binding site.

`colors` assigns interaction colors to the binding sites, which is what the interaction matrix
uses to decide which sites bond. It is also all that needs saying: the particle's symmetry
follows from it, since two edges are interchangeable exactly when a rotation carries one onto
the other and they are the same color (see [`siteorbits`](@ref)).

So the default `colors=1:n` gives every edge its own identity and a symmetry number of 1, while
`colors=fill(1, n)` makes all edges the same sticky stuff and gives the full `n`.
"""
function PolygonParticleSpecies(n::Integer, a::F=1.0; colors=1:n) where {F<:Real}
    r_in = convert(F, 0.5a * cot(π / n))
    r_out = convert(F, 0.5a * csc(π / n))

    tol = sqrt(eps(F)) * r_out
    poses = map(1:n) do i
        ψ = Angle2d{F}(-F(π) * (1 / 2 + 2 / n * (i - 1)))
        Pose(SVector{2,F}(pol2cart(r_in, rotation_angle(ψ))), ψ)
    end
    # A 2D site has no turn about its in-plane normal, so its gauge is 1 throughout, and with
    # it the stabiliser: the only rotation fixing a site is the identity.
    gauges = ones(Int, n)
    labels = siteorbits(poses, gauges, collect(colors))
    stabs = sitestabilisers(poses, gauges, labels)

    g, ranges = cycleencoding(n; labels)
    sites = BindingSite{Pose{2,F,Angle2d{F}},F}[]
    for (i, c) in enumerate(colors)
        push!(sites, BindingSite(poses[i], c, ranges[i], tol, tol / r_in, 1, stabs[i]))
    end

    corners = [
        SVector{2,F}(r_out * cos(-π / 2 - (2k - 1) * π / n), r_out * sin(-π / 2 - (2k - 1) * π / n)) for k in 1:n
    ]
    return _check_encoding(
        PolygonParticleSpecies{F,eltype(sites)}(g, sites, corners, r_in, r_out, true, tol)
    )
end

function Base.show(io::Core.IO, ps::PolygonParticleSpecies)
    return print(io, "$(dimension(ps))d PolygonParticleSpecies with $(nsites(ps)) sites")
end

function Base.copy(pps::PolygonParticleSpecies)
    return PolygonParticleSpecies(
        copy(pps.g), copy(pps.sites), copy(pps.corners), pps.rmin, pps.rmax, pps.isregular, pps.skin
    )
end


graphrep(p::PolygonParticleSpecies) = p.g
nsites(p::PolygonParticleSpecies) = length(p.sites)
bindingsites(p::PolygonParticleSpecies, i::Integer) = p.sites[i]
function can_skip_overlap_check(p1::SpeciesAndPose{<:PolygonParticleSpecies}, p2::SpeciesAndPose{<:PolygonParticleSpecies})
    spcs1, pose1 = p1
    spcs2, pose2 = p2

    norm(pose1.x - pose2.x) >= spcs1.rmax + spcs2.rmax && return true
    spcs1.isregular && spcs2.isregular || return false

    n = nsites(spcs1)
    n in (3, 4, 6) && nsites(spcs2) == n || return false

    sym_angle = 2π / n
    θ_rel = mod(rotation_angle(pose2.psi) - rotation_angle(pose1.psi), sym_angle)
    tol = (spcs1.skin + spcs2.skin) / spcs1.rmin
    face_to_face = mod(π, sym_angle)
    return θ_rel < tol || θ_rel > sym_angle - tol || abs(θ_rel - face_to_face) < tol
end

isconvex(::PolygonParticleSpecies) = true

function overlap(p1::SpeciesAndPose{<:PolygonParticleSpecies}, p2::SpeciesAndPose{<:PolygonParticleSpecies}; kwargs...)
    spcs1, pose1 = p1
    spcs2, pose2 = p2
    # Where the cheap test applies, the inradius decides; otherwise fall back to separating
    # axes, for which a 2D polygon's edge normals are a sufficient candidate set.
    can_skip_overlap_check(p1, p2) ||
        return sat_overlap(Iterators.flatten((edgenormals(spcs1.corners, pose1),
                                              edgenormals(spcs2.corners, pose2))),
                           spcs1.corners, pose1, spcs2.corners, pose2, spcs1.skin + spcs2.skin)
    return norm(pose1.x - pose2.x) < (spcs1.rmin + spcs2.rmin) - (spcs1.skin + spcs2.skin)
end

bounding_radius(ps::PolygonParticleSpecies) = ps.rmax

"""
    UnitTriangle

An equilateral triangle with unit-length edges and one binding site per edge.
"""
const UnitTriangle = PolygonParticleSpecies(3)

"""
    UnitSquare

A square with unit-length edges and one binding site per edge.
"""
const UnitSquare = PolygonParticleSpecies(4)

"""
    UnitHexagon

A regular hexagon with unit-length edges and one binding site per edge.
"""
const UnitHexagon = PolygonParticleSpecies(6)
