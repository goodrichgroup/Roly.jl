# The species interface as `docs/src/custom_species.md` documents it. This exists because that
# page's example was wrong for a long time without anything noticing: it labelled all four sites
# of a *rectangle* alike, claiming a 4-fold symmetry the shape does not have, which would have
# merged genuinely different structures for anyone who followed it.
using Roly
using Roly: BindingSite, ParticleSpecies, SpeciesAndPose, site_symmetry, setcolors!, color,
            cycleencoding, sitestabilisers
import Roly: graphrep, nsites, bindingsites, bounding_radius, isconvex
using NautyGraphs, StaticArrays, LinearAlgebra, Rotations

struct Rectangle{F,B<:BindingSite} <: ParticleSpecies{2,B}
    g::NautyDiGraph
    sites::Vector{B}
    corners::Vector{SVector{2,F}}
    width::F
    height::F
    skin::F
end

function Rectangle(width::Real, height::Real)
    F = float(promote_type(typeof(width), typeof(height)))
    w, h = F(width), F(height)
    tol = sqrt(eps(F)) * max(w, h)
    # Opposite edges are interchangeable, adjacent ones are not: a rectangle is 2-fold.
    sites = [
        BindingSite(Pose(SVector{2,F}(w / 2, 0), Angle2d{F}(0)), 1, 1:1, tol, tol),
        BindingSite(Pose(SVector{2,F}(0, -h / 2), Angle2d{F}(-F(π) / 2)), 2, 2:2, tol, tol),
        BindingSite(Pose(SVector{2,F}(-w / 2, 0), Angle2d{F}(F(π))), 1, 3:3, tol, tol),
        BindingSite(Pose(SVector{2,F}(0, h / 2), Angle2d{F}(F(π) / 2)), 2, 4:4, tol, tol),
    ]
    g, _ = cycleencoding(4; labels=[1, 2, 1, 2])
    corners = [SVector{2,F}(w / 2, h / 2), SVector{2,F}(w / 2, -h / 2),
               SVector{2,F}(-w / 2, -h / 2), SVector{2,F}(-w / 2, h / 2)]
    return Rectangle{F,eltype(sites)}(g, sites, corners, w, h, tol)
end

graphrep(ps::Rectangle) = ps.g
nsites(ps::Rectangle) = length(ps.sites)
bindingsites(ps::Rectangle, i::Integer) = ps.sites[i]
isconvex(::Rectangle) = true
bounding_radius(ps::Rectangle) = sqrt(ps.width^2 + ps.height^2) / 2
# Required: `BindingRules` copies each species before shifting its colors.
Base.copy(ps::Rectangle) =
    typeof(ps)(copy(ps.g), copy(ps.sites), copy(ps.corners), ps.width, ps.height, ps.skin)

# Required: separating axes, testing the edge normals of both rectangles.
function Roly.overlap(p1::SpeciesAndPose{<:Rectangle}, p2::SpeciesAndPose{<:Rectangle}; kwargs...)
    s1, pose1 = p1
    s2, pose2 = p2
    skin = s1.skin + s2.skin
    for (s, pose) in ((s1, pose1), (s2, pose2))
        n = length(s.corners)
        for i in 1:n
            e = pose.psi * (s.corners[mod1(i + 1, n)] - s.corners[i])
            axis = SVector(-e[2], e[1])
            lo1, hi1 = extrema(dot(axis, pose1 * c) for c in s1.corners)
            lo2, hi2 = extrema(dot(axis, pose2 * c) for c in s2.corners)
            (hi2 < lo1 + skin || hi1 < lo2 + skin) && return false
        end
    end
    return true
end

@testset "custom species" begin
    r = Rectangle(2.0, 1.0)
    @test nsites(r) == 4
    @test dimension(r) == 2

    # The graph has to claim exactly the symmetry the sites have, or structures get merged or
    # split. A rectangle is 2-fold; labelling all four sites alike would say 4.
    @test symmetrynumber(r) == site_symmetry(r) == 2
    @test Roly._check_encoding(r) === r

    # `setcolors!` needs no definition: the generic method finds the sites in the `sites` field
    # and re-derives the labelling and the stabilisers, which are a function of the coloring.
    setcolors!(r, [5, 6, 5, 6])
    @test [color(bindingsites(r, i)) for i in 1:4] == [5, 6, 5, 6]
    @test symmetrynumber(r) == site_symmetry(r) == 2
    # Colouring every edge alike cannot make a rectangle 4-fold, and the derivation knows it.
    setcolors!(r, fill(7, 4))
    @test symmetrynumber(r) == site_symmetry(r) == 2
    @test length(unique(Roly.labels(graphrep(r)))) == 2
    @test sitestabilisers(r) == fill(1, 4)

    # And it assembles: long edges bond long, short bond short, tiling the plane on a
    # rectangular grid. The counts are *not* the polyomino sequence, because a non-square
    # rectangle is only 2-fold, so shapes a quarter turn apart stay distinct. The first three
    # are checkable by hand: one monomer; two dominoes, one per pair of edges; and four
    # trominoes, the straight one in two orientations and the bent one in two more, the other
    # two bends being 180 degree turns of those.
    setcolors!(r, [1, 2, 1, 2])
    sys = BindingRules([1 1 1 1; 1 2 1 2], r)
    @test [polyenum(sys; maxsize=i)[1] for i in 1:5] == cumsum([1, 2, 4, 13, 35])
end
