# The species interface as `docs/src/custom_species.md` documents it.
using Roly
using Roly: BindingSite, ParticleSpecies, SpeciesAndPose, site_symmetry, setcolors!, color,
            cycleencoding, sitestabilisers, sat_overlap, edgenormals
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

function Rectangle(width::Real, height::Real; colors=1:4)
    F = float(promote_type(typeof(width), typeof(height)))
    w, h = F(width), F(height)
    tol = sqrt(eps(F)) * max(w, h)
    poses = [
        Pose(SVector{2,F}(w / 2, 0), Angle2d{F}(0)),
        Pose(SVector{2,F}(0, -h / 2), Angle2d{F}(-F(π) / 2)),
        Pose(SVector{2,F}(-w / 2, 0), Angle2d{F}(F(π))),
        Pose(SVector{2,F}(0, h / 2), Angle2d{F}(F(π) / 2)),
    ]

    labels = siteorbits(poses, ones(Int, 4), collect(colors))
    g, ranges = cycleencoding(4; labels)
    sites = [BindingSite(poses[i], colors[i], ranges[i], tol, tol) for i in 1:4]
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

function Roly.overlap(p1::SpeciesAndPose{<:Rectangle}, p2::SpeciesAndPose{<:Rectangle}; kwargs...)
    s1, pose1 = p1
    s2, pose2 = p2
    axes = Iterators.flatten((edgenormals(s1.corners, pose1), edgenormals(s2.corners, pose2)))
    return sat_overlap(axes, s1.corners, pose1, s2.corners, pose2, s1.skin + s2.skin)
end

@testset "custom species" begin
    r = Rectangle(2.0, 1.0; colors=[1, 2, 1, 2])
    @test nsites(r) == 4
    @test dimension(r) == 2

    # check automatic symmetry encoding
    @test symmetrynumber(r) == site_symmetry(r) == 2
    @test Roly._check_encoding(r) === r
    for cols in ([1, 1, 1, 1], [1, 2, 1, 2], [9, 9, 9, 9])
        @test symmetrynumber(Rectangle(2.0, 1.0; colors=cols)) == 2
    end
    # All four edges alike on a square really are interchangeable
    @test siteorbits([Pose(SVector(cos(t), sin(t)), Angle2d(t)) for t in (0, π/2, π, 3π/2)],
                     ones(Int, 4), fill(1, 4)) == [1, 1, 1, 1]
    # Distinct colors break the symmetry
    @test symmetrynumber(Rectangle(2.0, 1.0; colors=1:4)) == 1

    # `setcolors!` needs no definition: the generic method finds the sites in the `sites` field
    # and re-derives the labelling and the stabilisers
    setcolors!(r, [5, 6, 5, 6])
    @test [color(bindingsites(r, i)) for i in 1:4] == [5, 6, 5, 6]
    @test symmetrynumber(r) == site_symmetry(r) == 2
    # Colouring every edge alike cannot make a rectangle 4-fold symmetric
    setcolors!(r, fill(7, 4))
    @test symmetrynumber(r) == site_symmetry(r) == 2
    @test length(unique(Roly.labels(graphrep(r)))) == 2
    @test sitestabilisers(r) == fill(1, 4)

    # Check enumeration
    setcolors!(r, [1, 2, 1, 2])
    sys = BindingRules([1 1 1 1; 1 2 1 2], r)
    @test [polyenum(sys; maxsize=i)[1] for i in 1:5] == cumsum([1, 2, 4, 13, 35])
end
