# Custom particle species

A particle species is a shape together with a set of binding sites on it.
Every species is a subtype of `ParticleSpecies{D,B}`, where `D` is the spatial dimension and `B` is the concrete `BindingSite` type.

Roly ships three general implementations: [`PolygonParticleSpecies`](@ref) for regular polygons, [`PolyhedronParticleSpecies`](@ref) for any convex [`Polyhedron`](@ref) with one site per face, and [`PatchyParticleSpecies`](@ref) for disks and spheres with sites on the surface.
An arbitrary convex solid needs only its corners:

```julia
PolyhedronParticleSpecies(Polyhedron([SVector(x, y, z) for x in (-1.0, 1.0)
                                                       for y in (-2.0, 2.0)
                                                       for z in (-3.0, 3.0)]))
```

Define a new subtype only for a shape none of these covers.

## The interface

A species must define these methods:

| Method | Purpose |
|---|---|
| `graphrep(ps)` | The species' `NautyDiGraph`. See below. |
| `nsites(ps)` | Number of binding sites. |
| `bindingsites(ps, i)` | The `i`th `BindingSite`. |
| `bounding_radius(ps)` | Radius of a sphere at the pose origin enclosing the particle. |
| `overlap(p1::SpeciesAndPose, p2::SpeciesAndPose)` | Whether two particles at given poses overlap. |
| `Base.copy(ps)` | Deep copy. `BindingRules` copies each species to reassign its colors. |

These have defaults you may override:

| Method | Default |
|---|---|
| `isconvex(ps)` | `false`. Set `true` for convex shapes to speed up overlap checks. |
| `could_contact(p1, p2)` | Bounding-sphere pre-check. Override if the shape allows a tighter one. |
| `symmetrynumber(ps)` | Size of the graph's automorphism group, from nauty. |

## The graph representation

Each species carries a directed graph recording how its sites relate, ignoring geometry.
Roly uses it to detect isomorphic polyforms with [nauty](https://pallini.di.uniroma1.it/).
Every automorphism of the graph must be a symmetry of the particle.
The usual choice is a directed cycle with one vertex per site, site `i` on vertex `i`.
[`cycleencoding`](@ref) builds this for any number of sites.
Each `BindingSite` records its vertices in the `vertices` field, so `BindingSite(pose, color, i:i, ...)` puts one site on vertex `i`.

A site may instead span a contiguous range of vertices, which is how 3D species record the twist of a face: [`dartencoding`](@ref) gives each face its own directed cycle.
The graph must then keep the site boundaries fixed, or an automorphism can slide a site across two of them and the symmetry number comes out too large.

## What a binding site records

Besides its pose and color, a [`BindingSite`](@ref) carries three numbers that decide how a partner attaches:

| field | meaning |
|---|---|
| `gauge` | order of the site's own rotational symmetry about its normal. Always 1 in 2D. [`facegauge`](@ref) computes it for a polyhedron face. |
| `stab` | order of the site's stabiliser in the particle's rotation group, from [`sitestabilisers`](@ref). |
| `locking` | whether the site holds its partner in the orientation its frame names (the default) or admits every orientation the shape permits. |

[`nphases`](@ref) reads these to decide how many distinct bonds a pair of sites has, and no graph check catches getting them wrong.

In 2D both `gauge` and `stab` are 1, so the five-argument `BindingSite(pose, color, vertices, touching_tol, alignment_tol)` is right.
In 3D `gauge` is still 1 for a site on a single vertex, but `stab` need not be, since a rotation about a patch's axis can carry the particle onto itself.
Compute it with [`sitestabilisers`](@ref) and pass `BindingSite(pose, color, vertices, tol, tol, gauge, stab)`.
See [Orientation and phases](orientation.md).

## A worked example: rectangle

A non-square rectangle with a site at each edge midpoint, exercising the whole interface without special geometry.

```julia
using Roly
using Roly: BindingSite, ParticleSpecies, SpeciesAndPose
import Roly: graphrep, nsites, bindingsites,
             bounding_radius, isconvex, overlap
using NautyGraphs
using StaticArrays, LinearAlgebra, Rotations
```

### The struct

The struct holds what Roly needs plus your own geometry, here the two side lengths and the four corners in the particle's frame.

```julia
struct Rectangle{F,B<:BindingSite} <: ParticleSpecies{2,B}
    g::NautyDiGraph
    sites::Vector{B}
    corners::Vector{SVector{2,F}}
    width::F
    height::F
    skin::F
end
```

The parameters `{2,B}` fix the dimension to 2D and let Roly infer the pose type from the binding sites.

### The constructor

It builds the graph, puts one site at each edge midpoint facing outward, and computes the corners.

```julia
function Rectangle(width::Real, height::Real; colors=1:4)
    F = float(promote_type(typeof(width), typeof(height)))
    w, h = F(width), F(height)
    tol = sqrt(eps(F)) * max(w, h)

    # Site poses first, so the labelling can be derived from them.
    poses = [
        Pose(SVector{2,F}( w/2, 0),  Angle2d{F}(0)),
        Pose(SVector{2,F}(0,  -h/2), Angle2d{F}(-F(π)/2)),
        Pose(SVector{2,F}(-w/2, 0),  Angle2d{F}(F(π))),
        Pose(SVector{2,F}(0,   h/2), Angle2d{F}(F(π)/2)),
    ]
    # A 2D site has no turn about its in-plane normal, so every gauge is 1.
    labels = siteorbits(poses, ones(Int, 4), collect(colors))

    # Directed 4-cycle, one vertex per binding site; `ranges[i]` is `i:i` here.
    g, ranges = cycleencoding(4; labels)

    # Sites at the edge midpoints, ordered clockwise: right, bottom, left, top.
    sites = [BindingSite(poses[i], colors[i], ranges[i], tol, tol) for i in 1:4]

    corners = [
        SVector{2,F}( w/2,  h/2),
        SVector{2,F}( w/2, -h/2),
        SVector{2,F}(-w/2, -h/2),
        SVector{2,F}(-w/2,  h/2),
    ]

    return Rectangle{F,eltype(sites)}(g, sites, corners, w, h, tol)
end
```

`skin` is a small tolerance for distance comparisons, so sites that should touch count as touching despite floating-point noise.

The labels are derived, not written down.
[`siteorbits`](@ref) puts two sites in one orbit exactly when a rotation carries one onto the other **and** they have the same color.
For this rectangle it returns `[1, 2, 1, 2]` even when all four colors are equal, because opposite edges are interchangeable and adjacent ones are not.
The same call on a square with one color returns `[1, 1, 1, 1]`.

Every built-in species derives its labels this way.
Labels claiming more symmetry than the shape has make the graph merge structures that differ, and nothing downstream notices.
If you write them by hand, call [`check_encoding`](@ref) in your constructor.

### Interface methods

Four of them just read the struct:

```julia
graphrep(ps::Rectangle) = ps.g
nsites(ps::Rectangle) = length(ps.sites)
bindingsites(ps::Rectangle, i::Integer) = ps.sites[i]
bounding_radius(ps::Rectangle) = sqrt(ps.width^2 + ps.height^2) / 2
isconvex(::Rectangle) = true

Base.copy(ps::Rectangle) =
    typeof(ps)(copy(ps.g), copy(ps.sites), copy(ps.corners), ps.width, ps.height, ps.skin)
```

`setcolors!` needs no definition.
Recoloring must re-derive the labelling and the stabilisers too, since the coloring decides both which sites are interchangeable and how many ways a partner attaches.
The generic method does that, and finds the sites in the `sites` field.
A species storing them elsewhere defines its own method.

### Overlap

For convex shapes use the Separating Axis Theorem: two convex bodies are disjoint exactly when their projections onto some axis do not overlap, so a finite set of candidate axes suffices.
Only that set differs between dimensions, so [`sat_overlap`](@ref) takes it from you and does the rest.
In 2D it is the edge normals of both polygons; in 3D both solids' face normals plus the cross products of their edge directions.

```julia
function overlap(p1::SpeciesAndPose{<:Rectangle}, p2::SpeciesAndPose{<:Rectangle}; kwargs...)
    s1, pose1 = p1
    s2, pose2 = p2
    axes = Iterators.flatten((edgenormals(s1.corners, pose1), edgenormals(s2.corners, pose2)))
    return sat_overlap(axes, s1.corners, pose1, s2.corners, pose2, s1.skin + s2.skin)
end
```

Axes need not be normalized or even nonzero, since `sat_overlap` scales each one and skips degenerate ones.

## Using the new species

It now works like any built-in one:

```julia
rect = Rectangle(2.0, 1.0)
bonds = [1 1 1 3;  # right edge binds to left edge
         1 2 1 4]  # bottom edge binds to top edge
sys = BindingRules(bonds, rect)

polys = polygen(sys; maxsize=6)
```
