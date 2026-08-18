# Custom particle species

A particle species describes a shape together with a set of binding sites on it. All species are subtypes of `ParticleSpecies{D,B}`, where `D` is the spatial dimension and `B` is the concrete `BindingSite` type used by the species.

Roly ships with three general-purpose implementations: [`PolygonParticleSpecies`](@ref) (regular polygons), [`PolyhedronParticleSpecies`](@ref) (any convex [`Polyhedron`](@ref), with one binding site per face) and [`PatchyParticleSpecies`](@ref) (disks or spheres with binding sites on their surface). Between them they cover most rigid shapes — in particular, an arbitrary convex solid needs only its corners:

```julia
PolyhedronParticleSpecies(Polyhedron([SVector(x, y, z) for x in (-1.0, 1.0)
                                                       for y in (-2.0, 2.0)
                                                       for z in (-3.0, 3.0)]))
```

To model a shape none of these covers, define a new subtype of `ParticleSpecies` and implement the interface described below.

## The interface

A concrete `ParticleSpecies` must define the following methods:

| Method | Purpose |
|---|---|
| `graphrep(ps)` | Return the graph representation of the particle species (a `NautyDiGraph`). See the note below. |
| `nsites(ps)` | Number of binding sites. |
| `bindingsites(ps, i)` | Return the `i`th `BindingSite`. |
| `bounding_radius(ps)` | Radius of a sphere centered at the pose origin that fully encloses the particle. |
| `overlap(p1::SpeciesAndPose, p2::SpeciesAndPose)` | Return `true` if two particles at given poses overlap. |
| `Base.copy(ps)` | Deep copy. `BindingRules` copies each species during construction to reassign its colors. |

The following methods have generic fallbacks that you can override:

| Method | Default behavior |
|---|---|
| `isconvex(ps)` | Returns `false`. Set to `true` for convex shapes as a hint to overlap checks. |
| `could_contact(p1, p2)` | Cheap bounding-sphere pre-check. Override if the shape allows a tighter test. |
| `symmetrynumber(ps)` | Size of the graph automorphism group, computed via nauty. |

## The graph representation

Each particle species carries a directed graph (`graphrep(ps)`) that captures the *combinatorial* structure of its binding sites, independent of geometry. Roly uses it to detect polyform isomorphism via [nauty](https://pallini.di.uniroma1.it/) and to canonically order the vertices of an assembled polyform.

**Requirements.** Every automorphism of the graph must correspond to a symmetry of the particle, and the site partition must be respected. The most common choice is a directed cycle of `n` vertices, one vertex per binding site, with site `i` occupying vertex `i` — [`cycleencoding`](@ref) builds exactly this, for any `n >= 1`. Each `BindingSite` records which vertex range it covers via its `vertices` field, so `BindingSite(pose, color, i:i, ...)` places one site on vertex `i`.

**Sites spanning several vertices.** A site may occupy a contiguous *range* of vertices rather than a single one, which is how 3D species encode the twist of a face: [`dartencoding`](@ref) gives each face of a polyhedron a directed cycle of its own. When a site spans several vertices, make sure the graph structure still pins the site boundaries down — otherwise an automorphism can slide a site onto a straddling set of vertices and the symmetry number comes out too large. In `dartencoding` the bond-pairing edges are bidirectional while the face-cycle arcs are not, so the two classes cannot mix and the faces are preserved as blocks.

## What a binding site records

Beyond its pose and color, a [`BindingSite`](@ref) carries three numbers that decide how a partner may attach:

| field | meaning |
|---|---|
| `gauge` | order of the site's *own* rotational symmetry about its outward normal. 1 in 2D, where there is no such rotation; [`facegauge`](@ref) computes it for a polyhedron face. |
| `stab` | order of the site's stabiliser in the *particle's* rotation group, from [`sitestabilisers`](@ref). |
| `locking` | whether the site holds its partner in the orientation its frame names (the default) or admits every orientation its shape permits. |

They are what [`nphases`](@ref) reads to decide how many distinct bonds a pair of sites has, and no graph check catches getting them wrong.

In **2D** both are always 1: a site has no turn about its in-plane normal, and no rotation about the particle's origin fixes a site away from it. The five-argument `BindingSite(pose, color, vertices, touching_tol, alignment_tol)` therefore says exactly the right thing, which is why the example below uses it. In **3D** `gauge` is still 1 for a site occupying a single graph vertex — one vertex has no room to record a turn — but `stab` need not be, since a rotation about a patch's own axis can carry the whole particle onto itself. Compute it with [`sitestabilisers`](@ref) and pass it: `BindingSite(pose, color, vertices, tol, tol, gauge, stab)`. See [Orientation and phases](orientation.md) for what follows from these.

## A worked example: rectangle

Let's define a species for a non-square rectangle with binding sites at the midpoints of its four edges. This is a good template because it exercises every part of the interface without requiring specialized geometry.

```julia
using Roly
using Roly: BindingSite, ParticleSpecies, SpeciesAndPose
import Roly: graphrep, nsites, bindingsites,
             bounding_radius, isconvex, overlap
using NautyGraphs
using StaticArrays, LinearAlgebra, Rotations
```

### The struct

The struct stores everything Roly needs and any geometric parameters specific to the shape. Here we store the two side lengths and the four corners in the particle's local frame.

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

The type parameters `{2,B}` fix the spatial dimension to 2D and let Roly infer the pose type from the binding sites.

### The constructor

The constructor builds the graph, places one binding site at the midpoint of each edge (facing outward), and computes the corners. We use `Angle2d` for 2D rotations and `Pose` to describe each site's position and orientation.

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

The `skin` field is a small numerical tolerance used when comparing distances, so sites that should touch are recognized as touching despite floating-point noise.

Note that the labels are *derived*, not written down. [`siteorbits`](@ref) puts two sites in one orbit exactly when a rotation carries one onto the other **and** they are the same color, so it reads the symmetry off the geometry you already supplied. For this rectangle it returns `[1, 2, 1, 2]` even when all four colors are equal: opposite edges are interchangeable, adjacent ones are not, because a rectangle is 2-fold and not 4-fold. The same call on a square with one color returns `[1, 1, 1, 1]`.

This is worth doing rather than writing the labels by hand, and it is what all of Roly's own species do. A labelling that claims more symmetry than the shape has makes the graph merge structures that are genuinely different — silently, since nothing downstream re-examines it. If you do write labels yourself, call [`Roly._check_encoding`](@ref) once to confirm `symmetrynumber` and `site_symmetry` agree.

### Interface methods

The four "structural" methods just read from the struct:

```julia
graphrep(ps::Rectangle) = ps.g
nsites(ps::Rectangle) = length(ps.sites)
bindingsites(ps::Rectangle, i::Integer) = ps.sites[i]
bounding_radius(ps::Rectangle) = sqrt(ps.width^2 + ps.height^2) / 2
isconvex(::Rectangle) = true

Base.copy(ps::Rectangle) =
    typeof(ps)(copy(ps.g), copy(ps.sites), copy(ps.corners), ps.width, ps.height, ps.skin)
```

`setcolors!` needs no definition at all. Recoloring has to re-derive the labelling and the site stabilisers along with the colors, since a coloring is what decides both which sites are interchangeable and how many ways a partner can attach; the generic method does all of it, and finds the sites in the `sites` field. That field name is the only requirement — a species keeping them elsewhere defines its own method.

### Overlap

The overlap check decides whether two particles at given poses intersect. For convex shapes the standard approach is the *Separating Axis Theorem* (SAT): two convex bodies are disjoint if and only if their projections onto some axis do not overlap, so it suffices to test a finite set of candidate axes. Which set that is, is all that differs between dimensions — in 2D the edge normals of both polygons, in 3D both solids' face normals plus the cross products of their edge directions — so [`sat_overlap`](@ref) takes the candidates from you and does the rest:

```julia
function overlap(p1::SpeciesAndPose{<:Rectangle}, p2::SpeciesAndPose{<:Rectangle}; kwargs...)
    s1, pose1 = p1
    s2, pose2 = p2
    axes = Iterators.flatten((edgenormals(s1.corners, pose1), edgenormals(s2.corners, pose2)))
    return sat_overlap(axes, s1.corners, pose1, s2.corners, pose2, s1.skin + s2.skin)
end
```

Axes need not be normalized or even nonzero — `sat_overlap` scales each one and skips the degenerate ones. That is not cosmetic: `skin` is a length, so comparing it against projections along an unnormalized axis would scale the tolerance by that axis's magnitude, which is exactly the bug this page used to demonstrate.

## Using the new species

Once the interface is implemented, the new species can be used just like any built-in one:

```julia
rect = Rectangle(2.0, 1.0)
bonds = [1 1 1 3;  # right edge binds to left edge
         1 2 1 4]  # bottom edge binds to top edge
sys = BindingRules(bonds, rect)

polys = polygen(sys; maxsize=6)
```
