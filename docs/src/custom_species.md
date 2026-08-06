# Custom particle species

A particle species describes a shape together with a set of binding sites on it. All species are subtypes of `ParticleSpecies{D,B}`, where `D` is the spatial dimension and `B` is the concrete `BindingSite` type used by the species.

Roly ships with two general-purpose implementations, [`PolygonParticleSpecies`](@ref) (regular polygons) and [`PatchyParticleSpecies`](@ref) (disks or spheres with binding sites on their surface). To model a shape neither of these covers, you define a new subtype of `ParticleSpecies` and implement the interface described below.

## The interface

A concrete `ParticleSpecies` must define the following methods:

| Method | Purpose |
|---|---|
| `graphrep(ps)` | Return the graph representation of the particle species (a `NautyDiGraph`). See the note below. |
| `nsites(ps)` | Number of binding sites. |
| `bindingsites(ps, i)` | Return the `i`th `BindingSite`. |
| `setcolors!(ps, colors)` | Assign new colors to the binding sites. Called by `BindingRules` during construction. |
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

**Requirements.** The graph must have at least 3 vertices and contain a directed cycle. The most common choice is a directed cycle of `n` vertices, one vertex per binding site, with site `i` occupying vertex `i`. Each `BindingSite` records which vertex range it covers via its `vertices` field, so `BindingSite(pose, color, i:i, ...)` places one site on vertex `i`.

**Species with fewer than 3 sites.** If your particle has only 1 or 2 binding sites, use a 4-vertex directed cycle and group the vertices into pairs, each pair representing one binding site. The 2-patch branch of the built-in `PatchyDisk` in `src/species/patchyparticlespecies.jl` is an example of this pattern.

## A worked example: rectangle

Let's define a species for a non-square rectangle with binding sites at the midpoints of its four edges. This is a good template because it exercises every part of the interface without requiring specialized geometry.

```julia
using Roly
using Roly: BindingSite, ParticleSpecies, SpeciesAndPose
import Roly: graphrep, nsites, bindingsites, setcolors!,
             bounding_radius, isconvex, overlap
using NautyGraphs
using Graphs: cycle_digraph
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
function Rectangle(width::Real, height::Real)
    F = float(promote_type(typeof(width), typeof(height)))
    w, h = F(width), F(height)
    tol = sqrt(eps(F)) * max(w, h)

    # Four binding sites at edge midpoints, each facing outward along the edge normal.
    # Sites are ordered clockwise: right, bottom, left, top.
    sites = [
        BindingSite(Pose(SVector{2,F}( w/2, 0),   Angle2d{F}(0)),        1, 1:1, tol, tol),
        BindingSite(Pose(SVector{2,F}(0,  -h/2),  Angle2d{F}(-F(π)/2)),  2, 2:2, tol, tol),
        BindingSite(Pose(SVector{2,F}(-w/2, 0),   Angle2d{F}(F(π))),     3, 3:3, tol, tol),
        BindingSite(Pose(SVector{2,F}(0,   h/2),  Angle2d{F}(F(π)/2)),   4, 4:4, tol, tol),
    ]

    # Directed 4-cycle: one vertex per binding site, all with the same label.
    g = NautyDiGraph(cycle_digraph(4); vertex_labels=fill(Cint(1), 4))

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

`setcolors!` rebuilds each `BindingSite` with a new color while preserving pose, vertex range, and tolerances:

```julia
function setcolors!(ps::Rectangle, colors::AbstractVector{<:Integer})
    length(colors) == nsites(ps) || throw(ArgumentError("wrong number of colors"))
    for k in eachindex(ps.sites)
        s = ps.sites[k]
        ps.sites[k] = BindingSite(s.pose, colors[k], s.vertices,
                                  s.touching_tolerance, s.alignment_tolerance)
    end
    return nothing
end
```

### Overlap

The overlap check decides whether two particles at given poses intersect. For convex polygons the standard approach is the *Separating Axis Theorem* (SAT): two convex shapes are disjoint if and only if their projections onto some axis do not overlap; it suffices to test the axes perpendicular to their edges.

```julia
function overlap(p1::SpeciesAndPose{<:Rectangle}, p2::SpeciesAndPose{<:Rectangle}; kwargs...)
    spcs1, pose1 = p1
    spcs2, pose2 = p2

    skin_sum = spcs1.skin + spcs2.skin
    for (corners, pose) in ((spcs1.corners, pose1), (spcs2.corners, pose2))
        n = length(corners)
        for i in 1:n
            edge = pose.psi * (corners[mod1(i + 1, n)] - corners[i])
            normal = SVector(-edge[2], edge[1])
            lo1, hi1 = extrema(dot(normal, pose1 * c) for c in spcs1.corners)
            lo2, hi2 = extrema(dot(normal, pose2 * c) for c in spcs2.corners)
            (hi2 < lo1 + skin_sum || hi1 < lo2 + skin_sum) && return false
        end
    end
    return true
end
```

## Using the new species

Once the interface is implemented, the new species can be used just like any built-in one:

```julia
rect = Rectangle(2.0, 1.0)
bonds = [1 1 1 3;  # right edge binds to left edge
         1 2 1 4]  # bottom edge binds to top edge
sys = BindingRules(bonds, rect)

polys = polygen(sys; maxsize=6)
```
