# Workflow

## Defining binding rules

A `BindingRules` object holds a list of particle species and the bonds allowed between their binding sites.
A bond is written `[species_i site_i species_j site_j]`, meaning site `site_i` of species `species_i` may bind site `site_j` of species `species_j`.

```jldoctest workflow
julia> using Roly

julia> bonds = [1 3 2 3;
                2 2 3 2;
                2 1 4 1;
                3 1 4 1];

julia> sys = BindingRules(bonds, UnitTriangle)
2d BindingRules[n=4, k=4]
```

Pass a single species if all blocks have the same shape, otherwise a vector with one species per index used in `bonds`.

In 2D Roly ships [`UnitTriangle`](@ref), [`UnitSquare`](@ref), [`UnitHexagon`](@ref) and [`PolygonParticleSpecies`](@ref) for regular polygons, and [`PatchyDisk`](@ref) for a disk with sites on its rim.
In 3D it ships [`UnitTetrahedron`](@ref), [`UnitCube`](@ref), [`UnitOctahedron`](@ref), [`UnitDodecahedron`](@ref), [`UnitIcosahedron`](@ref), [`UnitPyramid`](@ref), [`UnitPrism`](@ref), [`UnitAntiprism`](@ref) and [`PolyhedronParticleSpecies`](@ref) for any convex polyhedron with one site per face, and [`PatchySphere`](@ref) for a sphere whose patches inherit a polyhedron's rotation group.
See [Custom particle species](custom_species.md) to define your own.

## Colors and symmetry

You describe a species by *coloring* its binding sites, and the bond table refers to those colors.
Which sites the particle cannot tell apart follows from the coloring and the geometry, and Roly derives it.

```jldoctest workflow
julia> symmetrynumber(PolyhedronParticleSpecies(Cube()))                  # every face distinct
1

julia> symmetrynumber(PolyhedronParticleSpecies(Cube(); colors=fill(1, 6)))  # all faces alike
24

julia> caps = [abs(facenormal(Cube(), i)[3]) > 0.5 ? 2 : 1 for i in 1:6];

julia> symmetrynumber(PolyhedronParticleSpecies(Cube(); colors=caps))     # caps apart from sides
8
```

Faces come in no particular order, so the third example picks the caps by their normals.
The answer 8 is `D_4`: telling two opposite faces apart leaves the 4-fold axis through them and the 2-fold axes across it.

Two more keywords say what a bond at a face *means*, rather than which bonds exist.
`locking` decides whether a site holds its partner in the orientation its frame names, and `twists` turns a site about its normal to pick which orientation that is.
See [Orientation and phases](orientation.md).

## Symmetry groups

Bodies are built by name: [`Cube`](@ref), [`Prism`](@ref), [`Antiprism`](@ref) and the rest all return a [`Polyhedron`](@ref).
[`rotationgroup`](@ref) lists the rotations a body has, and [`grouporder`](@ref) counts those of a named group.

```jldoctest workflow
julia> length(rotationgroup(Cube()))
24

julia> grouporder(Octahedral())
24

julia> grouporder(Dihedral(5))
10
```

The named groups are [`Cyclic`](@ref), [`Dihedral`](@ref), [`Tetrahedral`](@ref), [`Octahedral`](@ref) and [`Icosahedral`](@ref), the rotation-only point groups a rigid body can have.

## Sketching rules interactively

[`ruleeditor`](@ref) builds a bond table geometrically: place blocks on a lattice and it reads the rules off every pair of touching sites.

```julia
sys = ruleeditor(UnitSquare)  # also works for UnitTriangle and UnitHexagon
```

Arrow keys move the cursor, `Enter` places, `Space` erases, `r` and `R` rotate, digits `1` to `9` switch species, `q` accepts.

```
┌ Editor ──────────┐┌ Particles ──────────────────────────┐
│── Species ──     ││ ┌─────┐ ┌─────┐                     │
│▶ ■  species 1    ││ │     │ │     │                     │
│  ■  species 2    ││ │  →  │ │  ←  │                     │
│                  ││ │     │ │     │                     │
│── Keys ──────    ││ └─────┘ └─────┘                     │
│arrows  move      ││ ┌─────┐ ┌─────┐                     │
│enter   place     ││ │     │ │     │                     │
│space   erase     ││ │  ↑  │ │  ↑  │                     │
│r / R   rotate    ││ │     │ │     │                     │
│1-9     species   ││ └─────┘ └─────┘                     │
│c       clear     ││                                     │
│q       accept    ││                                     │
└──────────────────┘└─────────────────────────────────────┘
```

Pass `output=:bonds` or `output=:matrix` for a copy-pasteable result instead of a `BindingRules`:

```julia
bonds = ruleeditor(UnitSquare; output=:bonds)  # n×4 integer matrix
sys   = BindingRules(bonds, UnitSquare)         # reproduces the same rules
```

## Enumerating polyforms

`polyenum` walks every polyform the rules allow.

```jldoctest workflow
julia> result = polyenum(sys; maxsize=20, maxstrs=100_000);

julia> result.nstructures
16

julia> result.largest_size
5

julia> result.status
Finished::RSStatus = 0
```

Cap `maxsize` (particles per polyform) or `maxstrs` (total polyforms) when the rules allow unbounded growth.
`status` says why the run stopped: `Finished`, `MaxDepthReached`, `MaxVerticesReached` or `BreakTriggered`.

## Storing polyforms

`polygen` returns the polyforms in a `Vector`, sorted by size.

```jldoctest workflow
julia> polys = polygen(sys; maxsize=20);

julia> length(polys)
16
```

## Counting polyforms

`countpolyforms` counts without storing anything, switching to an unbiased sampled estimate when exact enumeration gets too expensive.

```jldoctest workflow
julia> c = countpolyforms(sys);

julia> c.n
16.0

julia> c.exact
true

julia> c.uncertainty
0.0
```

It requires an explicit `maxsize` when the rules allow polyforms of unbounded size.

## Applying constraints

`polyenum` takes a callback that runs at each polyform, receiving it and its size and returning one of three signals:

- `ACCEPT` counts the polyform and keeps exploring it,
- `REJECT` skips the polyform and everything grown from it,
- `BREAK` stops the enumeration.

```jldoctest workflow
julia> constraint(s, _) = composition(s)[4] <= 1 ? ACCEPT : REJECT;

julia> polyenum(constraint, sys).nstructures
14
```

`REJECT` prunes a whole subtree, so the constraint must be monotone.
If a polyform violates it, everything grown from it must violate it too.

## Visualizing polyforms

Load any [Makie](https://docs.makie.org) backend to activate the plotting extension.

```julia
using GLMakie      # or CairoMakie, WGLMakie, ...
render(polys[end]) # display a single polyform
```

[`render`](@ref) picks a 2D or 3D axis to match the polyform.
**Use GLMakie or WGLMakie for 3D**, since CairoMakie sorts primitives instead of depth-testing them and shows seams where faces meet.
2D output is fine in any backend.

A species renders too, which is the quickest way to see how its faces are colored and which way its sites face.

```julia
render(PolyhedronParticleSpecies(Prism(3); colors=[1, 2, 2, 2, 1]))
render(UnitCube; bindingrules=sys)   # sites no bond can use are drawn inert
```

[`polyformplot!`](@ref)`(ax, poly)` draws onto an existing Makie axis.
