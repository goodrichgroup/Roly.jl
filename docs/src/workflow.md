# Workflow

## Defining binding rules

A `BindingRules` object holds a list of particle species and a set of allowed bonds between their binding sites. A bond is specified in the form `[species_i site_i species_j site_j]`, meaning that site `site_i` of species `species_i` may bind to site `site_j` of species `species_j`.

```jldoctest workflow
julia> using Roly

julia> bonds = [1 3 2 3;
                2 2 3 2;
                2 1 4 1;
                3 1 4 1];

julia> sys = BindingRules(bonds, UnitTriangle)
2d BindingRules[n=4, k=4]
```

If all building blocks share the same shape, pass a single species. Otherwise, pass a vector of species — one per row-index used in `bonds`.

Roly ships with a few built-in species. In 2D: [`UnitTriangle`](@ref), [`UnitSquare`](@ref), [`UnitHexagon`](@ref) and [`PolygonParticleSpecies`](@ref) (regular polygons), and [`PatchyDisk`](@ref) (a disk with binding sites on its rim). In 3D: [`UnitTetrahedron`](@ref), [`UnitCube`](@ref), [`UnitOctahedron`](@ref), [`UnitDodecahedron`](@ref), [`UnitIcosahedron`](@ref), [`UnitPyramid`](@ref), [`UnitPrism`](@ref), [`UnitAntiprism`](@ref) and [`PolyhedronParticleSpecies`](@ref) (convex polyhedra with one binding site per face), and [`PatchySphere`](@ref) (a sphere whose patches inherit a polyhedron's rotation group). See the [custom species page](custom_species.md) to define your own.

## Colors, and the symmetry that follows from them

A species is described by *coloring* its binding sites. Colors are what the bond table refers to, and they are the whole statement: which sites the particle cannot tell apart follows from the coloring together with the geometry, and Roly derives it.

```jldoctest workflow
julia> symmetrynumber(PolyhedronParticleSpecies(Cube()))                  # every face distinct
1

julia> symmetrynumber(PolyhedronParticleSpecies(Cube(); colors=fill(1, 6)))  # all faces alike
24

julia> caps = [abs(facenormal(Cube(), i)[3]) > 0.5 ? 2 : 1 for i in 1:6];

julia> symmetrynumber(PolyhedronParticleSpecies(Cube(); colors=caps))     # caps apart from sides
8
```

A solid's faces come in no particular order, so the third example picks its two caps by their normals rather than by index. The result is `D_4`: distinguishing two opposite faces leaves the 4-fold axis through them and the 2-fold axes across it. There is no `labels` keyword to override this — a labelling claiming more symmetry than the shape has would merge structures that are genuinely different, silently, so it is derived by [`siteorbits`](@ref) rather than given.

Two further keywords control what a bond at a face *means*, rather than which bonds exist:

- `locking` — whether a site holds its partner in the orientation its frame names. The default, `true`, gives one attachment unless the particle's own symmetry makes the frame ambiguous. `false` admits every orientation the face geometrically permits.
- `twists` — an angle in radians turning a face's site about its own normal, which selects *which* relative orientation the bond means.

Both take one value for the species or one per face. See [Orientation and phases](orientation.md).

## Symmetry groups

A body is always built by naming it — [`Cube`](@ref), [`Prism`](@ref), [`Antiprism`](@ref) and the rest all return a [`Polyhedron`](@ref).
A rotation group cannot stand in for one, since duality gives the cube and the octahedron the same group, as it does the dodecahedron and the icosahedron.

The groups instead report symmetry: [`rotationgroup`](@ref) lists the rotations a body has, and [`grouporder`](@ref) counts those of a named group.

```jldoctest workflow
julia> length(rotationgroup(Cube())) == length(rotationgroup(Octahedron()))
true

julia> grouporder(Octahedral())
24

julia> grouporder(Dihedral(5))
10
```

The groups are [`Cyclic`](@ref), [`Dihedral`](@ref), [`Tetrahedral`](@ref), [`Octahedral`](@ref) and [`Icosahedral`](@ref) — the proper (rotation-only) point groups, which are the only ones a rigid body can have.

## Sketching rules interactively

Instead of writing the bonds matrix by hand, you can construct one geometrically with [`ruleeditor`](@ref) — a small terminal editor that lets you place building blocks on a lattice and infers the binding rules from every pair of touching sites:

```julia
sys = ruleeditor(UnitSquare)  # also works for UnitTriangle and UnitHexagon
```

Arrow keys move the cursor, `Enter` places a particle, `Space` erases, `r`/`R` rotate, digits `1`–`9` add and switch between species, `q` accepts. A typical session looks like this:

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

Pass `output=:bonds` or `output=:matrix` to get a copy-pasteable representation instead of a `BindingRules`, useful for pinning a specific design in code:

```julia
bonds = ruleeditor(UnitSquare; output=:bonds)  # n×4 integer matrix
sys   = BindingRules(bonds, UnitSquare)         # reproduces the same rules
```

## Enumerating polyforms

`polyenum` walks through every polyform allowed by a set of binding rules:

```jldoctest workflow
julia> result = polyenum(sys; maxsize=20, maxstrs=100_000);

julia> result.nstructures
16

julia> result.largest_size
5

julia> result.status
Finished::RSStatus = 0
```

For unbounded systems, always cap either `maxsize` (particles per polyform) or `maxstrs` (total number of polyforms). `status` reports whether the enumeration ran to completion (`Finished`, `MaxDepthReached`, `MaxVerticesReached`, or `BreakTriggered`).

## Storing polyforms

`polygen` returns all enumerated polyforms in a `Vector`, sorted by size:

```jldoctest workflow
julia> polys = polygen(sys; maxsize=20);

julia> length(polys)
16
```

## Counting polyforms

`countpolyforms` returns a count without storing individual polyforms. When exact enumeration is expensive, it switches to an unbiased importance-sampling estimate:

```jldoctest workflow
julia> c = countpolyforms(sys);

julia> c.n
16.0

julia> c.exact
true

julia> c.uncertainty
0.0
```

`countpolyforms` requires an explicit `maxsize` when the binding rules allow polyforms of unbounded size.

## Applying constraints

`polyenum` accepts a callback that decides what to do at each visited polyform. The callback receives the current polyform and its size, and returns `ACCEPT`, `REJECT`, or `BREAK`:

- `ACCEPT` — count this polyform and keep exploring its extensions,
- `REJECT` — skip this polyform *and everything that grows from it*,
- `BREAK` — stop the enumeration entirely.

```jldoctest workflow
julia> constraint(s, _) = composition(s)[4] <= 1 ? ACCEPT : REJECT;

julia> polyenum(constraint, sys).nstructures
14
```

Because `REJECT` prunes the entire subtree, the constraint must be monotone: if a polyform violates it, every larger polyform grown from it must also violate it.

## Visualizing polyforms

Roly ships a [Makie](https://docs.makie.org) extension. Load any Makie backend to activate it:

```julia
using GLMakie      # or CairoMakie, WGLMakie, ...
render(polys[end]) # display a single polyform
```

[`render`](@ref) picks a 2D or 3D axis to match the polyform's dimension. **For 3D, use a backend with a depth buffer — GLMakie or WGLMakie.** CairoMakie sorts primitives instead of depth-testing them, so a 3D polyform renders with visible artifacts where faces meet: faint seams inside transparent faces and a jagged fringe where a patch interpenetrates its body. 2D output is fine in any backend.

A species renders too, which is the quickest way to see how its faces are colored and which way its sites face:

```julia
render(PolyhedronParticleSpecies(Prism(3); colors=[1, 2, 2, 2, 1]))
render(UnitCube; bindingrules=sys)   # sites no bond can use are drawn inert
```

[`polyformplot!`](@ref)`(ax, poly)` draws onto an existing Makie axis.
