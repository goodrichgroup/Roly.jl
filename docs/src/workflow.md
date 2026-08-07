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

Roly ships with a few built-in species: [`UnitTriangle`](@ref), [`UnitSquare`](@ref), [`UnitHexagon`](@ref) (regular polygons), [`PolygonParticleSpecies`](@ref) (arbitrary regular polygons), and [`PatchyParticleSpecies`](@ref)/[`PatchyDisk`](@ref) (disks or spheres with binding sites on their surface). See the [custom species page](custom_species.md) to define your own.

## Sketching rules interactively

Instead of writing the bonds matrix by hand, you can construct one geometrically with [`ruleeditor`](@ref) — a small terminal editor that lets you place particles on a lattice and infers the binding rules from every pair of touching sites:

```julia
sys = ruleeditor(UnitSquare)  # also works for UnitTriangle and UnitHexagon
```

Each species is drawn in its own color (matching the Makie extension's palette) and the arrow inside each tile shows where site 1 points. Arrow keys move the cursor, `Enter` places a particle, `Space` erases, `r`/`R` rotate, digits `1`–`9` add and switch between species, `q` accepts. A typical session with two adjacent squares of species 1 and one of species 2 looks like:

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

Roly ships a [Makie](https://docs.makie.org) extension for 2D visualization. Load any Makie backend to activate it:

```julia
using GLMakie      # or CairoMakie, WGLMakie, ...
render(polys[end]) # display a single polyform
```

`polyformplot!(ax, poly)` draws a polyform onto an existing Makie axis.
