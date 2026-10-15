# Bounded and unbounded rules

Some rule sets run out of structures and others grow forever.
[`isbounded`](@ref) settles which, and [`tilings`](@ref) finds the periodic structures an unbounded set forms.

The figures below use CairoMakie because the documentation builds headless.
For 3D of your own, prefer GLMakie or WGLMakie, which depth-test instead of sorting.

```@setup growth
using Roly, CairoMakie, LinearAlgebra, Rotations
CairoMakie.activate!(type = "png")
```

## A bounded system runs out

Squares whose site 1 binds site 2 can only turn the same way at every bond, so a chain of them folds back on itself.

```@example growth
ring = BindingRules([1 1 1 2], UnitSquare)
isbounded(ring)
```

`true` here is a proof rather than a search giving up.
A chain long enough must repeat one of finitely many states, and the motion between the two repeats is then applied forever.
Here that motion is a quarter turn, and every copy of a pure rotation stays the same distance from its axis, so infinitely many of them cannot fit.
The whole system is therefore finite and [`countpolyforms`](@ref) counts it exactly.

```@example growth
countpolyforms(ring)
```

Four structures, and the largest closes into a ring with no open site left.

```@example growth
polys = polygen(ring)
biggest = polys[argmax(nparticles.(polys))]
nparticles(biggest), nbonds(biggest), length(opensitelocs(biggest))
```

```@example growth
render(biggest)
```

## An unbounded system, and its periodic structures

Binding opposite sides of a square instead gives the square lattice.

```@example growth
square = BindingRules([1 2 1 4; 1 3 1 1], UnitSquare)
isbounded(square)
```

A single square repeats under a pure translation, which is the simplest way to be unbounded.
It is also a tiling: one particle whose open sites all close onto translated copies of itself.
[`tilings`](@ref) enumerates every periodic structure whose cell holds at most `maxsize` particles.
Each one knows its cell, the translations its periodic bonds name, which of the rules' bonds it spends, and whether it leaves any site open.

```@example growth
ts = tilings(square; maxsize=1)
complete = filter(iscomplete, ts)
first(complete)
```

```@example growth
translations(first(complete))
```

Two translations, so the cell tiles the plane.
`tilings` also returns the partial closures, with one translation only, which are infinite strips rather than tilings, so filter on `iscomplete` when you want the plane covered.

```@example growth
length(ts), length(complete)
```

The cell is the repeat unit, and a tiling is drawn as a patch of it.

```@example growth
render(first(complete))
```

A cell may need several particles.
Squares whose site 1 binds site 2 and site 3 binds site 4 turn their neighbour by a quarter, so no translation carries one onto a bonded copy, and the cell holds two of them.

```@example growth
turn = BindingRules([1 1 1 2; 1 3 1 4], UnitSquare)
length(tilings(turn; maxsize=1)), length(tilings(turn; maxsize=2))
```

Every structure is reported once, by its smallest cell.
Raising `maxsize` adds cells that genuinely need more particles and never a supercell of one already found.

## Counting what cannot be enumerated

An unbounded system has infinitely many structures, so [`countpolyforms`](@ref) needs a `maxsize`.
Up to that size it enumerates exactly while the count is small, and switches to sampling the reverse-search tree beyond `exact_budget`.

```@example growth
countpolyforms(square; maxsize=9)
```

The `≈` and the error bar mark an estimate rather than a count.

## Unbounded without ever repeating

In the plane, growing forever and closing periodically are the same thing: a repeated motion is a translation or a rotation, and only a translation escapes.
In space a motion can be a screw, a turn about an axis together with a climb along it, and a screw escapes whether or not its turn ever adds up to a whole one.
Regular tetrahedra are the standard example.
Glued face to face they form the Boerdijk-Coxeter helix, whose twist per tetrahedron is `acos(-2/3)`, an irrational fraction of a turn.

```@example growth
tetrahedra = PolyhedronParticleSpecies(Tetrahedron(); colors=fill(1, 4))
helixrules = BindingRules([1 1 1 1], tetrahedra)
isbounded(helixrules)
```

Copies far enough apart along the axis can never touch, which is what makes the growth infinite, and no multiple of the angle is ever a full turn, so the chain never closes.
This is the case a search for periodic closures cannot decide, and the reason [`isbounded`](@ref) looks at the repeating motion itself.

The helix is visible in the structures.
Among the chains of eight tetrahedra, the longest reaches well short of eight straight steps, because it is curling.

```@example growth
function ischain(p)
    degree = zeros(Int, nparticles(p))
    for (a, b) in bonds(p)
        degree[a.particle] += 1
        degree[b.particle] += 1
    end
    return all(<=(2), degree)
end

chains = filter(ischain, polygen(helixrules; maxsize=8))
helix = chains[argmax(nparticles.(chains))]
xs = [p.pose.x for p in helix.particles]
norm(xs[end] - xs[1]) / (nparticles(helix) - 1)
```

```@example growth
render(helix)
```
