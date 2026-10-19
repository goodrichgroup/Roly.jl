# Workflow

```@meta
CurrentModule = Roly
```

## Defining binding rules

A `BindingRules` object holds a list of particle species and the bonds allowed between their binding sites.
A bond is written `[species_i site_i species_j site_j]`, meaning site `site_i` of species `species_i` may bind site `site_j` of species `species_j`.

```jldoctest workflow
julia> using Roly

julia> bonds = [1 3 2 3;
                2 2 3 2;
                2 1 4 1;
                3 1 4 1];

julia> rules = BindingRules(bonds, UnitTriangle)
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

julia> caps = [abs(Roly.facenormal(Cube(), i)[3]) > 0.5 ? 2 : 1 for i in 1:6];

julia> symmetrynumber(PolyhedronParticleSpecies(Cube(); colors=caps))     # caps apart from sides
8
```

Faces come in no particular order, so the third example picks the caps by their normals.
The answer 8 is `D_4`: telling two opposite faces apart leaves the 4-fold axis through them and the 2-fold axes across it.

Two more keywords say what a bond at a face *means*, rather than which bonds exist.
`locking` decides whether a site holds its partner in the orientation its frame names, and `twists` turns a site about its normal to pick which orientation that is.
See [Orientation and twists](orientation.md).

## Symmetry groups

Bodies are built by name: [`Cube`](@ref), [`Prism`](@ref), [`Antiprism`](@ref) and the rest all return a [`Polyhedron`](@ref).
[`rotationgroup`](@ref) lists the rotations a body has, and [`grouporder`](@ref) counts those of a named group.

```jldoctest workflow
julia> length(Roly.rotationgroup(Cube()))
24

julia> grouporder(Octahedral())
24

julia> grouporder(Dihedral(5))
10
```

The named groups are [`Cyclic`](@ref), [`Dihedral`](@ref), [`Tetrahedral`](@ref), [`Octahedral`](@ref) and [`Icosahedral`](@ref), the rotation-only point groups a rigid body can have.

## Enumerating polyforms

Once you have defined a set of binding rules, use `polyenum` to enumerate all allowed polyforms.

```jldoctest workflow
julia> result = polyenum(rules);

julia> result.n
16

julia> result.maxsize
5

julia> result.status
Complete::EnumerationStatus = 0
```

This system allows 16 polyforms, and `polyenum` enumerates all of them.
In general, binding rules may allow infinitely many polyforms.
Before an enumeration, Roly automatically checks if the set of allowed polyforms is finite, which can also be done manually with [`isbounded`](@ref); see [Bounded and unbounded rules](growth.md).
It is also possible to enumerate only to a certain size of polyforms (`maxsize`) or to limit the total number of polyforms enumerated (`maxpolys`).

```jldoctest workflow
julia> isbounded(rules)
true

julia> polyenum(rules; maxsize=3).status
MaxSizeReached::EnumerationStatus = 1
```

`status` says why the enumeration stopped: `Complete`, `MaxSizeReached`, `MaxPolysReached` or `Interrupted`.

## Storing polyforms

To store all polyforms in memory for further processing, use `polygen`, which returns a list sorted by size.

```jldoctest workflow
julia> polys = polygen(rules);

julia> length(polys)
16
```

## Counting polyforms

To count polyforms without storing them, or to estimate the number of polyforms when full enumeration is too expensive, use `polycount`.

```jldoctest workflow
julia> result = polycount(rules);

julia> result.n
16

julia> result.exact
true

julia> result.err
0.0
```

`polycount` enumerates exactly up to a configurable number of polyforms and switches to estimation beyond that.
Like `polyenum`, it first checks that the set of allowed polyforms is finite; an unbounded set of rules needs a `maxsize`.

## Incorporating constraints

It is often desirable to impose additional constraints on generated polyforms.
`polyenum` takes a callback that runs at each polyform, receiving it and its size and returning one of three signals:

- `ACCEPT` counts the polyform and keeps exploring it,
- `REJECT` skips the polyform and everything grown from it,
- `BREAK` stops the enumeration.

For example, to enumerate only polyforms with at most one particle of species 4:

```jldoctest workflow
julia> constraint(s, n) = composition(s)[4] <= 1 ? ACCEPT : REJECT;

julia> polyenum(constraint, rules).n
14
```

To ensure well-defined behavior, if a polyform `s` violates the constraint, all larger polyforms that can be generated by adding particles to `s` must also violate the constraint.

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
render(UnitCube; bindingrules=rules)   # sites no bond can use are drawn inert
```

[`polyformplot!`](@ref)`(ax, poly)` draws onto an existing Makie axis.
