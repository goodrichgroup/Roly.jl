# Roly.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://goodrichgroup.github.io/Roly.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://goodrichgroup.github.io/Roly.jl/dev/)
[![Build Status](https://github.com/goodrichgroup/Roly.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/goodrichgroup/Roly.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/goodrichgroup/Roly.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/goodrichgroup/Roly.jl)


Roly.jl (_<ins>R</ins>everse-Search P<ins>oly</ins>form Enumerator_) is a Julia package for the enumeration of arbitrary polyforms via [reverse search](https://en.wikipedia.org/wiki/Reverse-search_algorithm). It makes it possible to exhaustively enumerate polyforms (aggregates formed by connecting arbitrarily shaped building blocks at their binding sites) in 2D or 3D and provides an interface to define your own building block geometries and binding rules. Roly.jl is under active development, and breaking changes can occur at any time. Because of its dependencies, Roly.jl currently does not work on Windows.

## Installation
To install Roly.jl directly from your Julia REPL, first press `]` to enter Pkg mode, and then run
```
pkg> add https://github.com/goodrichgroup/Roly.jl
```

## Basic Usage
Enumeration in Roly.jl starts from a `BindingRules` object, which is a list of building block geometries together with an interaction matrix that specifies which binding sites are allowed to bind to each other.
The allowed polyforms can then be enumerated with `polyenum`, generated and stored with `polygen`, or counted (exactly or approximately) with `countpolyforms`.

### Defining Binding Rules
To illustrate the basic process, let's construct a system consisting of four species of triangular building blocks. Binding rules are defined as a list of bonds, where every bond is specified in the form `[species_i site_i species_j site_j]`. For example, `[1 3 2 3]` indicates that site 3 of species 1 is allowed to bind to site 3 of species 2. Roly already comes with definitions for simple polygonal particle geometries (e.g. `UnitTriangle`, `UnitSquare`, `UnitHexagon`), convex polyhedra (e.g. `UnitCube`, `UnitPrism(n)`), as well as patchy particles (e.g. `PatchyDisk`, `PatchySphere`). The `BindingRules` constructor takes either a list of geometries or a single geometry if all building blocks are identically shaped.
```julia
using Roly

bonds = [1 3 2 3;
         2 2 3 2;
         2 1 4 1;
         3 1 4 1]
sys = BindingRules(bonds, UnitTriangle)
```

The [documentation](https://goodrichgroup.github.io/Roly.jl/dev/workflow/) lists every built-in geometry and shows how coloring a species' binding sites sets its symmetry. In 3D a bond also fixes how the two blocks are turned relative to one another, which the [orientation page](https://goodrichgroup.github.io/Roly.jl/dev/orientation/) explains. To implement your own particle species, see [custom particle species](https://goodrichgroup.github.io/Roly.jl/dev/custom_species/).

### Sketching Binding Rules interactively
Instead of writing the bonds matrix out by hand, you can build one geometrically with `ruleeditor`, a terminal editor that lets you place building blocks on a lattice and reads the binding rules off every pair of touching sites:
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
sys = BindingRules(bonds, UnitSquare)           # reproduces the same rules
```

### Enumeration
Once you have defined a set of binding rules, use `polyenum` to enumerate all allowed polyforms:
```julia
result = polyenum(sys; maxsize=20, maxstrs=100_000)
result.nstructures   # number of polyforms found
result.largest_size  # size of the largest polyform found
result.status        # Finished, MaxDepthReached, MaxVerticesReached, or BreakTriggered
```
The simple system we have chosen here only allows 16 different polyforms to form. In general however, the number of polyforms might be unbounded, so it is advisable to always impose either a maximal size (`maxsize`) or a maximal count (`maxstrs`). To store all polyforms in memory for further processing, use `polygen`, which returns a list sorted by size:
```julia
strs = polygen(sys; maxsize=20, maxstrs=100_000)
```

### Counting
To count polyforms without storing them, or to estimate when full enumeration is too expensive, use `countpolyforms`:
```julia
c = countpolyforms(sys)
c.n            # count (exact or estimated mean)
c.exact        # true if the count is exact
c.uncertainty  # standard error of the estimate (0 if exact)
```
`countpolyforms` enumerates exactly up to a configurable budget and switches to importance-sampled estimation beyond it. Pass `maxsize` for systems that allow unbounded growth.

### Incorporating constraints
It is often desirable to impose additional constraints on generated polyforms. For example, to enumerate only polyforms with at most one particle of species 4:
```julia
constraint(s, n) = composition(s)[4] <= 1 ? ACCEPT : REJECT
polyenum(constraint, sys)
```
Warning: To ensure well-defined behavior, if a polyform `s` violates the constraint, all larger polyforms that can be generated by adding particles to `s` must also violate the constraint.

### Visualization
Roly.jl provides a [Makie](https://docs.makie.org) extension. Load any Makie backend to activate it:
```julia
using GLMakie  # or CairoMakie, WGLMakie, ...
render(s)      # display a single polyform, or a species
```
`render` picks a 2D or 3D axis to match. For 3D use GLMakie or WGLMakie, since CairoMakie sorts primitives rather than depth-testing them. `polyformplot!` can be used to draw onto an existing Makie axis.

## Citation
If you use Roly.jl in your work, please cite [our paper](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.134.058204) below:
```
@article{roly2025, 
         year = {2025}, 
         title = {{Accessing Semiaddressable Self-Assembly with Efficient Structure Enumeration}}, 
         author = {Hübl, Maximilian C. and Goodrich, Carl P.}, 
         journal = {Physical Review Letters}, 
         issn = {0031-9007}, 
         doi = {10.1103/physrevlett.134.058204}, 
         pages = {058204}, 
         number = {5}, 
         volume = {134}, 
         keywords = {}
}
```
