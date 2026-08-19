# Roly.jl

Roly.jl (*Reverse-search Polyform enumerator*) is a Julia package for enumerating polyforms: aggregates formed by connecting arbitrarily shaped building blocks at their binding sites, in 2D or 3D.

Given a set of building blocks and a list of allowed bonds between their binding sites, Roly.jl can:

- exhaustively enumerate every geometrically valid polyform, up to a chosen size,
- generate and store polyforms for further processing,
- estimate the number of polyforms when exact enumeration is too expensive,
- build binding rules interactively by placing particles on a lattice ([`ruleeditor`](workflow.md#Sketching-rules-interactively)).

## Installation

From the Julia REPL, press `]` to enter Pkg mode, then:

```
pkg> add https://github.com/goodrichgroup/Roly.jl
```

Roly.jl currently does not support Windows.

## Where to next

- [Workflow](workflow.md) walks through the basic usage: defining binding rules, coloring a species' binding sites, enumerating, generating, counting, and visualizing polyforms.
- [Custom particle species](custom_species.md) explains how to define your own building block geometries.
- [Orientation and phases](orientation.md) covers what a bond fixes about the two particles' relative orientation, and how to control it. Worth reading if a species does not enumerate what you expected.
- [API reference](api.md) lists every exported name.

## Citation

If you use Roly.jl in your work, please cite:

```
@article{roly2025,
    year = {2025},
    title = {{Accessing Semiaddressable Self-Assembly with Efficient Structure Enumeration}},
    author = {Hübl, Maximilian C. and Goodrich, Carl P.},
    journal = {Physical Review Letters},
    doi = {10.1103/physrevlett.134.058204},
    pages = {058204},
    number = {5},
    volume = {134}
}
```
