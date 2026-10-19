# Saving and loading

```@meta
CurrentModule = Roly
```

Binding rules and polyforms are written to files with [`tofile`](@ref) and read back with [`fromfile`](@ref).

```@example io
using Roly
rules = BindingRules([1 1 2 3], UnitSquare)
polys = filter(p -> nparticles(p) > 0, polygen(rules; maxsize=2))
file = tofile(joinpath(mktempdir(), "dimers.rly"), polys)
print(read(file, String))
```

```@example io
fromfile(file)
```

A file holds a set of binding rules, a single polyform, or a vector of polyforms.
A polyform refers to the species and sites of its rules by number, so a file with polyforms always holds their rules as well.

## Contents

Every file starts with a header:

- `contenttype`: `roly/rulesonly`, `roly/polysingle` or `roly/polylist`, which also marks the file as a Roly file
- `formatversion`: the version of the file format, which [`fromfile`](@ref) checks
- `rolyversion`: the version of Roly that wrote the file

The binding rules are stored as their species and as one row `[species1, site1, species2, site2]` per pair of binding sites that bond, the form the [`BindingRules`](@ref) constructor takes.
A polyform is stored as the species of its particles, in order, and one row `[particle1, site1, particle2, site2, twist]` per bond, where the twist is the one [`twist`](@ref) reports and always `0` in 2D.
Reading a polyform attaches its particles one bond at a time, so its absolute position and orientation are not kept.

## Formats

The format follows the file extension, or the keyword `format`.

| Format | Extensions     | Keyword          |
|:-------|:---------------|:-----------------|
| TOML   | `.rly`, `.toml` | `format=:toml`   |
| JLD2   | `.rlyb`, `.jld2` | `format=:jld2` |

Files with any other extension are written as TOML.
TOML is readable and editable by hand.
JLD2 is binary and needs `using JLD2`.
JLD2 files are HDF5 files laid out with plain groups and arrays, so other HDF5 tools can read them:
a dictionary is a group, a list of equal-length rows a matrix with one column per row, and any other list a group with members `1`, `2`, ....
The polyforms of a `polylist` file are stored as flat columns, `species` with `nparticles` and `bonds` with `nbonds`, which keeps files of many polyforms compact.

## Serializing species

[`tofile`](@ref) writes through [`todict`](@ref) and [`fromfile`](@ref) reads through [`fromdict`](@ref), which turn rules, polyforms and species into a `Dict{String,Any}` of strings, numbers and arrays, and back.
Every built-in species stores the arguments of its most general constructor, and is rebuilt by calling it.

A custom species opts in by defining both methods.
Its dictionary names its type, with the module, in the entry `type`:

```julia
Roly.todict(ps::MySpecies) = Dict{String,Any}("type" => "MyPackage.MySpecies", "size" => ps.size)
Roly.fromdict(::Type{<:MySpecies}, d) = MySpecies(d["size"])
```

A species without these methods is written with only its type, number of sites and colors, which is not enough to read the file back.
