# API reference

## Poses

```@docs
Pose
dimension
numtype
posetype
```

## Binding sites

```@docs
BindingSite
color
```

## Particle species

```@docs
ParticleSpecies
bindingsites
nsites
graphrep
isconvex
symmetrynumber
bounding_radius
could_contact
overlap
```

### Built-in species

```@docs
PolygonParticleSpecies
UnitTriangle
UnitSquare
UnitHexagon
PatchyParticleSpecies
PatchyDisk
```

## Binding rules

```@docs
BindingRules
interactionmatrix
nspecies
nbonds
ncolors
species
bonded_colors
bonded_sites
bonded_species
isinert
```

## Polyforms

```@docs
Polyform
nparticles
bindingrules
composition
bonds
bondindex
interior_edges
exterior_edges
```

## Enumeration

```@docs
polyenum
polygen
countpolyforms
PolyformCount
```
