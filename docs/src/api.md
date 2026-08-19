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

### Bonds and phases

How many distinct ways a partner may attach at a site, and which one a given attachment is.
See [Orientation and phases](orientation.md) for what the numbers mean.

```@docs
twistfreedom
bondperiod
nphases
phase
standard_offset
contact_pairing
```

## Polyhedra and graph encodings

```@docs
Polyhedron
corners
faces
facevertices
nfaces
nedges
facecentroid
facenormal
edgemidpoint
minedgelength
inradius
rotationgroup
geometriclabels
facegauge
siteorbits
site_symmetry
check_encoding
sitestabilisers
dartencoding
cycleencoding
```

### Rotation groups

```@docs
RotationGroup
Cyclic
Dihedral
Tetrahedral
Octahedral
Icosahedral
grouporder
```

### Solids

```@docs
Tetrahedron
Cube
Octahedron
Dodecahedron
Icosahedron
Pyramid
Prism
Antiprism
```

## Particle species

```@docs
ParticleSpecies
SpeciesAndPose
bindingsites
nsites
graphrep
isconvex
symmetrynumber
bounding_radius
could_contact
overlap
sat_overlap
edgenormals
```

### Built-in species

```@docs
PolygonParticleSpecies
UnitTriangle
UnitSquare
UnitHexagon
PolyhedronParticleSpecies
shape
UnitTetrahedron
UnitCube
UnitOctahedron
UnitDodecahedron
UnitIcosahedron
UnitPyramid
UnitPrism
UnitAntiprism
PatchyParticleSpecies
PatchyDisk
PatchySphere
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

An enumeration reports why it stopped as an `RSStatus`: `Finished`, `MaxDepthReached`, `MaxVerticesReached` or `BreakTriggered`.
A callback returns `ACCEPT`, `REJECT` or `BREAK`; see [Applying constraints](workflow.md#Applying-constraints).

## Visualization

Provided by the Makie extension, so a backend has to be loaded.

```@docs
render
polyformplot
polyformplot!
```

## Interactive rule editor

```@docs
ruleeditor
```
