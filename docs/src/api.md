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
rotationgroup
geometriclabels
dartencoding
cycleencoding
sitegraph
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
