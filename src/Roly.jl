module Roly

using LinearAlgebra, StaticArrays, SparseArrays, Rotations, Statistics, Random
using Base.Iterators, DataStructures
using Graphs, NautyGraphs
using ReverseSearch

# Geometry primitives
export Pose, dimension, numtype, posetype
export Rotation, Angle2d, RotXYZ, RotMatrix3, rotation_angle, rotation_axis, SVector

# Binding sites
export BindingSite, color, contact_pairing

# Polyhedra and graph encodings
export Polyhedron, corners, faces, facevertices, nfaces, nedges
export facecentroid, facecentroids, facenormal, facenormals, edgemidpoint
export dartencoding, cycleencoding
export rotationgroup, geometriclabels, facegauge, siteorbits, site_symmetry
export Tetrahedron, Cube, Octahedron, Dodecahedron, Icosahedron, Pyramid, Prism, Antiprism
export RotationGroup, Cyclic, Dihedral, Tetrahedral, Octahedral, Icosahedral, grouporder

# Particle species
export ParticleSpecies, SpeciesAndPose
export nsites, bindingsites, graphrep, isconvex, symmetrynumber
export bounding_radius, could_contact, overlap

# Assembly system
export BindingRules, interactionmatrix
export ncolors, nspecies, nbonds, bonded_colors, bonded_sites, bonded_species, isinert, species

# Polyforms
export Polyform, nparticles, bindingrules, composition
export bonds, bondindex, interior_edges, exterior_edges

# Enumeration
export ACCEPT, REJECT, BREAK
export RSStatus, Finished, MaxVerticesReached, MaxDepthReached, BreakTriggered
export polyenum, polygen, countpolyforms, PolyformCount

# Species
export PolygonParticleSpecies, UnitTriangle, UnitSquare, UnitHexagon
export PolyhedronParticleSpecies, shape
export UnitTetrahedron, UnitCube, UnitOctahedron, UnitDodecahedron, UnitIcosahedron
export UnitPyramid, UnitPrism, UnitAntiprism
export PatchyParticleSpecies, PatchyDisk, PatchySphere


include("utils.jl")
include("pose.jl")
include("bindingsite.jl")
include("particlespecies.jl")
include("encoding.jl")
include("bindingrules.jl")
include("particle.jl")
include("polyform.jl")
include("enumeration.jl")

include("species/polygonparticlespecies.jl")
include("species/polyhedronparticlespecies.jl")
include("species/patchyparticlespecies.jl")

export render, polyformplot, polyformplot!

function render end
function polyformplot end
function polyformplot! end

end