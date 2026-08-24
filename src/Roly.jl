module Roly

using LinearAlgebra, StaticArrays, SparseArrays, Rotations, Statistics, Random
using Base.Iterators, DataStructures
using Graphs, NautyGraphs
using ReverseSearch

# Geometry primitives
export Pose, dimension, numtype, posetype
export Rotation, Angle2d, RotXYZ, RotMatrix3, rotation_angle, rotation_axis, SVector

# Binding sites
export BindingSite, color

# Bodies and their symmetry groups
export Polyhedron
export Tetrahedron, Cube, Octahedron, Dodecahedron, Icosahedron, Pyramid, Prism, Antiprism
export RotationGroup, Cyclic, Dihedral, Tetrahedral, Octahedral, Icosahedral, grouporder

# Particle species
export ParticleSpecies, SpeciesAndPose
export nsites, bindingsites, graphrep, symmetrynumber

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
export PolygonParticleSpecies, UnitNgon, UnitTriangle, UnitSquare, UnitHexagon
export PolyhedronParticleSpecies
export UnitTetrahedron, UnitCube, UnitOctahedron, UnitDodecahedron, UnitIcosahedron
export UnitPyramid, UnitPrism, UnitAntiprism
export PatchyParticleSpecies, PatchyDisk, PatchySphere

# Public, but not exported: reach for these as `Roly.faces(p)`, or import them by name.
# They are stable API, but specific enough to a body, an encoding or a species that putting
# them in every user's namespace is not worth it.

# Polyhedron geometry
public corners, faces, facevertices, nfaces, nedges
public facecentroid, facecentroids, facenormal, facenormals, edgemidpoint
public inradius, minedgelength

# Graph encodings and the symmetry they record
public dartencoding, cycleencoding
public rotationgroup, geometriclabels, facegauge, siteorbits, sitestabilizers
public site_symmetry, check_encoding

# What a bond fixes about relative orientation
public contact_pairing, standard_offset, twistfreedom, bondperiod, nphases, phase

# The `ParticleSpecies` interface, implemented rather than called
public isconvex, bounding_radius, could_contact, overlap
public sat_overlap, edgenormals, polyhedron


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

include("ruleeditor.jl")
export ruleeditor

export render, polyformplot, polyformplot!

"""
    render(p; hidedecorations=true, kwargs...)

Draw a [`Polyform`](@ref) or [`ParticleSpecies`](@ref) into a fresh figure, picking a 2D or 3D
axis to match its dimension. Provided by the Makie extension, so it needs a backend loaded.

3D output needs a backend with a depth buffer, GLMakie or WGLMakie. CairoMakie sorts
primitives instead of depth-testing them, so 3D polyforms show artifacts where faces meet.

`bindingrules=rules` draws sites no bond can use as inert, which works for a bare species too.
"""
function render end

"""
    polyformplot(p; kwargs...)

Makie recipe behind [`render`](@ref), drawing a polyform or species into a new axis.
Provided by the Makie extension.
"""
function polyformplot end

"""
    polyformplot!(ax, p; kwargs...)

Draw a polyform or species onto an axis you already have, the mutating form of
[`polyformplot`](@ref). Provided by the Makie extension.
"""
function polyformplot! end

end