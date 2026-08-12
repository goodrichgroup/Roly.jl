module Roly

using LinearAlgebra, StaticArrays, SparseArrays, Rotations, Statistics, Random
using Base.Iterators, DataStructures
using Graphs, NautyGraphs
using ReverseSearch

# Geometry primitives
export Pose, dimension, numtype, posetype
export Rotation, Angle2d, RotXYZ, rotation_angle, rotation_axis, SVector

# Binding sites
export BindingSite, color

# Particle species
export ParticleSpecies, SpeciesAndPose
export nsites, bindingsites, graphrep, isconvex, symmetrynumber
export bounding_radius, could_contact, overlap

# Assembly system
export BindingRules, interactionmatrix
export ncolors, nspecies, nbonds, bonded_colors, bonded_sites, bonded_species, isinert, species

# Polyforms
export Polyform, nparticles, bindingrules, composition
export bonds, bondindex, interior_edges, exterior_edges, subpolyform

# Environments
export Environment, ParticleEnvironment, BondEnvironment
export particleenvironments, bondenvironments, crop
export tilings, isunitcell, tilelatticevectors, cantile

# Enumeration
export ACCEPT, REJECT, BREAK
export RSStatus, Finished, MaxVerticesReached, MaxDepthReached, BreakTriggered
export polyenum, polygen, countpolyforms, PolyformCount

# Species
export PolygonParticleSpecies, UnitTriangle, UnitSquare, UnitHexagon
export PatchyParticleSpecies, PatchyDisk


include("utils.jl")
include("pose.jl")
include("bindingsite.jl")
include("particlespecies.jl")
include("bindingrules.jl")
include("particle.jl")
include("polyform.jl")
include("enumeration.jl")
include("environments.jl")
include("tiling.jl")

include("species/polygonparticlespecies.jl")
include("species/patchyparticlespecies.jl")

include("ruleeditor.jl")
export ruleeditor

export render, polyformplot, polyformplot!

function render end
function polyformplot end
function polyformplot! end

end