module Roly

using LinearAlgebra, StaticArrays, SparseArrays, Rotations
using Base.Iterators, DataStructures
using Graphs, NautyGraphs
using ReverseSearch

export Rotation, Angle2d, RotXYZ, rotation_angle, rotation_axis, SVector
export Pose, dimension
export BindingSite, color
export ParticleSpecies, SpeciesAndPose, bindingsites,
       graphrep, isconvex, could_contact, overlap, setcolors!, nsites, symmetrynumber

export ACCEPT, REJECT, BREAK
export polyenum, polygen
export AssemblySystem, interactionmatrix,
       nspecies, nbonds, nsites, dimension,
       bonded_colors, bonded_sites, bonded_species, isinert, species,
       numtype, posetype

export Polyform, composition, interactionmatrix, assemblysystem, bonds, bondindex, nparticles, interior_edges, exterior_edges
export PolygonParticleSpecies, UnitTriangle, UnitSquare, UnitHexagon
export PatchyParticleSpecies, PatchyDisk


include("utils.jl")
include("pose.jl")
include("bindingsite.jl")
include("particlespecies.jl")
include("assemblysystem.jl")
include("particle.jl")
include("polyform.jl")
include("enumeration.jl")

include("species/polygonparticlespecies.jl")
include("species/patchyparticlespecies.jl")

export render, polyformplot, polyformplot!

function render end
function polyformplot end
function polyformplot! end

end