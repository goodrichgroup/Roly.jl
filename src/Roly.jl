module Roly

using LinearAlgebra, StaticArrays, SparseArrays, Rotations
using Base.Iterators, DataStructures
using Graphs, NautyGraphs
using ReverseSearch

export Rotation, Angle2d, RotXYZ, rotation_angle, rotation_axis, SVector
export Pose, tomatrix, dimension
export BindingSite, color
export ParticleSpecies, SpeciesAndPose, bindingsites,
       graphrep, isconvex, could_contact, overlap, setcolors!, nsites, symmetrynumber

export ACCEPT, REJECT, BREAK
export polyenum, polygen
export AssemblySystem, interactionmatrix, buildingblocks, 
       nspecies, nbonds, nsites, dimension,
       bonded_colors, bonded_sites, bonded_species, isinert

export Polyform, canonical_id, composition, compositions, interactionmatrix, isinert, assemblysystem
export PolygonParticleSpecies, UnitTriangle, UnitSquare, UnitHexagon
export MetaParticleSpecies
export cantile, lattice_vectors

include("utils.jl")
include("pose.jl")
include("bindingsite.jl")
include("particlespecies.jl")
include("assemblysystem.jl")
include("particle.jl")
include("polyform.jl")
include("enumeration.jl")
include("tiling.jl")

include("species/polygonparticlespecies.jl")
include("species/sphereparticlespecies.jl")
include("species/metaparticlespecies.jl")

export render, polyformplot, polyformplot!

function render end
function polyformplot end
function polyformplot! end

end