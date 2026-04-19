module Roly

using LinearAlgebra, StaticArrays, SparseArrays, Rotations
using Base.Iterators, DataStructures
using Graphs, NautyGraphs
using ReverseSearch
using CairoMakie

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


export Polyform, canonical_id, composition, compositions, interactionmatrix, isinert
export PolygonParticleSpecies, UnitTriangle, UnitSquare, UnitHexagon

include("utils.jl")
include("pose.jl")
include("bindingsite.jl")
include("particlespecies.jl")
include("assemblysystem.jl")
include("particle.jl")
include("polyform.jl")
include("enumeration.jl")

include("species/polygonparticlespecies.jl")
include("species/sphereparticlespecies.jl")
end