module Roly

using LinearAlgebra, StaticArrays, SparseArrays, Rotations
using Base.Iterators, DataStructures
using Graphs, NautyGraphs
using ReverseSearch

export ACCEPT, REJECT, BREAK
export polyenum, polygen
export AssemblySystem, interactionmatrix, buildingblocks, 
       nspecies, nbonds, nsites, dimension,
       bondlist, bonded_sites, bonded_species, isinert


export Polyform, canonical_id, composition, compositions, interactionmatrix, isinert
export PolygonGeometry, UnitTriangleGeometry, UnitSquareGeometry, UnitPentagonGeometry, UnitHexagonGeometry, UnitCubeGeometry

include("utils.jl")
include("geometry_utils.jl")
# include("geometry.jl")
# include("buildingblock.jl")
include("particlespecies.jl")
include("assemblysystem.jl")
include("particle.jl")
include("polyform_new.jl")
# include("polyform.jl")
# include("concatenation.jl")
# include("enumeration.jl")
end