module Roly

using LinearAlgebra, StaticArrays, SparseArrays, Rotations
using Base.Iterators, DataStructures
using Graphs, NautyGraphs
using ReverseSearch

export ACCEPT, REJECT, BREAK
export polyenum, polygen
export AssemblySystem, BuildingBlock, Polyform, canonical_id, composition, compositions, intmat, isinert
export PolygonGeometry, UnitTriangleGeometry, UnitSquareGeometry, UnitPentagonGeometry, UnitHexagonGeometry, UnitCubeGeometry

include("utils.jl")
include("geometry_utils.jl")
include("geometry.jl")
include("buildingblock.jl")
# include("polyform.jl")
include("assemblysystem.jl")
# include("concatenation.jl")
# include("enumeration.jl")
end