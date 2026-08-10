using Test
using Roly
using Rotations, StaticArrays, LinearAlgebra, Random
using NautyGraphs

@testset verbose=true "Roly" begin
    include("pose.jl")
    include("bindingsite.jl")
    include("encoding.jl")
    include("particle.jl")
    include("bindingrules.jl")
    include("polyform.jl")
    include("enumeration.jl")
    include("utils.jl")

    @testset verbose=true "species" begin
        include("species/polygonparticlespecies.jl")
        include("species/polyhedronparticlespecies.jl")
        include("species/patchyparticlespecies.jl")
    end
end;


