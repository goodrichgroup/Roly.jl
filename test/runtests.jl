using Test
using Roly
using Rotations, StaticArrays, LinearAlgebra
using NautyGraphs

@testset verbose=true "Roly" begin
    include("pose.jl")
    include("bindingsite.jl")
    include("particle.jl")
    include("bindingrules.jl")
    include("polyform.jl")
    include("enumeration.jl")
    include("utils.jl")
    # include("tiling.jl")

    @testset verbose=true "species" begin
        include("species/polygonparticlespecies.jl")
        include("species/patchyparticlespecies.jl")
    end
end;


