using Test
using Roly
using Rotations, StaticArrays, LinearAlgebra, Random
using NautyGraphs
using JLD2
using CairoMakie

include("testutils.jl")

@testset verbose=true "Roly" begin
    include("pose.jl")
    include("bindingsite.jl")
    include("encoding.jl")
    include("particle.jl")
    include("bindingrules.jl")
    include("polyform.jl")
    include("enumeration.jl")
    include("utils.jl")
    include("conventions.jl")
    include("plotting.jl")
    include("custom_species.jl")
    include("symmetry.jl")
    include("io.jl")

    @testset verbose=true "species" begin
        include("species/polygonparticlespecies.jl")
        include("species/polyhedronparticlespecies.jl")
        include("species/patchyparticlespecies.jl")
        include("species/metaparticlespecies.jl")
    end
end;


