using Test
using Roly
using Rotations, StaticArrays

@testset verbose=true "Roly" begin
    # include("geometry.jl")
    # include("assembly_system.jl")
    # include("polyform.jl")
    # include("utils.jl")
    # include("enumeration.jl")
    include("pose.jl")
    include("bindingsite.jl")
end;