using Roly: Particle, graphvertices, leading_vertex, species_index

@testset "particle" begin
    sys = AssemblySystem([1 1 1 3; 1 2 1 4], UnitSquare)

    part = Particle(sys, 1; leading_vertex=1)
    @test leading_vertex(part) == 1
    @test species_index(part) == 1

    part2 = Particle(sys, 1; leading_vertex=1)
    gvs = graphvertices(part2, sys)
    @test first(gvs) == 1

    part_shifted = Particle(sys, 1; leading_vertex=length(gvs) + 1)
    gvs2 = graphvertices(part_shifted, sys)
    @test first(gvs2) == length(gvs) + 1
    @test length(gvs2) == length(gvs)

    @test nsites(part, sys) == 4

    site1 = bindingsites(part, sys, 1)
    @test first(site1.vertices) == 1
    site1_shifted = bindingsites(part_shifted, sys, 1)
    @test first(site1_shifted.vertices) == length(gvs) + 1

    @test length(collect(bindingsites(part, sys))) == nsites(part, sys)

    v = SVector(3.0, 0.0)
    part_far = part + v
    @test part_far.pose.x ≈ v

    rot = Angle2d(π/2)
    @test (rot * part).pose.psi ≈ rot * part.pose.psi
    @test (part * rot).pose.psi ≈ part.pose.psi * rot

    @test overlap(part, part2, sys)
    @test !overlap(part, part_far, sys)
    @test !could_contact(part, part_far, sys)

    io = IOBuffer()
    show(io, part)
    @test contains(String(take!(io)), "Particle")
    show(io, MIME"text/plain"(), part)
    @test contains(String(take!(io)), "species_index")
end
