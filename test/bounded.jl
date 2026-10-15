@testset "bounded" begin
    chainlike = BindingRules([1 1 1 3], UnitSquare)
    squarerules = BindingRules([1 2 1 4; 1 3 1 1], UnitSquare)
    dimerrules = BindingRules([1 1 2 1], UnitSquare)
    ring = BindingRules([1 1 1 2], UnitSquare)
    turnrules = BindingRules([1 1 1 2; 1 3 1 4], UnitSquare)

    # a chain that marches off under a pure translation, from one particle
    @test !isbounded(chainlike)
    @test !isbounded(squarerules)
    # a system whose structures all close
    @test isbounded(dimerrules)
    # rules that only fold rings back on themselves: a quarter turn repeated is a rotation, and its
    # copies stay in a bounded annulus
    @test isbounded(ring)
    # two quarter turns compose to a translation over two particles
    @test !isbounded(turnrules)
    # a cell holding two species, which no single monomer's copies make
    @test !isbounded(BindingRules([1 4 2 2; 2 1 1 3], UnitSquare))

    # the motion behind the answer
    w = Roly._repeatingchain(chainlike)
    @test w !== nothing && w.period == 1
    @test Roly._repeatingchain(dimerrules) === nothing
    @test Roly._repeatingchain(ring) === nothing
    # the motion accumulates a full 2pi over the two turns, which has to be wrapped before it is
    # judged to be turning at all
    wt = Roly._repeatingchain(turnrules)
    @test wt !== nothing && wt.period == 2
    @test isapprox(rem(rotation_angle(wt.generator.psi), 2pi, RoundNearest), 0; atol=1e-8)
    @test norm(wt.generator.x) ≈ 2

    # the default length is derived, not guessed: a chain past it must repeat a state
    @test Roly._maxchainlength(chainlike) >= nsites(UnitSquare) ÷ 2
    @test all(Roly._maxchainlength(r) >= 1 for r in (chainlike, dimerrules, squarerules))
    # `maxlength` caps the chains tried: this repeat needs two particles, and a budget of one does
    # not reach it, so the answer below the default is a budget running out and not a proof
    @test isbounded(turnrules; maxlength=1)
    @test !isbounded(turnrules; maxlength=4)

    # 3D: a stack of cubes is caught the same way
    cube = PolyhedronParticleSpecies(Cube(); colors=fill(1, 6))
    @test !isbounded(BindingRules([1 1 1 1], cube))

    # and so is growth that never repeats periodically. Regular tetrahedra glued face to face form
    # a helix whose twist per tetrahedron is an irrational fraction of a turn, so no stretch of it
    # ever closes under a translation. The repeat is a screw, and its climb is what escapes
    tetra = PolyhedronParticleSpecies(Tetrahedron(); colors=fill(1, 4))
    helix = BindingRules([1 1 1 1], tetra)
    @test !isbounded(helix)
    wh = Roly._repeatingchain(helix)
    @test isapprox(rotation_angle(wh.generator.psi), acos(-2 / 3); atol=1e-6)
    @test abs(dot(wh.generator.x, rotation_axis(wh.generator.psi))) > 0.1
end
