@testset "tiling" begin
    # the chain monomer is a unit cell of the infinite chain: one lattice vector closes its
    # 1-3 bond onto the translated copy
    chainlike = BindingRules([1 1 1 3], UnitSquare)
    chainmono = first(polygen(chainlike; maxsize=1))
    @test isunitcell(chainmono)
    @test length(tilelatticevectors(chainmono)) == 1
    @test all(t.complete && t.bondtypes == [1] for t in tilings(chainmono))

    # the square-lattice monomer tiles the plane with two lattice vectors; partial closures along
    # a single axis (the infinite chains) are enumerated alongside the complete tilings
    squarerules = BindingRules([1 2 1 4; 1 3 1 1], UnitSquare)
    sqmono = first(polygen(squarerules; maxsize=1))
    @test isunitcell(sqmono)
    @test length(tilelatticevectors(sqmono)) == 2
    ts = tilings(sqmono)
    @test count(t -> t.complete, ts) > 0
    @test all(length(t.bondtypes) == 2 for t in ts if t.complete)
    @test any(t -> !t.complete && length(t.bondtypes) == 1, ts)

    # a closed dimer has no open sites left to close and is not a unit cell
    dimerrules = BindingRules([1 1 2 1], UnitSquare)
    dimer = last(polygen(dimerrules; maxsize=2))
    @test !isunitcell(dimer)
    @test isempty(tilings(dimer))

    # cantile finds a unit cell among the chain structures, and none in a closing system
    @test cantile(chainlike; maxtilesize=2) !== nothing
    @test cantile(dimerrules; maxtilesize=3) === nothing
end
