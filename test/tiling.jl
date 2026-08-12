@testset "tiling" begin
    # Known-broken contract for the periodic-witness machinery, pending its port from git
    # history (dcc134c..0b53e61, removed in afb02b1). These pin the intended API; an unexpected
    # pass means the port landed and the test should be promoted to @test.

    # a single chain particle is a unit cell of the infinite chain: one lattice vector closes
    # its 1-3 bond onto the translated copy
    chainlike = BindingRules([1 1 1 3], UnitSquare)
    chainmono = first(polygen(chainlike; maxsize=1))
    @test_broken isunitcell(chainmono)
    @test_broken length(tilelatticevectors(chainmono)) == 1

    # the square-lattice monomer tiles the plane with two lattice vectors
    squarerules = BindingRules([1 2 1 4; 1 3 1 1], UnitSquare)
    sqmono = first(polygen(squarerules; maxsize=1))
    @test_broken isunitcell(sqmono)
    @test_broken length(tilelatticevectors(sqmono)) == 2

    # a closed dimer has no open sites left to close and is not a unit cell
    dimerrules = BindingRules([1 1 2 1], UnitSquare)
    dimer = last(polygen(dimerrules; maxsize=2))
    @test_broken !isunitcell(dimer)

    # cantile searches super-structures of a seed for a unit cell
    @test_broken cantile(chainmono; maxtilesize=4)
end
