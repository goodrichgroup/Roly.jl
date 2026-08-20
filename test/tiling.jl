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

    @testset "unbounded rules via periodic chains" begin
        # the chain monomer repeats on its own, so one particle already proves it
        hit = canchain(chainlike)
        @test hit !== nothing && nparticles(hit) == 1
        @test !isempty(tilings(hit))              # what was found really does close
        @test isunbounded(chainlike; maxlength=1)

        # the square lattice is caught the same way, from a single particle
        @test isunbounded(squarerules)
        # a system whose structures all close is not
        @test canchain(dimerrules) === nothing
        @test !isunbounded(dimerrules)
        # rules that only fold rings back on themselves are bounded too
        @test !isunbounded(BindingRules([1 1 1 2], UnitSquare); maxlength=6)

        # the default length is derived, not guessed: a chain past it must repeat a state
        @test chainstatebound(chainlike) >= nsites(UnitSquare) ÷ 2
        @test all(chainstatebound(r) >= 1 for (r, _) in
                  ((chainlike, 0), (dimerrules, 0), (squarerules, 0)))

        # `maxlength` bounds the chain: this repeat needs two particles, not one
        turnrules = BindingRules([1 1 1 2; 1 3 1 4], UnitSquare)
        @test canchain(turnrules; maxlength=1) === nothing
        @test nparticles(canchain(turnrules; maxlength=4)) == 2
        # and the derived default reaches it without being told
        @test isunbounded(turnrules)
        # and a cell of two copies reaches the same fact from a one-particle chain
        @test nparticles(canchain(turnrules; maxlength=1, maxblock=2)) == 1
    end

    @testset "cells of several copies" begin
        # every closure above is a cell of one copy, and says so
        @test all(t.order == 1 for t in tilings(chainmono))
        @test all(t.order == 1 for t in tilings(sqmono; maxblock=1))

        # site 1 binds site 2, a quarter turn, so no translation carries a single square onto a
        # bonded copy -- the cell has to hold two of them, related by that turn, and only block
        # growth can find it
        turn = BindingRules([1 1 1 2; 1 3 1 4], UnitSquare)
        mono = first(polygen(turn; maxsize=1))
        @test isempty(tilings(mono; maxblock=1))
        blocked = tilings(mono; maxblock=2)
        @test !isempty(blocked)
        @test all(t.order == 2 for t in blocked)
        # the bond joining the two copies belongs to the cell, so it is counted there
        @test all(!isempty(t.bondtypes) for t in blocked)

        # a supercell of a tiling is still a tiling, so raising the bound only adds cells
        for k in 1:3
            ts = tilings(chainmono; maxblock=k)
            @test sort(unique(t.order for t in ts)) == collect(1:k)
            @test all(t.complete for t in ts)
            # a cell of k copies closes k bonds
            @test all(length(t.bondtypes) == t.order for t in ts)
        end
    end
end
