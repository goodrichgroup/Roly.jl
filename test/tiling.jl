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
        @test canchain(chainlike; maxlength=1) !== nothing

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
        @test nparticles(canchain(turnrules; maxlength=1, maxorder=2)) == 1
    end

    @testset "unbounded rules via a repeating motion" begin
        # a repeat whose motion is a pure translation: the chain marches off forever
        w = growthwitness(chainlike)
        @test w !== nothing && w.period == 1
        @test isunbounded(squarerules)
        # structures that all close have no such motion
        @test growthwitness(dimerrules) === nothing
        # nor do rules that only fold back on themselves: a quarter turn repeated is a rotation,
        # and its copies stay in a bounded annulus
        @test !isunbounded(BindingRules([1 1 1 2], UnitSquare))

        # two quarter turns compose to a translation over two particles. The motion accumulates a
        # full 2pi, which has to be wrapped before it is judged to be turning at all
        turnrules = BindingRules([1 1 1 2; 1 3 1 4], UnitSquare)
        wt = growthwitness(turnrules)
        @test wt !== nothing && wt.period == 2
        @test isapprox(rem(rotation_angle(wt.generator.psi), 2pi, RoundNearest), 0; atol=1e-8)
        @test norm(wt.generator.x) ≈ 2

        # 3D: a stack of cubes is caught the same way
        cube = PolyhedronParticleSpecies(Cube(); colors=fill(1, 6))
        @test isunbounded(BindingRules([1 1 1 1], cube))
    end

    @testset "cells of several copies" begin
        # every closure above is a cell of one copy, and says so
        @test all(t.order == 1 for t in tilings(chainmono))
        @test all(t.order == 1 for t in tilings(sqmono; maxorder=1))

        # site 1 binds site 2, a quarter turn, so no translation carries a single square onto a
        # bonded copy -- the cell has to hold two of them, related by that turn, and only
        # growing the cell as a meta-polyform can find it
        turn = BindingRules([1 1 1 2; 1 3 1 4], UnitSquare)
        mono = first(polygen(turn; maxsize=1))
        @test isempty(tilings(mono; maxorder=1))
        turned = tilings(mono; maxorder=2)
        @test !isempty(turned)
        @test all(t.order == 2 for t in turned)
        # the bond joining the two copies belongs to the cell, so it is counted there
        @test all(!isempty(t.bondtypes) for t in turned)

        # a supercell of a tiling is still a tiling, so raising the bound only adds cells
        for k in 1:3
            ts = tilings(chainmono; maxorder=k)
            @test sort(unique(t.order for t in ts)) == collect(1:k)
            @test all(t.complete for t in ts)
            # a cell of k copies closes k bonds
            @test all(length(t.bondtypes) == t.order for t in ts)
        end
    end

    @testset "3D" begin
        using Roly: nfaces, facenormal, twistfreedom

        # every face binds every face, so the cube monomer is a unit cell of space itself
        cubic = BindingRules([1 1 1 1], PolyhedronParticleSpecies(Cube(); colors=fill(1, 6)))
        cubemono = first(polygen(cubic; maxsize=1))
        @test isunitcell(cubemono)
        @test length(tilelatticevectors(cubemono)) == 3
        cts = tilings(cubemono)
        # a complete tiling closes all three axes, each through the system's one bond type
        complete = filter(t -> t.complete, cts)
        @test !isempty(complete)
        @test all(t -> length(t.vectors) == 3 && t.bondtypes == [1, 1, 1], complete)
        # the partial closures come back alongside: a column with one vector, a sheet with two
        @test sort(unique(length(t.vectors) for t in cts)) == [1, 2, 3]
        @test all(!t.complete for t in cts if length(t.vectors) < 3)
        @test all(t.order == 1 for t in cts)

        # four inert sides leave one axis to close, and one vector closes it
        column = BindingRules([1 1 1 1], PolyhedronParticleSpecies(Cube(); colors=[1, 2, 2, 2, 2, 1]))
        colmono = first(polygen(column; maxsize=1))
        @test length(opensites(colmono)) == 2
        @test all(t -> t.complete && length(t.vectors) == 1, tilings(colmono))

        # distinct colors do not stop a particle tiling. `_canonical_faces` starts two faces that
        # face each other at corresponding corners, so a bond between them turns the neighbour not
        # at all, and a six-color cube tiles the cubic lattice with no twists asked for
        pure(rules) = all(polygen(rules; maxsize=2)) do p
            nparticles(p) < 2 && return true
            g = p.particles[2].pose * inv(p.particles[1].pose)
            return isapprox(rem(rotation_angle(g.psi), 2pi, RoundNearest), 0; atol=1e-8)
        end
        opposites = [1 1 1 6; 1 2 1 5; 1 3 1 4]
        keyedcube = BindingRules(opposites, PolyhedronParticleSpecies(Cube()))
        @test pure(keyedcube)
        kcube = filter(t -> t.complete, tilings(first(polygen(keyedcube; maxsize=1))))
        @test !isempty(kcube)
        # the lattice is the three unit axes, and each bond type is spent once -- where the
        # one-color cube spends its single type three times
        @test all(t -> sort(norm.(t.vectors)) ≈ [1, 1, 1], kcube)
        @test all(t -> sort(t.bondtypes) == [1, 2, 3], kcube)
        # and squaring the faces up costs nothing, since the turns it takes are whole steps of a
        # square face's own symmetry: a cube colored alike keeps its full group
        @test symmetrynumber(PolyhedronParticleSpecies(Cube(); colors=fill(1, 6))) == 24

        # a solid whose opposite faces are anti-aligned mates the other way. Half a step is not
        # something a threefold face can absorb, but a further half turn is, so bonded octahedra
        # meet turned by 180 degrees: square, and still free, but no translation carries one onto
        # the other and there is no cell of one copy
        octopposite(i) = findfirst(j -> isapprox(dot(facenormal(Octahedron(), i),
                                                     facenormal(Octahedron(), j)), -1; atol=1e-8),
                                   1:nfaces(Octahedron()))
        octpairs = unique([minmax(i, octopposite(i)) for i in 1:nfaces(Octahedron())])
        octa = BindingRules(reduce(vcat, [[1 a 1 b] for (a, b) in octpairs]),
                            PolyhedronParticleSpecies(Octahedron()))
        @test all(polygen(octa; maxsize=2)) do p
            nparticles(p) < 2 && return true
            g = p.particles[2].pose * inv(p.particles[1].pose)
            return isapprox(rotation_angle(g.psi), pi; atol=1e-8)
        end
        @test symmetrynumber(PolyhedronParticleSpecies(Octahedron(); colors=fill(1, 8))) == 24
        @test isempty(tilings(first(polygen(octa; maxsize=1))))

        # unkeyed faces reach a lattice without any of this. A non-locking site admits every twist
        # its own symmetry allows, four for a square face, so the translation is among them
        free = PolyhedronParticleSpecies(Cube(); locking=false)
        @test twistfreedom(bindingsite(free, 1), bindingsite(free, 6)) == 4
        @test symmetrynumber(PolyhedronParticleSpecies(Cube(); colors=fill(1, 6), locking=false)) == 24
        freerules = BindingRules(opposites, free)
        @test count(p -> nparticles(p) == 2, polygen(freerules; maxsize=2)) == 12
        @test isunitcell(first(polygen(freerules; maxsize=1)))

        # and patches say the same thing as faces, with the twist named directly rather than taken
        # from the corner order
        axes = [SVector(-0.5, 0.0, 0.0), SVector(0.0, -0.5, 0.0), SVector(0.0, 0.0, -0.5),
                SVector(0.0, 0.0, 0.5), SVector(0.0, 0.5, 0.0), SVector(0.5, 0.0, 0.0)]
        keyed = PatchyParticleSpecies(NautyDiGraph(6), 0.5, axes, [0.0, 0.0, 0.0, pi, 0.0, 0.0]; colors=1:6)
        keyedpatches = BindingRules(opposites, keyed)
        @test allunique(color(bindingsite(keyed, i)) for i in 1:6)
        @test pure(keyedpatches)
        kpatch = filter(t -> t.complete, tilings(first(polygen(keyedpatches; maxsize=1))))
        @test !isempty(kpatch)
        @test all(t -> sort(norm.(t.vectors)) ≈ [1, 1, 1], kpatch)
        @test all(t -> sort(t.bondtypes) == [1, 2, 3], kpatch)

        # supercells behave as they do in the plane
        @test sort(unique(t.order for t in tilings(cubemono; maxorder=2))) == [1, 2]
        @test cantile(cubic; maxtilesize=1) !== nothing
        @test cantile(octa; maxtilesize=2) === nothing
    end
end
