@testset "enumeration" begin
    I16 = BindingRules(
        [1 3 2 3;
         2 1 4 1;
         2 2 3 2;
         3 1 4 1],
        UnitTriangle)

    I137 = BindingRules(
        [1 2 2 1;
         1 3 3 1;
         1 3 3 2;
         1 3 4 3;
         2 1 3 3;
         3 2 4 2;
         3 1 4 2;
         1 2 5 1], 
        UnitTriangle)

    Icyc = BindingRules(
        [1 3 2 3;
         2 2 3 2;
         3 3 1 1;
         3 3 4 3;
         4 2 1 1],
        UnitTriangle)

    nstrs_16_enum = polyenum(I16)[1]
    nstrs_137_enum = polyenum(I137)[1]
    nstrs_cyc_enum = polyenum(Icyc)[1]

    @test nstrs_16_enum == 16
    @test nstrs_137_enum == 137
    @test nstrs_cyc_enum == 283

    # nstrs_16_gen = polygen(I16) |> length
    # nstrs_137_gen = polygen(I137) |> length
    # nstrs_cyc_gen = polygen(Icyc) |> length

    # @test nstrs_16_gen == 16
    # @test nstrs_137_gen == 137
    # @test nstrs_cyc_gen == 283

    # All four edges the same sticky stuff, so any edge binds any edge. A square is 4-fold, so
    # that is the same square lattice the old edge-1-to-edge-3 pairing described.
    I_polymino = BindingRules([1 1 1 1], PolygonParticleSpecies(4, 1.0; colors=fill(1, 4)))

    # Number of one-sided polyminoes (https://oeis.org/A000988) [starting from 1]
    n_polyminoes = [1, 1, 2, 7, 18, 60, 196, 704, 2500]
    n_polyminoes_cumulative = cumsum(n_polyminoes)

    I_polyiamond = BindingRules([1 1 1 1], PolygonParticleSpecies(3, 1.0; colors=fill(1, 3)))
    
    # Number of one-sided polyiamonds (https://oeis.org/A006534)
    n_polyiamonds = [1, 1, 1, 4, 6, 19, 43, 120, 307, 866]
    n_polyiamonds_cumulative = cumsum(n_polyiamonds)

    # 3D lattice animals. Each of these needs ring closures to work — the first polycube with
    # a ring is the 2x2 square tetracube, the first polyiamond the six-triangle hexagon — so
    # they exercise bond registration and dart pairing end to end.
    #
    # Every one is stated in colors alone: the sticky faces get one color, the rest another,
    # and the symmetry follows. Nothing is said about labels, and nothing needs to be — see
    # `_check_labelling` for why saying it would be worse than redundant.
    #
    # The planar cases come out *free* rather than one-sided: a flat assembly can be turned
    # over by a rotation about an axis in its plane, which is a rigid motion in 3D, so mirror
    # images are the same structure. That is the 3D answer, and differs from the 2D species
    # above, where reflections are not available.
    sidefaces(p) = [i for i in 1:nfaces(p) if abs(Roly.facenormal(p, i)[3]) < 1e-8]
    function rules(shp, sticky)
        colors = [i in sticky ? 1 : 2 for i in 1:nfaces(shp)]
        return BindingRules([1 first(sticky) 1 first(sticky)],
                            PolyhedronParticleSpecies(shp; colors))
    end

    # Whether a bond admits more than one registration is what separates these two groups, and
    # it is `gauge ÷ stab` of the sticky face: the turns the face has that the whole particle
    # does not. Every solid here has gauge == stab, so each bond lands in one registration and
    # the assembly is confined to its lattice.
    lattice_animals = [
        # Cubes bonded on all six faces fill space: polycubes up to rotation. A cube face is
        # 4-fold and so is the cube about it, so a neighbour can only be attached one way.
        ("polycubes (https://oeis.org/A000162)", Cube(), 1:nfaces(Cube()),
         [1, 1, 2, 8, 29, 166, 1023]),
        # Square prisms with h != a tile a plane on their four sides. The side faces are
        # rectangles, 2-fold about their normals, and the prism is 2-fold about them too.
        ("polyominoes (https://oeis.org/A000105)", Prism(4, 1.0; h=2.0),
         sidefaces(Prism(4, 1.0; h=2.0)), [1, 1, 2, 5, 12, 35, 108]),
        # Triangular prisms tile a plane, with neighbouring triangles related by a π rotation
        # rather than a translation. That the ring of six closes is the test.
        ("polyiamonds (https://oeis.org/A000577)", Prism(3, 1.0; h=2.0),
         sidefaces(Prism(3, 1.0; h=2.0)), [1, 1, 1, 3, 4, 12, 24]),
        # Hexagonal prisms tile a plane by translation: free polyhexes.
        ("polyhexes (https://oeis.org/A000228)", Prism(6, 1.0; h=2.0),
         sidefaces(Prism(6, 1.0; h=2.0)), [1, 1, 3, 7, 22, 82, 333]),
    ]

    for (name, shp, sticky, want) in lattice_animals
        sys = rules(shp, sticky)
        ps = species(sys, 1)
        # The premise of the case: one registration per bond, so it cannot leave its lattice.
        @test all(i -> Roly.nregistrations(Roly.bindingsites(ps, i), Roly.bindingsites(ps, i)) == 1,
                  sticky)
        @test [polyenum(sys; maxsize=i)[1] for i in eachindex(want)] == cumsum(want)
    end

    # The same two tilings from prisms whose side faces are *squares* rather than rectangles.
    # The face is then 4-fold about its normal where the prism is only 2-fold about it, so the
    # two solids differ in `gauge` but not in `stab` — and since registrations come from the
    # particle's symmetry rather than the face's, they must enumerate identically. Two
    # differently-shaped solids modelling the same tiling agreeing is the strongest check that
    # the twist references are being pinned at the right strength.
    for (shp, want) in [(Prism(3), [1, 2, 3, 6, 10, 22]), (Prism(6), [1, 2, 5, 12, 34])]
        sys = rules(shp, sidefaces(shp))
        b = Roly.bindingsites(species(sys, 1), first(sidefaces(shp)))
        @test (b.gauge, b.stab) == (4, 2)
        @test Roly.nregistrations(b, b) == 1
        @test [polyenum(sys; maxsize=i)[1] for i in eachindex(want)] == want
    end

    # Freeing those side faces from their orientation is how the out-of-plane lattice is asked
    # for. The square face then contributes its own 4-fold symmetry instead of the prism's
    # 2-fold, so a bond admits two registrations — the in-plane one and one standing the
    # neighbour on its side — and the assemblies leave the plane.
    for (shp, want) in [(Prism(3), [1, 3, 6, 22, 73, 357]), (Prism(6), [1, 3, 12, 81, 812])]
        sides = sidefaces(shp)
        colors = [i in sides ? 1 : 2 for i in 1:nfaces(shp)]
        ps = PolyhedronParticleSpecies(shp; colors, locking=[!(i in sides) for i in 1:nfaces(shp)])
        b = Roly.bindingsites(ps, first(sides))
        @test Roly.twistfreedom(b) == b.gauge == 4
        @test Roly.nregistrations(b, b) == 2
        sys = BindingRules([1 first(sides) 1 first(sides)], ps)
        @test [polyenum(sys; maxsize=i)[1] for i in eachindex(want)] == want
    end

    nstrs_polymino_enum = [polyenum(I_polymino, maxsize=i)[1] for i in 1:length(n_polyminoes_cumulative)]
    nstrs_polyiamond_enum = [polyenum(I_polyiamond, maxsize=i)[1] for i in 1:length(n_polyiamonds_cumulative)]

    @test nstrs_polymino_enum == n_polyminoes_cumulative
    @test nstrs_polyiamond_enum == n_polyiamonds_cumulative

    # nstrs_polymino_gen = [length(polygen(I_polymino, maxsize=i)) for i in 1:length(n_polyminoes_cumulative)]
    # nstrs_polyiamond_gen = [length(polygen(I_polyiamond, maxsize=i)) for i in 1:length(n_polyiamonds_cumulative)]

    # @test nstrs_polymino_gen == n_polyminoes_cumulative
    # @test nstrs_polyiamond_gen == n_polyiamonds_cumulative

    # Count polyforms
    @test polyenum(I16).status == Finished
    @test polyenum(I_polymino; maxsize=6).status == MaxDepthReached
    @test polyenum(I_polymino; maxstrs=50).status == MaxVerticesReached

    # Small enough to enumerate exactly, so no estimation should happen at all.
    c = countpolyforms(I137)
    @test c.exact
    @test c.n == 137
    @test c.uncertainty == 0.0
    @test !c.size_truncated

    # Capped by maxsize, but still exactly enumerated up to that cap.
    c = countpolyforms(I_polymino; maxsize=6)
    @test c.exact
    @test c.n == n_polyminoes_cumulative[6]
    @test c.size_truncated
    @test c.largest_size == 6

    # Polyominoes are unbounded, so a count without a maxsize should error.
    @test_throws ArgumentError countpolyforms(I_polymino; exact_budget=500)

    for maxsamples in (300, 3000)
        c = countpolyforms(I_polymino; maxsize=9, exact_budget=500, maxsamples, rng=Xoshiro(7))
        @test c.n >= 285 # we should get at least this many from the exact budget
        @test abs(c.n - n_polyminoes_cumulative[9]) / n_polyminoes_cumulative[9] < 0.5
    end

    # Estimation. The tolerance is deliberately loose: sweeping 200 seeds of this configuration
    # gives a median relative error of 1.7% and a worst case of 6.9%, so 25% is not a fit to the
    # seed below but headroom against the estimator's heavy tail.
    c = countpolyforms(I_polymino; maxsize=9, exact_budget=500, ntrials=5, rng=Xoshiro(1))
    @test !c.exact
    @test c.size_truncated
    @test c.largest_size == 9
    @test c.n >= 285 # we should get at least this many from the exact budget
    @test abs(c.n - n_polyminoes_cumulative[9]) / n_polyminoes_cumulative[9] < 0.25
    @test length(c.trials) == 5
    @test sum(c.trials) / 5 ≈ c.n

    # Counts per size
    counts, _, status = Roly._count_upto_budget(Icyc; maxsize=12, budget=40)
    @test status == MaxVerticesReached
    @test counts == [4, 9, 15, 22, 31]
    @test let counts_per_size = diff([0; counts]); counts_per_size[end] / counts_per_size[end-1] end ≈ 9 / 7

    # A thinning probability with p*B < 1 makes the estimate collapse towards zero while
    # reporting a small spread, so the default heuristic must stay on the p*B > 1 side.
    for (sys, budget) in ((I_polymino, 500), (Icyc, 40))
        persize = diff([0; Roly._count_upto_budget(sys; maxsize=12, budget)[1]])
        branching = persize[end] / persize[end-1]
        @test clamp(1.2 / branching, eps(), 0.95) * branching >= 1
    end

    # A budget too small to reach size 2 is raised to one that can, instead of failing.
    c = countpolyforms(I_polymino; maxsize=9, exact_budget=1, rng=Xoshiro(5))
    @test c.n >= 4
end;
