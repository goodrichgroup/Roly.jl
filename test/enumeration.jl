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

    I_polymino = BindingRules(
        [1 1 1 3;
         1 2 1 4], 
         PolygonParticleSpecies(4, 1.0; labels=[1,1,1,1]))

    # Number of one-sided polyminoes (https://oeis.org/A000988) [starting from 1]
    n_polyminoes = [1, 1, 2, 7, 18, 60, 196, 704, 2500]
    n_polyminoes_cumulative = cumsum(n_polyminoes)

    I_polyiamond = BindingRules(
            [1 1 1 1;
             1 2 1 2;
             1 3 1 3], 
             PolygonParticleSpecies(3, 1.0; labels=[1,1,1]))
    
    # Number of one-sided polyiamonds (https://oeis.org/A006534)
    n_polyiamonds = [1, 1, 1, 4, 6, 19, 43, 120, 307, 866]
    n_polyiamonds_cumulative = cumsum(n_polyiamonds)

    # I_polycube = BindingRules(
    #     [1 1 1 1; 
    #      1 2 1 2;
    #      1 3 1 3;
    #      1 4 1 4;
    #      1 5 1 5;
    #      1 6 1 6], UnitCubeGeometry, [ones(Int, 24)])

    # # Number of one-sided polycubes (https://oeis.org/A006534)
    # n_polycubes = [1, 1, 2, 8, 29, 166, 1023]
    # n_polycubes_cumulative = cumsum(n_polycubes)

    nstrs_polymino_enum = [polyenum(I_polymino, maxsize=i)[1] for i in 1:length(n_polyminoes_cumulative)]
    nstrs_polyiamond_enum = [polyenum(I_polyiamond, maxsize=i)[1] for i in 1:length(n_polyiamonds_cumulative)]
    # nstrs_polycube_enum = [polyenum(I_polycube, maxsize=i)[1] for i in 1:length(n_polycubes_cumulative)]

    @test nstrs_polymino_enum == n_polyminoes_cumulative
    @test nstrs_polyiamond_enum == n_polyiamonds_cumulative
    # @test nstrs_polycube_enum == n_polycubes_cumulative

    # nstrs_polymino_gen = [length(polygen(I_polymino, maxsize=i)) for i in 1:length(n_polyminoes_cumulative)]
    # nstrs_polyiamond_gen = [length(polygen(I_polyiamond, maxsize=i)) for i in 1:length(n_polyiamonds_cumulative)]
    # nstrs_polycube_gen = [length(polygen(I_polycube, maxsize=i)) for i in 1:length(n_polycubes_cumulative)]

    # @test nstrs_polymino_gen == n_polyminoes_cumulative
    # @test nstrs_polyiamond_gen == n_polyiamonds_cumulative
    # @test nstrs_polycube_gen == n_polycubes_cumulative

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
