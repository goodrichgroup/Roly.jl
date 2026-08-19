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


    I_polymino = BindingRules([1 1 1 1], PolygonParticleSpecies(4, 1.0; colors=fill(1, 4)))

    # Number of one-sided polyminoes (https://oeis.org/A000988) [starting from 1]
    n_polyminoes = [1, 1, 2, 7, 18, 60, 196, 704, 2500]
    n_polyminoes_cumulative = cumsum(n_polyminoes)

    I_polyiamond = BindingRules([1 1 1 1], PolygonParticleSpecies(3, 1.0; colors=fill(1, 3)))
    
    # Number of one-sided polyiamonds (https://oeis.org/A006534)
    n_polyiamonds = [1, 1, 1, 4, 6, 19, 43, 120, 307, 866]
    n_polyiamonds_cumulative = cumsum(n_polyiamonds)

    # 3D lattice animals
    sidefaces(p) = [i for i in 1:nfaces(p) if abs(Roly.facenormal(p, i)[3]) < 1e-8]
    function rules(shp, sticky)
        colors = [i in sticky ? 1 : 2 for i in 1:nfaces(shp)]
        return BindingRules([1 first(sticky) 1 first(sticky)],
                            PolyhedronParticleSpecies(shp; colors))
    end

    lattice_animals = [
        ("polycubes (https://oeis.org/A000162)", Cube(), 1:nfaces(Cube()),
         [1, 1, 2, 8, 29, 166, 1023]),
        ("polyominoes (https://oeis.org/A000105)", Prism(4, 1.0; h=2.0),
         sidefaces(Prism(4, 1.0; h=2.0)), [1, 1, 2, 5, 12, 35, 108]),
        ("polyiamonds (https://oeis.org/A000577)", Prism(3, 1.0; h=2.0),
         sidefaces(Prism(3, 1.0; h=2.0)), [1, 1, 1, 3, 4, 12, 24]),
        ("polyhexes (https://oeis.org/A000228)", Prism(6, 1.0; h=2.0),
         sidefaces(Prism(6, 1.0; h=2.0)), [1, 1, 3, 7, 22, 82, 333]),
    ]

    for (name, shp, sticky, want) in lattice_animals
        sys = rules(shp, sticky)
        ps = species(sys, 1)
        # one phase per bond, so it cannot leave the lattice.
        @test all(i -> Roly.nphases(Roly.bindingsites(ps, i), Roly.bindingsites(ps, i)) == 1,
                  sticky)
        @test [polyenum(sys; maxsize=i)[1] for i in eachindex(want)] == cumsum(want)
    end

    # The same two tilings from prisms whose side faces are *squares* rather than rectangles.
    # The face is then 4-fold about its normal where the prism is only 2-fold about it, so the
    # two solids differ in `gauge` but not in `stab`.
    for (shp, want) in [(Prism(3), [1, 2, 3, 6, 10, 22]), (Prism(6), [1, 2, 5, 12, 34])]
        sys = rules(shp, sidefaces(shp))
        b = Roly.bindingsites(species(sys, 1), first(sidefaces(shp)))
        @test (b.gauge, b.stab) == (4, 2)
        @test Roly.nphases(b, b) == 1
        @test [polyenum(sys; maxsize=i)[1] for i in eachindex(want)] == want
    end

    # setting locking=false, allows out of plane binding
    for (shp, want) in [(Prism(3), [1, 3, 6, 22, 73, 357]), (Prism(6), [1, 3, 12, 81, 812])]
        sides = sidefaces(shp)
        colors = [i in sides ? 1 : 2 for i in 1:nfaces(shp)]
        ps = PolyhedronParticleSpecies(shp; colors, locking=[!(i in sides) for i in 1:nfaces(shp)])
        b = Roly.bindingsites(ps, first(sides))
        @test Roly.twistfreedom(b) == b.gauge == 4
        @test Roly.nphases(b, b) == 2
        sys = BindingRules([1 first(sides) 1 first(sides)], ps)
        @test [polyenum(sys; maxsize=i)[1] for i in eachindex(want)] == want
    end

    # 2D hexagons (one-sided) vs 3D hexagons (free)
    hexagons = BindingRules([1 1 1 1], PolygonParticleSpecies(6, 1.0; colors=fill(1, 6)))
    n_polyhexes_onesided = [1, 1, 3, 10, 33, 147]     # https://oeis.org/A006535
    @test [polyenum(hexagons; maxsize=i)[1] for i in eachindex(n_polyhexes_onesided)] ==
          cumsum(n_polyhexes_onesided)

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

    # Estimation. The tolerance is deliberately loose
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

    # replacable species test

    # `k` copies of one species, every bonding site compatible across every pair of them, so
    # the shapes are exactly the single-species shapes and the only new freedom is which
    # species sits at each position. 

    # Each species below has its two bonding sites facing opposite ways, so assemblies are
    # straight chains: one shape per size, and the only automorphism that moves particles is
    # reversal. Whether reversal is available is decided by `alike` — with the two sites the
    # same color it maps the chain onto itself, with different colors it would have to swap
    # two colors and is not an automorphism at all. So
    #
    #   distinct sites:  k^n            (trivial action on particles, every assignment distinct)
    #   alike sites:     (k^n + k^⌈n/2⌉)/2   (Burnside over {id, reversal}, which has ⌈n/2⌉ cycles)
    #
    # Four geometries, so the identity is checked against both graph encodings and in both
    # dimensions, and across species symmetry numbers 1, 2, 3, 4, 6 and 8 — including the
    # chiral C₃ and C₄, which no other enumeration test reaches.
    opposite(p, i) = findfirst(j -> isapprox(dot(Roly.facenormal(p, i), Roly.facenormal(p, j)), -1;
                                             atol=1e-8), 1:nfaces(p))
    # Each entry builds a species whose two bonding sites carry colors `a` and `b`, and states
    # the symmetry numbers it should have with those two distinct and with them alike.
    function prismchain(a, b)
        p = Prism(3, 1.0; h=2.0)
        caps = [i for i in 1:nfaces(p) if abs(Roly.facenormal(p, i)[3]) > 1e-8]
        return PolyhedronParticleSpecies(
            p; colors=[i == caps[1] ? a : i == caps[2] ? b : 3 for i in 1:nfaces(p)]
        )
    end
    function cubechain(a, b)
        p = Cube()
        return PolyhedronParticleSpecies(
            p; colors=[i == 1 ? a : i == opposite(p, 1) ? b : 3 for i in 1:nfaces(p)]
        )
    end
    chains = [
        ("disk", (a, b) -> PatchyDisk([0.0, π]; colors=[a, b]), (1, 2)),
        ("square", (a, b) -> PolygonParticleSpecies(4, 1.0; colors=[a, 3, b, 3]), (1, 2)),
        ("prism", prismchain, (3, 6)),
        ("cube", cubechain, (4, 8)),
    ]

    N = 5
    for (name, build, symnums) in chains
        chainspecies(alike) = alike ? build(1, 1) : build(1, 2)
        for (alike, want_sym) in zip((false, true), symnums)
            ps = chainspecies(alike)
            @test symmetrynumber(ps) == want_sym
            # The two bonding sites, found by color so each geometry can index its own faces.
            s1 = findfirst(i -> color(bindingsites(ps, i)) == 1, 1:nsites(ps))
            s2 = alike ? s1 : findfirst(i -> color(bindingsites(ps, i)) == 2, 1:nsites(ps))

            for k in 1:3
                bonds = reduce(vcat, [[s s1 t s2] for s in 1:k for t in 1:k])
                sys = BindingRules(bonds, [chainspecies(alike) for _ in 1:k])
                @test nspecies(sys) == k
                want = alike ? [(k^n + k^cld(n, 2)) ÷ 2 for n in 1:N] : [k^n for n in 1:N]
                @test [polyenum(sys; maxsize=n)[1] for n in 1:N] == cumsum(want)
            end
        end
    end


    # Reference enumerations
   
    # `data/enumeration_reference.txt` holds 2725 assembly systems together with the size counts
    # each produced under commit 1bfe470, the last revision of the old implementation.
    # One system per line, `geometry|bond rows|expected counts`, the bond rows comma-separated and
    # each row `species1 site1 species2 site2`.
    referencespecies(s) = s === :square   ? UnitSquare :
                          s === :triangle ? UnitTriangle :
                          s === :hexagon  ? UnitHexagon :
                          error("unknown geometry $s")

    function sizecounts(sys)
        counts = Int[]
        polyenum(sys) do _, n
            while length(counts) < n
                push!(counts, 0)
            end
            counts[n] += 1
            return ACCEPT
        end
        return counts
    end

    mismatches = String[]
    nchecked = 0
    for (i, line) in enumerate(eachline(joinpath(@__DIR__, "data", "enumeration_reference.txt")))
        isempty(line) && continue
        geom, table, want = split(line, '|')
        bonds = reduce(vcat, (reshape(parse.(Int, split(r, ' ')), 1, 4)
                              for r in split(table, ',')))
        expected = parse.(Int, split(want, ','))
        counts = sizecounts(BindingRules(bonds, referencespecies(Symbol(geom))))
        counts == expected || push!(mismatches, "system $i: got $counts, want $expected")
        nchecked += 1
    end

    @test mismatches == String[]
    # Guards against a truncated or unreadable data file passing silently.
    @test nchecked == 2725
end;
