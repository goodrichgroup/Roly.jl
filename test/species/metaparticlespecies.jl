@testset "MetaParticleSpecies" begin
    using Roly: setcolors!, opensites, overlap, polyform
    using Graphs: ne

    # one species chaining through opposite sites, so a dimer keeps one open end of each kind
    chainlike = BindingRules([1 1 1 3], UnitSquare)
    chainstrs = polygen(chainlike; maxsize=4)
    dimer = chainstrs[findfirst(s -> nparticles(s) == 2, chainstrs)]
    open = opensites(dimer)
    @test length(open) == 2

    mp = MetaParticleSpecies(dimer)
    @test dimension(mp) == 2
    @test nsites(mp) == length(open)
    @test polyform(mp) == dimer
    @test [color(bindingsite(mp, i)) for i in 1:nsites(mp)] == [color(s) for s in open]
    @test [bindingsite(mp, i).pose for i in 1:nsites(mp)] == [s.pose for s in open]
    @test bounding_radius(mp) >= maximum(norm(p.pose.x) for p in dimer.particles)
    # `_ndistincttwists` needs the twist freedoms, which are 1 in 2D
    @test stabilizerorders(mp) == ones(Int, nsites(mp))
    @test _ndistincttwists(bindingsite(mp, 1), bindingsite(mp, 2)) == 1

    @testset "the encoding is the polyform's own graph" begin
        # bound sites stay as interior vertices, so the graph is larger than the site count
        @test nv(graphrep(mp)) == nv(graphrep(dimer))
        @test nv(graphrep(mp)) > nsites(mp)
        @test ne(graphrep(mp)) == ne(graphrep(dimer))

        # a symmetric polyform: two squares joined site1-site1, swapping them is an automorphism
        selfrules = BindingRules([1 1 1 1; 1 2 1 2], UnitSquare)
        selfstrs = polygen(selfrules; maxsize=2)
        sym = selfstrs[findfirst(s -> nparticles(s) == 2, selfstrs)]
        @test symmetrynumber(sym) == 2
        nopen = length(opensites(sym))
        @test nopen == 2

        # the species inherits the polyform's symmetry exactly -- an encoding over the sites alone
        # could not have known it
        @test symmetrynumber(MetaParticleSpecies(sym)) == symmetrynumber(sym)
        @test symmetrynumber(MetaParticleSpecies(sym; colors=fill(9, nopen))) == 2
        # colours a bond can tell apart break it
        @test symmetrynumber(MetaParticleSpecies(sym; colors=1:nopen)) == 1
    end

    @testset "a meta-assembly is the assembly its particles would have formed" begin
        metasys = BindingRules([1 1 1 2], MetaParticleSpecies(dimer))
        mstrs = polygen(metasys; maxsize=2)
        pair = mstrs[findfirst(s -> nparticles(s) == 2, mstrs)]
        quad = chainstrs[findfirst(s -> nparticles(s) == 4, chainstrs)]

        gm, gp = graphrep(pair), graphrep(quad)
        @test nv(gm) == nv(gp) && ne(gm) == ne(gp)
        # identical up to the labeling, which differs by construction: open-site vertices carry
        # their colors and every species' labels are shifted into a global range
        flatten(g) = (h = copy(g); setlabels!(h, fill(Cint(1), nv(h))); nauty(h; canonize=true); h)
        @test flatten(gm) == flatten(gp)
    end

    @testset "choosing which sites to expose" begin
        # every unbound site is offered, bound ones never are
        @test length(exposedsites(ParticleSite, dimer)) == 2 * nsites(UnitSquare) - 2
        @test all(((p, k),) -> 1 <= p <= nparticles(dimer), exposedsites(ParticleSite, dimer))
        # the open ones are the non-inert ones, and are what the default exposes, in this order
        @test length(opensites(ParticleSite, dimer)) == length(open)
        @test opensites(dimer) == filter(s -> !Roly.isinert(chainlike, color(s)), exposedsites(dimer))
        @test [color(bindingsite(mp, i)) for i in 1:nsites(mp)] == [color(s) for s in opensites(dimer)]
        # the two views agree site for site
        @test exposedsites(dimer) ==
              [bindingsite(dimer.particles[p], chainlike, k) for (p, k) in exposedsites(ParticleSite, dimer)]

        # a site that is inert inside the polyform becomes usable by being named and given a
        # colour the new rules speak
        inert = setdiff(exposedsites(ParticleSite, dimer), opensites(ParticleSite, dimer))
        @test !isempty(inert)
        activated = MetaParticleSpecies(dimer, inert[1:2]; colors=[1, 2])
        @test nsites(activated) == 2
        @test [color(bindingsite(activated, i)) for i in 1:2] == [1, 2]
        # and those sites really do bind under rules written over the new colors
        @test [polyenum(BindingRules([1 1 1 2], activated); maxsize=k)[1] for k in 1:3] ==
              [1, 2, 3]

        # exposure order is the caller's
        pair = opensites(ParticleSite, dimer)
        @test [color(bindingsite(MetaParticleSpecies(dimer, reverse(pair)), i)) for i in 1:2] ==
              reverse([color(s) for s in opensites(dimer)])

        bound = [ParticleSite(p, k) for p in 1:nparticles(dimer) for k in 1:nsites(UnitSquare)
                 if ParticleSite(p, k) ∉ Set(exposedsites(ParticleSite, dimer))]
        @test length(bound) == 2
        @test_throws ArgumentError MetaParticleSpecies(dimer, bound[1:1])
        @test_throws ArgumentError MetaParticleSpecies(dimer, [(1, 99)])
        @test_throws ArgumentError MetaParticleSpecies(dimer, [(99, 1)])
        @test_throws ArgumentError MetaParticleSpecies(dimer, [pair[1], pair[1]])
        @test_throws ArgumentError MetaParticleSpecies(dimer, Tuple{Int,Int}[])
    end

    @testset "exposure decides the symmetry" begin
        # the symmetric stack again: exposing both equivalent ends keeps the swap, exposing one
        # of them cannot
        selfrules = BindingRules([1 1 1 1; 1 2 1 2], UnitSquare)
        sym = let strs = polygen(selfrules; maxsize=2)
            strs[findfirst(s -> nparticles(s) == 2, strs)]
        end
        live = opensites(ParticleSite, sym)
        @test length(live) == 2
        @test symmetrynumber(MetaParticleSpecies(sym, live)) == 2
        @test symmetrynumber(MetaParticleSpecies(sym, live[1:1])) == 1
    end

    @testset "keywords and rejections" begin
        recolored = MetaParticleSpecies(dimer; colors=[4, 9])
        @test [color(bindingsite(recolored, i)) for i in 1:2] == [4, 9]
        @test_throws DimensionMismatch MetaParticleSpecies(dimer; colors=[1])

        # a saturated polyform exposes nothing: the 1-2-3 chain's ends are inert
        chain = BindingRules([1 1 2 3; 2 1 3 3], UnitSquare)
        trimer = let strs = polygen(chain; maxsize=3)
            strs[findfirst(s -> nparticles(s) == 3, strs)]
        end
        @test isempty(opensites(trimer))
        @test_throws ArgumentError MetaParticleSpecies(trimer)

    end

    @testset "3D: symmetry and twist freedom come from the polyform" begin
        # every face alike, so a cube is 24-fold and each face 4-fold about its normal
        cube = PolyhedronParticleSpecies(Cube(); colors=fill(1, 6))
        @test symmetrynumber(cube) == 24
        @test stabilizerorders(cube) == fill(4, 6)

        stack = BindingRules([1 1 1 1], cube)
        cubedimer = let strs = polygen(stack; maxsize=2)
            strs[findfirst(s -> nparticles(s) == 2, strs)]
        end
        # 4 turns about the stacking axis, times swapping the two cubes
        @test symmetrynumber(cubedimer) == 8

        cmp = MetaParticleSpecies(cubedimer)
        @test dimension(cmp) == 3
        @test symmetrynumber(cmp) == symmetrynumber(cubedimer)
        @test nsites(cmp) == 10
        # a dart-encoded face spans four vertices, which stay contiguous and in order
        @test all(length(bindingsite(cmp, i).vertices) == 4 for i in 1:nsites(cmp))
        @test nv(graphrep(cmp)) == nv(graphrep(cubedimer))

        # the two faces opposite the bond keep the whole 4-fold twist, since turning about them
        # carries the stack onto itself; the eight side faces keep none of it
        stabs = stabilizerorders(cmp)
        @test count(==(4), stabs) == 2
        @test count(==(1), stabs) == 8
        @test all(bindingsite(cmp, i).sitesym == 4 for i in 1:nsites(cmp))

        # a cube with distinctly coloured faces has no symmetry to inherit
        plaindimer = let strs = polygen(BindingRules([1 1 1 2], UnitCube); maxsize=2)
            strs[findfirst(s -> nparticles(s) == 2, strs)]
        end
        @test stabilizerorders(MetaParticleSpecies(plaindimer)) == [1, 1]

        # and blocks assemble out of blocks in 3D too
        metasys = BindingRules([1 1 1 1], cmp)
        @test polyenum(metasys; maxsize=2)[1] >= 2
    end

    @testset "setcolors! relabels the sites and leaves the interior alone" begin
        ps = MetaParticleSpecies(dimer)
        sitevs = Set(v for i in 1:nsites(ps) for v in bindingsite(ps, i).vertices)
        interior = [v for v in 1:nv(graphrep(ps)) if v ∉ sitevs]
        before = labels(graphrep(ps))[interior]

        setcolors!(ps, [7, 8])
        @test [color(bindingsite(ps, i)) for i in 1:2] == [7, 8]
        @test labels(graphrep(ps))[interior] == before
        # site labels sit above every interior one, so a color can never alias interior structure
        @test all(labels(graphrep(ps))[v] > maximum(before) for v in sitevs)
        @test stabilizerorders(ps) == [1, 1]
        @test_throws ArgumentError setcolors!(ps, [1, 2, 3])

        # equal colors on the two ends make them interchangeable only if the polyform agrees;
        # this chain dimer is not symmetric, so it stays rigid
        setcolors!(ps, [5, 5])
        @test symmetrynumber(ps) == symmetrynumber(dimer)
    end

    @testset "overlap delegates to the constituents" begin
        P = Roly.posetype(mp)
        @test overlap(mp => one(P), mp => one(P))
        @test !overlap(mp => one(P), mp => Pose(SVector(50.0, 0.0), Angle2d(0.0)))
        @test !overlap(mp => one(P), mp => Pose(SVector(2.0, 0.0), Angle2d(0.0)))
        @test overlap(mp => one(P), mp => Pose(SVector(0.5, 0.0), Angle2d(0.0)))
    end

    @testset "assembling blocks out of blocks" begin
        metasys = BindingRules([1 1 1 2], MetaParticleSpecies(dimer))
        # the two ends are distinguishable, so blocks chain in exactly one way per length
        @test [polyenum(metasys; maxsize=k)[1] for k in 1:4] == [1, 2, 3, 4]
        blocks = polygen(metasys; maxsize=3)
        @test sort(nparticles.(blocks)) == [1, 2, 3]
        triple = blocks[findfirst(s -> nparticles(s) == 3, blocks)]
        @test nparticles(triple) * nparticles(dimer) == 6
    end

    @testset "a block assembles like the shape it makes" begin
        using Roly: metarules, unwrap, metabonds, nbonds, nfaces, facenormal

        persize(rules, K) = [count(p -> nparticles(p) == k, polygen(rules; maxsize=K)) for k in 1:K]
        outward(s) = s.pose.psi[:, 1]
        # exposed sites that face opposite ways. Unique when the block has one site per direction,
        # which is what lets a bond through it fix the neighbour completely
        function facing(sites)
            pairs = NTuple{2,Int}[]
            for i in eachindex(sites), j in eachindex(sites)
                i < j && isapprox(outward(sites[j]), -outward(sites[i]); atol=1e-8) || continue
                push!(pairs, (i, j))
            end
            return pairs
        end
        # the site directly across the block, with no sideways offset between the two. Needed
        # where a side carries several sites, so that facing alone leaves a choice
        function across(sites)
            return filter(facing(sites)) do (i, j)
                n = outward(sites[i])
                d = sites[i].pose.x - sites[j].pose.x
                return norm(d - dot(d, n) * n) < 1e-8
            end
        end
        # a distinct color per exposed site, bonding only the given pairs
        function keyed(poly, pairing)
            sites = opensites(poly)
            ps = MetaParticleSpecies(poly; colors=1:length(sites))
            pairs = pairing(sites)
            @test sort(collect(Iterators.flatten(pairs))) == 1:length(sites)
            return BindingRules(reduce(vcat, [[1 a 1 b] for (a, b) in pairs]), ps)
        end

        ### wrapping a single particle changes nothing
        for (rules, K) in ((BindingRules([1 1 1 3; 1 2 1 4], UnitSquare), 5),
                           (BindingRules([1 1 1 1], PolygonParticleSpecies(3; colors=fill(1, 3))), 5),
                           (BindingRules([1 1 1 1],
                                         PolyhedronParticleSpecies(Cube(); colors=fill(1, 6))), 4))
            wrapped = metarules(MetaParticleSpecies(Polyform(rules, 1)))
            @test persize(wrapped, K) == persize(rules, K)
        end

        # Which lattice-animal count a system gives is set by its coloring and its dimension. A
        # distinct color on every site leaves the tile no symmetry, so structures are counted up
        # to translation alone -- the *fixed* sequences. One color throughout leaves the tile its
        # own symmetry; in 2D the encoding is a digraph, so a reflection is not an automorphism
        # and the counts are *one-sided*, while in 3D flipping a flat arrangement over is a proper
        # rotation, which makes them *free*. All four below extend correctly by one more term.

        ### two triangles make a rhombus, which tiles the plane the way a square does
        iamonds = BindingRules([1 1 1 1], PolygonParticleSpecies(3; colors=fill(1, 3)))
        rhombus = only(p for p in polygen(iamonds; maxsize=2) if nparticles(p) == 2)
        @test length(opensites(rhombus)) == 4
        # fixed polyominoes, https://oeis.org/A001168
        @test persize(keyed(rhombus, facing), 5) ==
              persize(BindingRules([1 1 1 3; 1 2 1 4], UnitSquare), 5) == [1, 2, 6, 19, 63]

        ### six triangles make a hexagon, which tiles like one under either coloring
        ring = only(p for p in polygen(iamonds; maxsize=6)
                    if nparticles(p) == 6 && length(opensites(p)) == 6)
        # a distinct color per site leaves the block no symmetry: fixed polyhexes,
        # https://oeis.org/A001207
        @test persize(keyed(ring, facing), 4) ==
              persize(BindingRules([1 1 1 4; 1 2 1 5; 1 3 1 6], UnitHexagon), 4) == [1, 3, 11, 44]
        # the colors it inherits are all one, so it keeps the hexagon's full 6-fold symmetry, and
        # the counts are one-sided: https://oeis.org/A006535
        @test symmetrynumber(MetaParticleSpecies(ring)) == 6
        @test persize(metarules(MetaParticleSpecies(ring)), 5) ==
              persize(BindingRules([1 1 1 1], PolygonParticleSpecies(6; colors=fill(1, 6))), 5) ==
              [1, 1, 3, 10, 33]

        ### the same construction in 3D: six triangular prisms make a hexagonal prism, again
        ### fully symmetric. A reflection of a flat arrangement is a rotation about an in-plane
        ### axis here, so the very same coloring gives the free counts rather than the one-sided
        sides(p) = [i for i in 1:nfaces(p) if abs(facenormal(p, i)[3]) < 1e-8]
        sticky(shp, s) = PolyhedronParticleSpecies(shp; colors=[i in s ? 1 : 2 for i in 1:nfaces(shp)])
        tri3, hex3 = Prism(3, 1.0; h=2.0), Prism(6, 1.0; h=2.0)
        prismrules(shp) = (s = sides(shp); BindingRules([1 first(s) 1 first(s)], sticky(shp, s)))
        ring3 = only(p for p in polygen(prismrules(tri3); maxsize=6)
                     if nparticles(p) == 6 && length(opensites(p)) == 6)
        mp3 = MetaParticleSpecies(ring3)
        @test nsites(mp3) == 6
        # the ring is as symmetric as the hexagonal prism it makes, D_6 of order 12
        @test symmetrynumber(mp3) == symmetrynumber(ring3) == 12
        # free polyhexes, https://oeis.org/A000228
        @test persize(metarules(mp3), 4) == persize(prismrules(hex3), 4) == [1, 1, 3, 7]

        # the same split without a block, on plain squares: one-sided in 2D, free in 3D
        # https://oeis.org/A000988 and https://oeis.org/A000105
        @test persize(BindingRules([1 1 1 1], PolygonParticleSpecies(4; colors=fill(1, 4))), 5) ==
              [1, 1, 2, 7, 18]
        @test persize(prismrules(Prism(4, 1.0; h=2.0)), 5) == [1, 1, 2, 5, 12]

        ### a 2x2 block of squares, whose sides carry two sites each
        sqrules = BindingRules([1 1 1 3; 1 2 1 4], UnitSquare)
        block = only(p for p in polygen(sqrules; maxsize=4)
                     if nparticles(p) == 4 && length(opensites(p)) == 8)
        # keeping the colors it inherits lets a block meet its neighbour half a block over,
        # sharing one edge instead of two, which no single square can do
        loose = polygen(metarules(MetaParticleSpecies(block)); maxsize=2)
        @test count(m -> nparticles(m) == 2, loose) > 1
        @test any(m -> length(metabonds(m)) == 1, loose)
        # every one of them is still an assembly of squares
        direct = Set(graphrep(p) for p in polygen(sqrules; maxsize=8))
        @test all(graphrep(unwrap(m)) in direct for m in loose)
        @test all(nparticles(unwrap(m)) == 4nparticles(m) for m in loose)
        @test all(nbonds(unwrap(m)) == 4nparticles(m) + length(metabonds(m)) for m in loose)
        # giving the two sites of a side different colors leaves the offset nothing to bond to,
        # and the block tiles exactly as a square does
        @test persize(keyed(block, across), 5) == [1, 2, 6, 19, 63]

        ### two species, which must not be conflated when the block is unwrapped
        two = BindingRules([1 1 2 3; 1 2 2 4], [UnitSquare, UnitSquare])
        direct2 = Set(graphrep(p) for p in polygen(two; maxsize=6))
        for seed in (p for p in polygen(two; maxsize=2) if nparticles(p) == 2)
            metas = polygen(metarules(MetaParticleSpecies(seed)); maxsize=3)
            @test !isempty(metas)
            @test all(graphrep(unwrap(m)) in direct2 for m in metas)
            @test all(composition(unwrap(m))[1:2] == nparticles(m) .* composition(seed)[1:2]
                      for m in metas)
        end

        ### 3D polycubes, where a block keeps the symmetry of the polyform it wraps
        cuberules = BindingRules([1 1 1 1], PolyhedronParticleSpecies(Cube(); colors=fill(1, 6)))
        direct3 = Set(graphrep(p) for p in polygen(cuberules; maxsize=6))
        for seed in (p for p in polygen(cuberules; maxsize=3) if nparticles(p) > 1)
            mp = MetaParticleSpecies(seed)
            @test symmetrynumber(mp) == symmetrynumber(seed)
            metas = polygen(metarules(mp); maxsize=nparticles(seed) == 2 ? 3 : 2)
            @test all(graphrep(unwrap(m)) in direct3 for m in metas)
        end
    end

    @testset "copy is independent" begin
        ps = MetaParticleSpecies(dimer)
        before = copy(labels(graphrep(ps)))
        cp = copy(ps)
        setcolors!(cp, [11, 12])
        @test [color(bindingsite(ps, i)) for i in 1:2] != [11, 12]
        @test [color(bindingsite(cp, i)) for i in 1:2] == [11, 12]
        @test labels(graphrep(ps)) == before
    end

    @testset "the encoding is checked like any other species" begin
        # `check_encoding` runs in the constructor now, against the polyform's own rotations
        # rather than a group derived from the sites, which could not see the polyform behind them
        @test Roly.check_encoding(mp) === mp
        @test length(Roly.permutationgroup(mp)) == symmetrynumber(mp)

        selfrules = BindingRules([1 1 1 1; 1 2 1 2], UnitSquare)
        sym = polygen(selfrules; maxsize=2)[end]
        @test length(Roly.rotationgroup(sym)) == 2

        # the polyform's rotations are only the candidates. Coloring the two exposed sites apart
        # costs the meta-species a symmetry the polyform keeps, so its group is a subgroup
        for (cols, order) in (([9, 9], 2), ([4, 5], 1))
            ps = MetaParticleSpecies(sym; colors=cols)
            @test Roly.check_encoding(ps) === ps
            @test length(Roly.rotationgroup(ps)) == order
            @test length(Roly.permutationgroup(ps)) == symmetrynumber(ps) == order
        end

        # and so does exposing only one of the two
        one = MetaParticleSpecies(sym, opensites(ParticleSite, sym)[1:1])
        @test length(Roly.rotationgroup(one)) == symmetrynumber(one) == 1
    end

    @testset "recoloring tracks the orbits, not the colors" begin
        # this dimer has no symmetry relating its two open sites, so they sit in different orbits
        # whatever they are colored, and recoloring cannot move the labeling
        ps = MetaParticleSpecies(dimer)
        before = copy(labels(graphrep(ps)))
        setcolors!(ps, [11, 12])
        @test labels(graphrep(ps)) == before
        setcolors!(ps, [7, 7])
        @test labels(graphrep(ps)) == before

        # where a symmetry does relate the two sites, they share an orbit and a label until a
        # coloring tells them apart
        selfrules = BindingRules([1 1 1 1; 1 2 1 2], UnitSquare)
        sym = polygen(selfrules; maxsize=2)[end]
        sitelabels(q) = [labels(graphrep(q))[first(bindingsite(q, i).vertices)] for i in 1:nsites(q)]
        together = MetaParticleSpecies(sym; colors=[9, 9])
        @test length(unique(sitelabels(together))) == 1
        setcolors!(together, [4, 5])
        @test length(unique(sitelabels(together))) == 2
        @test symmetrynumber(together) == 1
    end
end
