@testset "MetaParticleSpecies" begin
    using Roly: setcolors!, collect_open_bindingsites, overlap, cluster
    using Graphs: ne

    # one species chaining through opposite sites, so a dimer keeps one open end of each kind
    chainlike = BindingRules([1 1 1 3], UnitSquare)
    chainstrs = polygen(chainlike; maxsize=4)
    dimer = chainstrs[findfirst(s -> nparticles(s) == 2, chainstrs)]
    open = collect_open_bindingsites(dimer)
    @test length(open) == 2

    mp = MetaParticleSpecies(dimer)
    @test dimension(mp) == 2
    @test nsites(mp) == length(open)
    @test cluster(mp) == dimer
    @test [color(bindingsite(mp, i)) for i in 1:nsites(mp)] == [color(s) for s in open]
    @test [bindingsite(mp, i).pose for i in 1:nsites(mp)] == [s.pose for s in open]
    @test bounding_radius(mp) >= maximum(norm(p.pose.x) for p in dimer.particles)
    # `_ndistincttwists` needs the twist freedoms, which are 1 in 2D
    @test stabilizerorders(mp) == ones(Int, nsites(mp))
    @test _ndistincttwists(bindingsite(mp, 1), bindingsite(mp, 2)) == 1

    @testset "the encoding is the cluster's own graph" begin
        # bound sites stay as interior vertices, so the graph is larger than the site count
        @test nv(graphrep(mp)) == nv(graphrep(dimer))
        @test nv(graphrep(mp)) > nsites(mp)
        @test ne(graphrep(mp)) == ne(graphrep(dimer))

        # a symmetric cluster: two squares joined site1-site1, swapping them is an automorphism
        selfrules = BindingRules([1 1 1 1; 1 2 1 2], UnitSquare)
        selfstrs = polygen(selfrules; maxsize=2)
        sym = selfstrs[findfirst(s -> nparticles(s) == 2, selfstrs)]
        @test symmetrynumber(sym) == 2
        nopen = length(collect_open_bindingsites(sym))
        @test nopen == 2

        # the species inherits the cluster's symmetry exactly -- an encoding over the sites alone
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
        rows = exposablesites(dimer)
        # every unbound site is offered, bound ones never are
        @test length(rows) == 2 * nsites(UnitSquare) - 2
        @test count(r -> !r.inert, rows) == length(open)
        @test all(r -> 1 <= r.particle <= nparticles(dimer), rows)
        # the default exposure is exactly the non-inert rows, in this order
        @test [color(bindingsite(mp, i)) for i in 1:nsites(mp)] ==
              [r.color for r in rows if !r.inert]

        # a site that is inert inside the cluster becomes usable by being named and given a
        # colour the new rules speak
        inert = [(r.particle, r.site) for r in rows if r.inert]
        @test !isempty(inert)
        activated = MetaParticleSpecies(dimer, inert[1:2]; colors=[1, 2])
        @test nsites(activated) == 2
        @test [color(bindingsite(activated, i)) for i in 1:2] == [1, 2]
        # and those sites really do bind under rules written over the new colors
        @test [polyenum(BindingRules([1 1 1 2], activated); maxsize=k)[1] for k in 1:3] ==
              [1, 2, 3]

        # exposure order is the caller's
        pair = [(r.particle, r.site) for r in rows if !r.inert]
        @test [color(bindingsite(MetaParticleSpecies(dimer, reverse(pair)), i)) for i in 1:2] ==
              reverse([r.color for r in rows if !r.inert])

        bound = [(p, k) for p in 1:nparticles(dimer) for k in 1:nsites(UnitSquare)
                 if (p, k) ∉ Set((r.particle, r.site) for r in rows)]
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
        live = [(r.particle, r.site) for r in exposablesites(sym) if !r.inert]
        @test length(live) == 2
        @test symmetrynumber(MetaParticleSpecies(sym, live)) == 2
        @test symmetrynumber(MetaParticleSpecies(sym, live[1:1])) == 1
    end

    @testset "keywords and rejections" begin
        recolored = MetaParticleSpecies(dimer; colors=[4, 9])
        @test [color(bindingsite(recolored, i)) for i in 1:2] == [4, 9]
        @test_throws DimensionMismatch MetaParticleSpecies(dimer; colors=[1])

        # a saturated cluster exposes nothing: the 1-2-3 chain's ends are inert
        chain = BindingRules([1 1 2 3; 2 1 3 3], UnitSquare)
        trimer = let strs = polygen(chain; maxsize=3)
            strs[findfirst(s -> nparticles(s) == 3, strs)]
        end
        @test isempty(collect_open_bindingsites(trimer))
        @test_throws ArgumentError MetaParticleSpecies(trimer)

    end

    @testset "3D: symmetry and twist freedom come from the cluster" begin
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

        # equal colors on the two ends make them interchangeable only if the cluster agrees;
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

    @testset "copy is independent" begin
        ps = MetaParticleSpecies(dimer)
        cp = copy(ps)
        setcolors!(cp, [11, 12])
        @test [color(bindingsite(ps, i)) for i in 1:2] != [11, 12]
        @test [color(bindingsite(cp, i)) for i in 1:2] == [11, 12]
        @test labels(graphrep(ps)) != labels(graphrep(cp))
    end
end
