@testset "MetaParticleSpecies" begin
    using Roly: setcolors!, collect_open_bindingsites, overlap, cluster

    # one species chaining through opposite sites, so a dimer keeps one open end of each kind
    chainlike = BindingRules([1 1 1 3], UnitSquare)
    dimer = let strs = polygen(chainlike; maxsize=2)
        strs[findfirst(s -> nparticles(s) == 2, strs)]
    end
    open = collect_open_bindingsites(dimer)
    @test length(open) == 2

    mp = MetaParticleSpecies(dimer)
    @test dimension(mp) == 2
    @test nsites(mp) == length(open)
    @test cluster(mp) == dimer
    # the sites are the cluster's open sites, in order, with their colors
    @test [color(bindingsites(mp, i)) for i in 1:nsites(mp)] == [color(s) for s in open]
    @test [bindingsites(mp, i).pose for i in 1:nsites(mp)] == [s.pose for s in open]
    # the bounding radius has to contain every constituent particle
    @test bounding_radius(mp) >= maximum(norm(p.pose.x) for p in dimer.particles)

    # a cluster's shape is not its sites, so no symmetry is claimed: labels distinct, sigma 1,
    # and stabilizers 1 so that `nphases` keeps every attachment phase
    @test Roly.labels(graphrep(mp)) == collect(1:nsites(mp))
    @test symmetrynumber(mp) == 1
    @test sitestabilizers(mp) == ones(Int, nsites(mp))
    @test all(nphases(bindingsites(mp, i), bindingsites(mp, j)) == 1
              for i in 1:nsites(mp), j in 1:nsites(mp))

    @testset "keywords and rejections" begin
        recolored = MetaParticleSpecies(dimer; colors=[4, 9])
        @test [color(bindingsites(recolored, i)) for i in 1:2] == [4, 9]
        # a declared symmetry is taken as given rather than re-derived
        tied = MetaParticleSpecies(dimer; labels=[1, 1])
        @test Roly.labels(graphrep(tied)) == [1, 1]

        @test_throws DimensionMismatch MetaParticleSpecies(dimer; colors=[1])
        @test_throws DimensionMismatch MetaParticleSpecies(dimer; labels=[1, 2, 3])

        # a saturated cluster exposes nothing: the 1-2-3 chain's ends are inert
        chain = BindingRules([1 1 2 3; 2 1 3 3], UnitSquare)
        trimer = let strs = polygen(chain; maxsize=3)
            strs[findfirst(s -> nparticles(s) == 3, strs)]
        end
        @test isempty(collect_open_bindingsites(trimer))
        @test_throws ArgumentError MetaParticleSpecies(trimer)

        # 3D needs a dart encoding: a face's twist freedom cannot be declared away
        cubes = BindingRules([1 1 1 2], UnitCube)
        cubedimer = let strs = polygen(cubes; maxsize=2)
            strs[findfirst(s -> nparticles(s) == 2, strs)]
        end
        @test_throws ArgumentError MetaParticleSpecies(cubedimer)
    end

    @testset "setcolors! keeps the declared labeling" begin
        ps = MetaParticleSpecies(dimer)
        before = copy(Roly.labels(graphrep(ps)))
        setcolors!(ps, [7, 8])
        @test [color(bindingsites(ps, i)) for i in 1:2] == [7, 8]
        # the generic method would re-derive labels and stabilizers from the site poses here
        @test Roly.labels(graphrep(ps)) == before
        @test sitestabilizers(ps) == [1, 1]
        # even a coloring that makes the two site poses look alike claims no new symmetry
        setcolors!(ps, [5, 5])
        @test Roly.labels(graphrep(ps)) == before
        @test symmetrynumber(ps) == 1
        @test_throws ArgumentError setcolors!(ps, [1, 2, 3])
    end

    @testset "overlap delegates to the constituents" begin
        P = Roly.posetype(mp)
        @test overlap(mp => one(P), mp => one(P))
        far = Pose(SVector(50.0, 0.0), Angle2d(0.0))
        @test !overlap(mp => one(P), mp => far)
        # a block one lattice step away along the chain axis clears its neighbour
        step = Pose(SVector(2.0, 0.0), Angle2d(0.0))
        @test !overlap(mp => one(P), mp => step)
        # while half a step in overlaps it
        @test overlap(mp => one(P), mp => Pose(SVector(0.5, 0.0), Angle2d(0.0)))
    end

    @testset "assembling blocks out of blocks" begin
        ps = MetaParticleSpecies(dimer)
        # the two ends are distinguishable, so blocks chain in exactly one way per length
        metasys = BindingRules([1 1 1 2], ps)
        @test [polyenum(metasys; maxsize=k)[1] for k in 1:4] == [1, 2, 3, 4]
        blocks = polygen(metasys; maxsize=3)
        @test sort(nparticles.(blocks)) == [1, 2, 3]
        # a meta-cluster counts blocks, and each block carries the wrapped cluster's particles
        triple = blocks[findfirst(s -> nparticles(s) == 3, blocks)]
        @test nparticles(triple) * nparticles(dimer) == 6
    end

    @testset "copy is independent" begin
        ps = MetaParticleSpecies(dimer)
        cp = copy(ps)
        setcolors!(cp, [11, 12])
        @test [color(bindingsites(ps, i)) for i in 1:2] != [11, 12]
        @test [color(bindingsites(cp, i)) for i in 1:2] == [11, 12]
    end
end
