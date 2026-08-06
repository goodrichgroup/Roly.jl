using Roly:
    PatchyParticleSpecies,
    PatchyDisk,
    nsites,
    dimension,
    isconvex,
    numtype,
    bindingsites,
    graphrep,
    setcolors!,
    color,
    could_contact,
    overlap,
    symmetrynumber

@testset "PatchyParticleSpecies" begin
    # basic construction for n > 2
    for n in (3, 4, 6)
        ps = PatchyDisk(range(0, 2π; length=n + 1)[1:n])
        @test nsites(ps) == n
        @test dimension(ps) == 2
        @test isconvex(ps)
        @test numtype(ps) == Float64
        @test nv(graphrep(ps)) == n
        for i in 1:n
            @test norm(bindingsites(ps, i).pose.x) ≈ 1.0 atol = 1e-10
        end
    end

    ps32 = PatchyDisk([0.0f0, Float32(2π / 3), Float32(4π / 3)])
    @test numtype(ps32) == Float32

    # n=1: single site, 1-vertex graph
    ps1 = PatchyDisk([0.0])
    @test nsites(ps1) == 1
    @test dimension(ps1) == 2
    @test nv(graphrep(ps1)) == 1
    @test norm(bindingsites(ps1, 1).pose.x) ≈ 1.0 atol = 1e-10

    # n=2: uses a 4-vertex graph to encode orientation
    ps2 = PatchyDisk([0.0, π])
    @test nsites(ps2) == 2
    @test dimension(ps2) == 2
    @test nv(graphrep(ps2)) == 4
    @test bindingsites(ps2, 1).vertices == 1:2
    @test bindingsites(ps2, 2).vertices == 3:4
    @test bindingsites(ps2, 1).pose.x ≈ SVector(1.0, 0.0) atol = 1e-10
    @test bindingsites(ps2, 2).pose.x ≈ SVector(-1.0, 0.0) atol = 1e-10

    # custom radius
    r = 2.5
    psr = PatchyDisk([0.0, π / 2, π, 3π / 2], r)
    @test nsites(psr) == 4
    for i in 1:4
        @test norm(bindingsites(psr, i).pose.x) ≈ r atol = 1e-10
    end

    # colors
    ps_col = PatchyDisk([0.0, π / 2, π, 3π / 2]; colors=[1, 2, 1, 2])
    @test color(bindingsites(ps_col, 1)) == 1
    @test color(bindingsites(ps_col, 2)) == 2
    @test color(bindingsites(ps_col, 3)) == 1
    @test color(bindingsites(ps_col, 4)) == 2

    # symmetrynumber: all distinct → 1, all equal → n
    @test symmetrynumber(PatchyDisk([0.0, 2π / 3, 4π / 3])) == 1
    @test symmetrynumber(PatchyDisk([0.0, 2π / 3, 4π / 3]; labels=[1, 1, 1])) == 3
    @test symmetrynumber(PatchyDisk([0.0, π / 2, π, 3π / 2]; labels=[1, 1, 1, 1])) == 4

    # copy and setcolors!
    ps = PatchyDisk([0.0, 2π / 3, 4π / 3])
    ps_c = copy(ps)
    @test nsites(ps_c) == nsites(ps)
    setcolors!(ps_c, [10, 20, 30])
    @test color(bindingsites(ps_c, 1)) == 10
    @test color(bindingsites(ps, 1)) != 10
    @test_throws ArgumentError setcolors!(ps, [1, 2])

    io = IOBuffer()
    show(io, ps)
    @test contains(String(take!(io)), "PatchyParticleSpecies")

    # could_contact and overlap
    pd = PatchyDisk([0.0, π])
    id = Pose{2}()
    far = Pose(SVector(100.0, 0.0), Angle2d(0.0))
    bonded = Pose(SVector(2.0, 0.0), Angle2d(0.0))   # centers at 2r, bonded config
    close = Pose(SVector(1.5, 0.0), Angle2d(0.0))

    @test could_contact(pd => id, pd => far)
    @test overlap(pd => id, pd => id)
    @test overlap(pd => id, pd => close)
    @test !overlap(pd => id, pd => bonded)
    @test !overlap(pd => id, pd => far)
end
