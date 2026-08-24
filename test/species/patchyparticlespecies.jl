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
    symmetrynumber,
    site_symmetry


using Roly: PatchySphere, Polyhedron, Tetrahedron, Cube, Dodecahedron, Prism,
            nfaces, facecentroid, geometriclabels, rotationgroup, symmetrynumber,
            graphrep, edgemidpoint

using Roly: PatchyParticleSpecies, BindingRules, Polyform, raise!, collect_compatible_pairs,
            nparticles

using Rotations: RotMatrix3, rotation_angle
using Graphs: cycle_digraph
using NautyGraphs: NautyDiGraph
using LinearAlgebra: normalize, dot, det, norm

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
    @test nv(graphrep(ps2)) == 2
    @test bindingsites(ps2, 1).vertices == 1:1
    @test bindingsites(ps2, 2).vertices == 2:2

    # Two equivalent patches give a symmetry number of 2.
    @test symmetrynumber(PatchyDisk([0.0, π]; colors=[1, 1])) == 2
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

    # symmetrynumber: all distinct -> 1, all equal -> n
    @test symmetrynumber(PatchyDisk([0.0, 2π / 3, 4π / 3])) == 1
    @test symmetrynumber(PatchyDisk([0.0, 2π / 3, 4π / 3]; colors=[1, 1, 1])) == 3
    @test symmetrynumber(PatchyDisk([0.0, π / 2, π, 3π / 2]; colors=[1, 1, 1, 1])) == 4

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


    # spheres

    for (name, shp, np, order) in [("T", Tetrahedron(), 4, 12), ("O", Cube(), 6, 24),
                                   ("I", Dodecahedron(), 12, 60), ("D5", Prism(5), 7, 10)]
        ps = PatchySphere(shp, 2.0)
        @test dimension(ps) == 3
        @test nsites(ps) == np
        @test numtype(ps) === Float64

        # Patches sit on the sphere, along the face centroid directions.
        for i in 1:np
            b = bindingsites(ps, i)
            @test isapprox(norm(b.pose.x), 2.0)
            @test isapprox(normalize(b.pose.x), normalize(facecentroid(shp, i)); atol=1e-10)
            # On a sphere the patch normal is radial, and local z is the tangential part
            # of the direction to the face's first edge.
            @test isapprox(b.pose.psi[:, 1], normalize(b.pose.x); atol=1e-10)
            @test isapprox(dot(b.pose.psi[:, 3], b.pose.psi[:, 1]), 0; atol=1e-10)
            @test isapprox(det(b.pose.psi), 1; atol=1e-10)
        end

        # Same labeling rules and the same graph as the polyhedron species.
        @test symmetrynumber(ps) == 1
        @test symmetrynumber(PatchySphere(shp, 2.0; colors=geometriclabels(shp))) == order
        @test order == length(rotationgroup(shp))
        @test nv(graphrep(ps)) == np
        @test nv(graphrep(dartsphere(shp, 2.0))) == 2 * Roly.nedges(shp)
    end

    # Naming a rotation group resolves to the `Polyhedron` realizes it with.
    @test nsites(PatchySphere(Tetrahedron())) == 4
    @test nsites(PatchySphere(Cube())) == 6
    @test nsites(PatchySphere(Dodecahedron())) == 12
    @test nsites(PatchySphere(Prism(5))) == 7
    @test nsites(PatchySphere(Pyramid(6))) == 7
    @test_throws ArgumentError PatchySphere(Cube(); colors=1:5)

    # Spheres overlap by radius alone.
    ps = PatchySphere(Cube(), 1.0)
    id = Pose{3,Float64,RotMatrix3{Float64}}(SVector(0.0, 0.0, 0.0), one(RotMatrix3{Float64}))
    apart(d) = Pose{3,Float64,RotMatrix3{Float64}}(SVector(d, 0.0, 0.0), one(RotMatrix3{Float64}))
    @test overlap(ps => id, ps => id)
    @test overlap(ps => id, ps => apart(1.5))
    @test !overlap(ps => id, ps => apart(2.0))
    # Patchy particles skip the bounding-sphere pre-check, since it would just repeat the
    # overlap test, so could_contact is unconditionally true.
    @test could_contact(ps => id, ps => apart(50.0))

    # encoding 

    # Evenly spaced identical patches really do have that symmetry
    @test symmetrynumber(PatchyDisk([0.0, 2π / 3, 4π / 3]; colors=[1, 1, 1])) == 3
    @test site_symmetry(PatchyDisk([0.0, 2π / 3, 4π / 3]; colors=[1, 1, 1])) == 3
    @test symmetrynumber(PatchyDisk([0.0, π]; colors=[1, 1])) == 2
    @test site_symmetry(PatchyDisk([0.0, π]; colors=[1, 1])) == 2

    # Unevenly spaced ones do not, no matter what colors
    @test symmetrynumber(PatchyDisk([0.0, 0.5, 3.0]; colors=[1, 1, 1])) == 1
    @test length(unique(labels(graphrep(PatchyDisk([0.0, 0.5, 3.0]; colors=[1, 1, 1]))))) == 3
    @test symmetrynumber(PatchyDisk([0.0, 0.5, 3.0])) == 1

    # The general constructor takes the graph's structure but derives its labels, so a
    # labeling that does not describe the arrangement cannot be supplied in the first place.
    # Three unevenly spaced patches come out distinct however the graph was labeled.
    uneven = [SVector(cos(t), sin(t)) for t in (0.0, 0.5, 3.0)]
    for labs in (Cint[1, 1, 1], Cint[1, 2, 3])
        ps = PatchyParticleSpecies(NautyDiGraph(cycle_digraph(3); vertex_labels=labs), 1.0, uneven)
        @test symmetrynumber(ps) == site_symmetry(ps) == 1
        @test labels(graphrep(ps)) == Cint[1, 2, 3]
    end
    # Evenly spaced patches of one color are equivalent
    even = [SVector(cos(t), sin(t)) for t in (0.0, 2π / 3, 4π / 3)]
    ps = PatchyParticleSpecies(NautyDiGraph(cycle_digraph(3)), 1.0, even; colors=fill(1, 3))
    @test symmetrynumber(ps) == site_symmetry(ps) == 3


    # twists

    # A patch is one graph vertex, so it pins no turn about its own normal and its frame fully 
    # determines the bond
    pos = [SVector(1.0, 0.0, 0.0), SVector(-1.0, 0.0, 0.0),
           SVector(0.0, 1.0, 0.0), SVector(0.0, -1.0, 0.0)]
    function dimer(twists)
        ps = PatchyParticleSpecies(NautyDiGraph(cycle_digraph(4)), 1.0, pos, twists;
                                   colors=[1, 2, 3, 4])
        @test all(i -> Roly.bindingsites(ps, i).sitesym == 1, 1:4)
        rules = BindingRules([1 1 1 2], ps)
        poly = Polyform(rules, 1)
        for (site, loc, r) in collect_compatible_pairs(poly)
            trial = copy(poly)
            ismissing(raise!(trial, site, loc, r)) && continue
            @test r == 0                      # sitesym 1 on both sides, so one twist
            return trial.particles[2].pose
        end
        return nothing
    end

    straight = dimer(zeros(4))
    turned = dimer([0.0, π / 2, 0.0, 0.0])
    @test !isnothing(straight) && !isnothing(turned)
    # Same place, since the twist is about the bond axis...
    @test isapprox(straight.x, turned.x; atol=1e-8)
    # ...and a quarter turn apart, which is the point of it.
    @test !isapprox(straight.psi, turned.psi; atol=1e-8)
    @test isapprox(rotation_angle(RotMatrix3(straight.psi * inv(turned.psi))), π / 2; atol=1e-8)
end
