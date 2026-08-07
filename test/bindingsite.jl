using Roly: BindingSite, shift_vertices, shift_color, isaligned, istouching, isincontact, color,
            standard_offset, contact_pairing
using Graphs, NautyGraphs

@testset "binding site" begin
    tol = sqrt(eps(Float64))
    b = BindingSite(Pose{2}(), 1, 1:1, tol, tol)

    @test color(b) == 1
    @test b == BindingSite(Pose{2}(), 1, 1:1, tol, tol)
    @test b != shift_color(b, 1)
    @test b != shift_vertices(b, 1)
    @test hash(b) == hash(BindingSite(Pose{2}(), 1, 1:1, 0.0, 0.0))

    @test shift_vertices(b, 2).vertices == 3:3
    @test color(shift_color(b, 3)) == 4
    @test b < shift_vertices(b, 1)

    p = Pose(SVector(1.0, 0.0), Angle2d(0.0))
    @test (p * b).pose == p * b.pose
    @test (b * Angle2d(π/4)).pose == b.pose * Angle2d(π/4)

    io = IOBuffer()
    show(io, b)
    @test contains(String(take!(io)), "BindingSite")
    show(io, MIME"text/plain"(), b)
    @test contains(String(take!(io)), "color")

    # 2D contact geometry: sites at same position, rotations differing by π
    b1 = BindingSite(Pose(SVector(1.0, 0.0), Angle2d(0.0)), 1, 1:1, tol, tol)
    b2 = BindingSite(Pose(SVector(1.0, 0.0), Angle2d(Float64(π))), 1, 2:2, tol, tol)
    b3 = BindingSite(Pose(SVector(2.0, 0.0), Angle2d(Float64(π))), 1, 3:3, tol, tol)

    @test istouching(b1, b2)
    @test isaligned(b1, b2)
    @test isincontact(b1, b2)

    @test !istouching(b1, b3)
    @test isaligned(b1, b3)
    @test !isincontact(b1, b3)

    # 3D contact geometry
    b1_3d = BindingSite(Pose(SVector(1.0, 0.0, 0.0), one(RotXYZ)), 1, 1:1, tol, tol)
    b2_3d = BindingSite(Pose(SVector(1.0, 0.0, 0.0), standard_offset(b1_3d).psi), 1, 2:2, tol, tol)
    b3_3d = BindingSite(Pose(SVector(2.0, 0.0, 0.0), standard_offset(b1_3d).psi), 1, 3:3, tol, tol)

    @test istouching(b1_3d, b2_3d)
    @test isaligned(b1_3d, b2_3d)
    @test isincontact(b1_3d, b2_3d)

    @test !istouching(b1_3d, b3_3d)
    @test isaligned(b1_3d, b3_3d)
    @test !isincontact(b1_3d, b3_3d)

    ### contact_pairing
    # Equal-sized sites are matched by the full counter-rotating bijection, anchored on the
    # first vertex of each range: that is the pair whose polyhedron edges coincide.
    @test collect(contact_pairing(1:1, 5:5)) == [1 => 5]
    @test collect(contact_pairing(1:2, 5:6)) == [1 => 5, 2 => 6]
    @test collect(contact_pairing(1:4, 11:14)) == [1 => 11, 2 => 14, 3 => 13, 4 => 12]

    # Sites of different sizes are joined only where their vertices land at the same angle:
    # gcd(k1, k2) pairs. Coprime sizes leave just the anchor.
    @test collect(contact_pairing(1:6, 11:13)) == [1 => 11, 3 => 13, 5 => 12]
    @test collect(contact_pairing(1:4, 11:16)) == [1 => 11, 3 => 14]
    @test collect(contact_pairing(1:4, 11:15)) == [1 => 11]
    for (k1, k2) in [(1, 1), (2, 4), (6, 3), (4, 6), (12, 8), (4, 5), (5, 5), (7, 3)]
        @test length(collect(contact_pairing(1:k1, 11:(10 + k2)))) == gcd(k1, k2)
        @test first(contact_pairing(1:k1, 11:(10 + k2))) == (1 => 11)
    end

    # The point of the partial matching: a bond keeps the symmetry the two sites have in
    # common. Joining a k1-fold site to a k2-fold one must leave gcd(k1, k2) turns, since a
    # rotation about the bond axis has to be a symmetry of both.
    function dimer_symmetrynumber(k1, k2)
        g = NautyDiGraph(k1 + k2; vertex_labels=Cint[fill(1, k1); fill(2, k2)])
        for i in 1:k1
            add_edge!(g, i, mod1(i + 1, k1))
        end
        for j in 1:k2
            add_edge!(g, k1 + j, k1 + mod1(j + 1, k2))
        end
        for (a, b) in contact_pairing(1:k1, (k1 + 1):(k1 + k2))
            add_edge!(g, a, b)
            add_edge!(g, b, a)
        end
        return convert(Int, nauty(g)[2].n)
    end
    for (k1, k2) in [(3, 3), (4, 4), (6, 6), (6, 3), (3, 6), (4, 6), (6, 4), (12, 8), (8, 6), (4, 5)]
        @test dimer_symmetrynumber(k1, k2) == gcd(k1, k2)
    end
end
