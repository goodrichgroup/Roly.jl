using Roly: BindingSite, shift_vertices, shift_color, isaligned, istouching, isincontact, color, standard_offset

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
end
