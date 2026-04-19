using Roly: BindingSite, shift_vertices, shift_color, isaligned, istouching, isincontact, color, standard_offset

@testset "binding site" begin
    tol = sqrt(eps(Float64))
    b = BindingSite(Pose{2}(), 1, 1:1, tol, tol)
    @test color(b) == 1

    b2 = shift_vertices(b, 1)
    @test b2.vertices == 2:2

    b3 = shift_color(b, 1)
    @test color(b3) == 2

    ### 

    b1 = BindingSite(Pose(SVector(1.0, 0.0, 0.0), one(RotXYZ)), 1, 1:1, tol, tol)
    b2 = BindingSite(Pose(SVector(2.0, 0.0, 0.0), one(RotXYZ)) * standard_offset(b1).psi, 1, 1:1, tol, tol)
    b3 = BindingSite(Pose(SVector(1.0, 0.0, 0.0), RotXYZ(π, 0, 0)), 1, 1:1, tol, tol)
    b4 = BindingSite(Pose(SVector(1.0, 0.0, 0.0), standard_offset(b1).psi), 1, 1:1, tol, tol)

    @test isaligned(b1, b2)
    @test !istouching(b1, b2)
    @test !isincontact(b1, b2)

    @test !isaligned(b1, b3)
    @test istouching(b1, b3)
    @test !isincontact(b1, b3)

    @test isaligned(b1, b4)
    @test istouching(b1, b4)
    @test isincontact(b1, b4)
end
