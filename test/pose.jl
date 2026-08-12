using StaticArrays, Rotations

@testset "pose" begin
    p2 = Pose{2,Float64}()
    @test eltype(p2) === Float64
    @test eltype(typeof(p2)) === Float64
    @test typeof(p2.psi) === Angle2d{Float64}

    @test isapprox(p2.x, zero(SVector{2,Float64}); atol=1e-12)
    @test isapprox(p2.psi, one(Angle2d{Float64}); atol=1e-12)

    p3 = Pose{3,Float32}()
    @test eltype(p3) === Float32
    @test eltype(typeof(p3)) === Float32
    @test typeof(p3.psi) === RotMatrix3{Float32}

    @test isapprox(p3.x, zero(SVector{3,Float64}); atol=1e-12)
    @test isapprox(p3.psi,  one(RotMatrix3{Float64}); atol=1e-12)

    @test Pose{2}() == p2
    @test Pose{3}() == p3

    @test dimension(p2) == 2
    @test dimension(p3) == 3

    ####

    X = SVector(1.0, 0.0, 0.0)
    Y = SVector(0.0, 1.0, 0.0)
    Z = SVector(0.0, 0.0, 1.0)
    RY = RotXYZ(0, π/2, 0)
    RZ = RotXYZ(0, 0, π)

    p1 = Pose(X, one(RotXYZ))
    p2 = Pose(X, one(RotXYZ))

    @test isapprox(p1 * p1, Pose(2X, one(RotXYZ)); atol=1e-12)
    @test isapprox(RY * p1, Pose(RY * X, RY*one(RotXYZ)); atol=1e-12)
    @test isapprox(p1 * RY, Pose(X, RY*one(RotXYZ)); atol=1e-12)

    p3 = Pose(Y, one(RotXYZ))
    p4 = Pose(X, RZ)

    @test isapprox(p3 * p4, Pose(X + Y, RZ); atol=1e-12)
    @test isapprox(p4 * p3, Pose(-Y + X, RZ); atol=1e-12)

    p5 = Pose(X, RY)
    p6 = Pose(Z, RZ)

    @test isapprox(p5 * p6, Pose(X + X, RY * RZ); atol=1e-12)
    @test isapprox(p6 * p5, Pose(-X + Z, RZ * RY); atol=1e-12)
    @test isapprox(p5 / p6, Pose(-X + X, RY * RZ'); atol=1e-12)
    @test isapprox(p6 / p5, Pose(-Z + Z, RZ * RY'); atol=1e-12)
    @test isapprox(p5 \ p6, Pose(-X - Z, RY' * RZ); atol=1e-12)
    @test isapprox(p6 \ p5, Pose(-X - Z, RZ' * RY); atol=1e-12)
end
@testset "pose rotation types convert" begin
    # Rotation types are not closed under multiplication, so a pose can arrive parameterised
    # differently from the field it is being stored into. 3D species all use `RotMatrix3` so
    # this does not come up in the hot path, but `Pose` should still behave as a parametric
    # type rather than erroring.
    p = Pose(SVector(1.0, 2.0, 3.0), RotXYZ(0.1, 0.2, 0.3))
    q = convert(Pose{3,Float64,RotMatrix3{Float64}}, p)
    @test q isa Pose{3,Float64,RotMatrix3{Float64}}
    @test isapprox(q.psi, p.psi; atol=1e-12)
    @test q.x == p.x
    @test convert(Pose{3,Float64,RotXYZ{Float64}}, q) isa Pose{3,Float64,RotXYZ{Float64}}
    @test convert(typeof(p), p) === p
end
