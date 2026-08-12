using Roly: is_cutset, blockdiag!

@testset "utils" begin
    # is_cutset
    g = NautyGraph(3)
    add_edge!(g, 1, 2)
    add_edge!(g, 2, 3)

    @test is_cutset(g, [2]) == true
    @test is_cutset(g, [1]) == false

    add_edge!(g, 1, 3)

    @test is_cutset(g, [2]) == false
    @test is_cutset(g, [1, 2]) == false
    @test is_cutset(g, [1, 2, 3]) == false

    # blockdiag!
    h = NautyGraph(; vertex_labels=[1,2,3,4])
    add_edge!(h, 1, 3)
    add_edge!(h, 2, 4)

    g_old = copy(g)
    b = blockdiag!(g, h)

    @test b[1:3] == g_old
    @test b[4:7] == h
    @test labels(b) == vcat(labels(g_old), labels(h))
    for i in 1:3, j in 4:7
        @test has_edge(g, i, j) == has_edge(g, j, i) == false
    end
end
@testset "separating axes" begin
    using Roly: sat_overlap, edgenormals
    id = Pose{2}()
    at(x, y, θ=0.0) = Pose(SVector(x, y), Angle2d(θ))
    square(a) = [SVector(a/2, a/2), SVector(a/2, -a/2), SVector(-a/2, -a/2), SVector(-a/2, a/2)]
    unit = square(1.0)
    ov(p1, p2, skin=0.0; cs=unit) =
        sat_overlap(Iterators.flatten((edgenormals(cs, p1), edgenormals(cs, p2))),
                    cs, p1, cs, p2, skin)

    @test ov(id, id)
    @test ov(id, at(0.5, 0.0))
    @test !ov(id, at(1.5, 0.0))
    @test !ov(id, at(5.0, 5.0))
    # Exactly touching: the test is strict, so zero clearance still counts as overlapping,
    # and any positive skin separates them.
    @test ov(id, at(1.0, 0.0))
    @test !ov(id, at(1.0, 0.0), 1e-9)
    # A 45 degree turn reaches further along x, so a gap that was clear no longer is.
    @test !ov(id, at(1.2, 0.0))
    @test ov(id, at(1.2, 0.0, π / 4))

    # `skin` is a clearance in length units, and only means one because the axes are
    # normalised: an unnormalised edge normal would scale the tolerance by the edge length.
    # Two squares overlapping by 0.1 must therefore give the same answer at the same skin
    # whatever their size, which they would not if the axis magnitude leaked in.
    for (a, d) in ((1.0, 0.9), (2.0, 1.9), (5.0, 4.9))
        cs = square(a)
        @test ov(id, at(d, 0.0), 0.05; cs)        # overlap 0.1 exceeds the clearance
        @test !ov(id, at(d, 0.0), 0.2; cs)        # ...and does not
    end

    # Degenerate candidates carry no information and must not be read as separating.
    @test sat_overlap((SVector(0.0, 0.0),), unit, id, unit, id, 0.0)
end
