using Roly
using Roly: Polyhedron, Tetrahedron, Cube, Octahedron, Dodecahedron, Icosahedron,
            Pyramid, Prism, Antiprism,
            corners, ncorners, faces, facevertices, nfaces, nedges, facedegree,
            facecentroid, facecentroids, facenormal, facenormals, edgemidpoint,
            minedgelength, bounding_radius, inradius,
            rotationgroup, geometriclabels, dartencoding, cycleencoding,
            facegauge, RotationGroup, Cyclic, Dihedral, Tetrahedral, Octahedral,
            Icosahedral, grouporder
using Graphs, NautyGraphs, LinearAlgebra, StaticArrays

@testset "encoding" begin
    symnum(g) = convert(Int, nauty(g)[2].n)
    solids = [
        ("Tetrahedron", Tetrahedron(), 4, 6, 4, 12),
        ("Cube", Cube(), 8, 12, 6, 24),
        ("Octahedron", Octahedron(), 6, 12, 8, 24),
        ("Dodecahedron", Dodecahedron(), 20, 30, 12, 60),
        ("Icosahedron", Icosahedron(), 12, 30, 20, 60),
        # Pyramid(n) is C_n. Note Pyramid(3) is combinatorially a tetrahedron but its default
        # height makes the lateral edges longer than the base edges, so its true group is C_3.
        ("Pyramid(3)", Pyramid(3), 4, 6, 4, 3),
        ("Pyramid(5)", Pyramid(5), 6, 10, 6, 5),
        ("Pyramid(7)", Pyramid(7), 8, 14, 8, 7),
        # Prism(n) is D_n, except Prism(4) whose default height makes it a cube.
        ("Prism(3)", Prism(3), 6, 9, 5, 6),
        ("Prism(4)", Prism(4), 8, 12, 6, 24),
        ("Prism(4,h=2)", Prism(4, 1.0; h=2.0), 8, 12, 6, 8),
        ("Prism(5)", Prism(5), 10, 15, 7, 10),
        ("Prism(8)", Prism(8), 16, 24, 10, 16),
        # Antiprism(n) is D_n, except the uniform Antiprism(3) which is a regular octahedron.
        ("Antiprism(3)", Antiprism(3), 6, 12, 8, 24),
        ("Antiprism(4)", Antiprism(4), 8, 16, 10, 8),
        ("Antiprism(5)", Antiprism(5), 10, 20, 12, 10),
    ]

    for (name, p, ncrn, nedg, nfac, _) in solids
        @test ncorners(p) == ncrn
        @test nedges(p) == nedg
        @test nfaces(p) == nfac
        # Euler characteristic of a sphere.
        @test ncorners(p) - nedges(p) + nfaces(p) == 2
        # Every face's outward normal points away from the origin.
        @test all(dot(facenormal(p, i), facecentroid(p, i)) > 0 for i in 1:nfac)
        @test length(facecentroids(p)) == nfac
        @test length(facenormals(p)) == nfac
        @test all(isapprox(norm(n), 1) for n in facenormals(p))
        @test 0 < inradius(p) <= bounding_radius(p)
        # The corner centroid sits at the origin, so the rotation group acts about it.
        @test isapprox(sum(corners(p)) / ncrn, zero(SVector{3,Float64}); atol=1e-12)
    end

    @test isapprox(minedgelength(Cube(2.5)), 2.5)
    @test isapprox(bounding_radius(Cube(2.0)), sqrt(3.0))
    @test facedegree(Prism(5), 1) in (5, 4)
    @test length(facevertices(Cube(), 1)) == 4
    # Prism(4) with the default height is a cube, and Antiprism(3) a regular octahedron.
    @test isapprox(minedgelength(Prism(4)), 1.0)
    @test length(unique(round.([norm(corners(Antiprism(3))[f[k]] - corners(Antiprism(3))[f[mod1(k+1, 3)]])
                                for f in faces(Antiprism(3)) for k in 1:3]; digits=8))) == 1

    p = Cube()
    cs, fs = corners(p), faces(p)

    # The library solids all pass.
    for (_, q, _, _, _, _) in solids
        @test Polyhedron(corners(q), faces(q)) isa Polyhedron
    end

    # A single reversed face breaks orientability and must be rejected.
    bad = [copy(f) for f in fs]
    reverse!(bad[1])
    @test_throws ArgumentError Polyhedron(cs, bad)

    # A missing face leaves edges shared by only one face.
    @test_throws ArgumentError Polyhedron(cs, [copy(f) for f in fs[1:(end - 1)]])
    # Degenerate and out-of-range faces.
    @test_throws ArgumentError Polyhedron(cs, [[1, 2], [1, 2, 3], [1, 2, 4], [2, 3, 4]])
    @test_throws ArgumentError Polyhedron(cs, [[1, 1, 2], [1, 2, 3], [1, 2, 4], [2, 3, 4]])
    @test_throws ArgumentError Polyhedron(cs, [[1, 2, 99], [1, 2, 3], [1, 2, 4], [2, 3, 4]])

    for (name, p, _, _, _, order) in solids
        group = rotationgroup(p)
        @test length(group) == order
        # Every element is a proper rotation.
        @test all(isapprox(det(R), 1; atol=1e-8) for R in group)
        # Closed under composition, so it really is a group.
        @test all(any(isapprox(R1 * R2, R; atol=1e-8) for R in group) for R1 in group, R2 in group)
    end

    for (name, p, _, nedg, nfac, order) in solids
        g, ranges = dartencoding(p)

        # One vertex per dart, ranges contiguous and matching the face degrees.
        @test nv(g) == 2nedg
        @test length(ranges) == nfac
        @test all(length(ranges[i]) == facedegree(p, i) for i in 1:nfac)
        @test reduce(vcat, collect.(ranges)) == 1:2nedg
        # One arc per dart for the face cycles, plus two per edge for the pairings.
        @test ne(g) == 2nedg + 2nedg

        # All labels distinct: symmetry number 1, as in 2D.
        @test symnum(g) == 1
        # Labels by geometric orbit: exactly the solid's rotation group.
        @test symnum(dartencoding(p; labels=geometriclabels(p))[1]) == order

        @test_throws ArgumentError dartencoding(p; labels=fill(1, nfac + 1))
    end

    # Labels grouping faces that are not one geometric orbit recover the combinatorial
    # symmetry of the face lattice, which can exceed the geometric group.
    @test symnum(dartencoding(Pyramid(3); labels=fill(1, 4))[1]) == 12
    @test length(rotationgroup(Pyramid(3))) == 3

    # Partial labelling selects the subgroup preserving it: distinguishing the two caps of a
    # cube from its four sides leaves the 4-fold axis and the 2-fold axes through it, D_4.
    # The faces come out of the hull search in no particular order, so pick them by normal.
    caps = [abs(n[3]) > 0.5 ? 2 : 1 for n in facenormals(Cube())]
    @test sort(caps) == [1, 1, 1, 1, 2, 2]
    @test symnum(dartencoding(Cube(); labels=caps)[1]) == 8
    @test geometriclabels(Prism(4, 1.0; h=2.0)) |> x -> length(unique(x)) == 2
    # One face singled out leaves only the rotations about its normal.
    @test symnum(dartencoding(Cube(); labels=[2, 1, 1, 1, 1, 1])[1]) == 4
    @test symnum(dartencoding(Dodecahedron(); labels=[2; fill(1, 11)])[1]) == 5
    @test symnum(dartencoding(Icosahedron(); labels=[2; fill(1, 19)])[1]) == 3

    # A face's own rotational symmetry about its normal. Equal to the degree only for a regular
    # face, and a divisor of it otherwise, so it cannot be read off the dart count.
    for p in (Tetrahedron(), Cube(), Octahedron(), Dodecahedron(), Icosahedron(), Antiprism(4))
        @test facegauge(p) == [facedegree(p, i) for i in 1:nfaces(p)]
    end
    # A pyramid's sides are isosceles triangles: no rotational symmetry at all, despite degree 3.
    @test facegauge(Pyramid(5)) == [5, 1, 1, 1, 1, 1]
    # The case the whole registration story turns on: a prism's side faces are squares when
    # h == a and mere rectangles otherwise, so their gauge halves while the degree stays 4.
    @test facegauge(Prism(3)) == [3, 4, 4, 4, 3]
    @test facegauge(Prism(3, 1.0; h=2.0)) == [3, 2, 2, 2, 3]
    @test facegauge(Prism(6, 1.0; h=2.0)) == [6, 2, 2, 2, 2, 2, 2, 6]
    # A box with three distinct edge lengths has rectangular faces throughout.
    box3 = Polyhedron([SVector(x, y, z) for x in (-1.0, 1.0) for y in (-2.0, 2.0) for z in (-3.0, 3.0)])
    @test facegauge(box3) == fill(2, 6)
    # It is a property of the face, not of the solid: a triangular prism is only 2-fold about a
    # square side face, but the face is still a square.
    @test facegauge(Prism(3), 2) == 4

    g, ranges = cycleencoding(6)
    @test nv(g) == 6
    @test ranges == [i:i for i in 1:6]
    @test symnum(g) == 1
    @test symnum(cycleencoding(6; labels=fill(1, 6))[1]) == 6

    # Two sites give a plain 2-cycle, one vertex each: a pair of opposite arcs inside a
    # particle is fine now that bonds are recognised by joining different particles.
    g2, ranges2 = cycleencoding(2)
    @test nv(g2) == 2
    @test ne(g2) == 2
    @test ranges2 == [1:1, 2:2]
    @test symnum(g2) == 1
    @test symnum(cycleencoding(2; labels=[1, 1])[1]) == 2
    # And a single site is a single vertex.
    g1, ranges1 = cycleencoding(1)
    @test nv(g1) == 1
    @test ranges1 == [1:1]
    @test symnum(g1) == 1

    @test_throws ArgumentError cycleencoding(4; labels=[1, 2, 3])
    @test_throws ArgumentError cycleencoding(0)

    # The sparse encoding agrees with the dart encoding exactly when labels are all
    # distinct, which is the condition under which it may be substituted.
    for (_, p, _, _, nfac, _) in solids
        @test symnum(cycleencoding(nfac)[1]) == symnum(dartencoding(p)[1]) == 1
    end

    # Naming a rotation group gives a solid realizing it, and the group's order is what the
    # solid's own rotations and its dart encoding both come out at.
    for group in [Tetrahedral(), Octahedral(), Icosahedral(),
                  Cyclic(3), Cyclic(6), Dihedral(3), Dihedral(5)]
        p = Polyhedron(group)
        @test length(rotationgroup(p)) == grouporder(group)
        @test symnum(dartencoding(group; labels=geometriclabels(p))[1]) == grouporder(group)
    end
    @test grouporder.([Cyclic(4), Dihedral(4), Tetrahedral(), Octahedral(), Icosahedral()]) ==
          [4, 8, 12, 24, 60]
    @test sprint(show, Cyclic(5)) == "Cyclic(5)"
    @test sprint(show, Dihedral(5)) == "Dihedral(5)"
    # A group whose smallest realization is a solid: a pyramid or prism needs at least 3 sides.
    @test_throws ArgumentError Polyhedron(Cyclic(2))
    @test_throws ArgumentError Polyhedron(Dihedral(2))
    @test nfaces(Polyhedron(Octahedral(); a=2.0)) == 6

    # Faces derived from the corners alone match the closed-surface requirement, and a
    # user-supplied solid needs nothing but its corners.
    p = Polyhedron(corners(Cube()))
    @test nfaces(p) == 6
    @test all(facedegree(p, i) == 4 for i in 1:6)
    @test length(rotationgroup(p)) == 24

    # A corner sitting mid-edge is not dropped: it lies on the planes of both faces meeting
    # there, so both gain a degree and the solid gains an edge. The shape is unchanged and the
    # dart count is not, which is what a combinatorial encoding of a subdivided face means.
    cs = collect(corners(Cube(2.0)))
    f = faces(Cube(2.0))[1]
    push!(cs, (cs[f[1]] + cs[f[2]]) / 2)
    sub = Polyhedron(cs)
    @test nfaces(sub) == 6
    @test sort([facedegree(sub, i) for i in 1:6]) == [4, 4, 4, 4, 5, 5]
    @test nedges(sub) == nedges(Cube()) + 1

    # An irregular but convex solid.
    box = Polyhedron([SVector(x, y, z) for x in (-1.0, 1.0) for y in (-2.0, 2.0) for z in (-3.0, 3.0)])
    @test nfaces(box) == 6
    @test length(rotationgroup(box)) == 4       # D_2
    @test symnum(dartencoding(box; labels=geometriclabels(box))[1]) == 4

    p = Cube(1.0f0)
    @test eltype(p) === Float32
    @test eltype(first(corners(p))) === Float32
    @test eltype(facecentroid(p, 1)) === Float32
    @test length(rotationgroup(p)) == 24
end

@testset "corner normalisation" begin
    # Corners are recentred on construction, and everything downstream may assume it:
    # `bounding_radius` and `inradius` measure from the origin and feed the overlap fast path,
    # so an off-centre solid would read both wrong with nothing to say so.
    cube = Cube()
    offset = SVector(3.0, -1.0, 7.0)
    shifted = Polyhedron([c + offset for c in corners(cube)], faces(cube))

    @test isapprox(sum(corners(shifted)) / ncorners(shifted), zero(SVector{3,Float64}); atol=1e-12)
    @test isapprox(bounding_radius(shifted), bounding_radius(cube); atol=1e-12)
    @test isapprox(inradius(shifted), inradius(cube); atol=1e-12)
    @test length(rotationgroup(shifted)) == length(rotationgroup(cube)) == 24
    # Same solid, so the same species down to the site frames.
    a = PolyhedronParticleSpecies(cube; colors=fill(1, 6))
    b = PolyhedronParticleSpecies(shifted; colors=fill(1, 6))
    @test symmetrynumber(a) == symmetrynumber(b) == 24
    for i in 1:6
        @test isapprox(Roly.bindingsites(a, i).pose.x, Roly.bindingsites(b, i).pose.x; atol=1e-12)
        @test isapprox(Roly.bindingsites(a, i).pose.psi, Roly.bindingsites(b, i).pose.psi; atol=1e-12)
    end

    # A corner used by no face contributes nothing to the shape while still counting towards
    # the bounding radius and dragging the centroid, so it is refused rather than carried.
    @test_throws ArgumentError Polyhedron([corners(cube); [SVector(0.0, 0.0, 0.0)]], faces(cube))
    # And the derived-faces path cannot hide one either: an interior point lies on no
    # supporting plane, so it lands in no face and the same check catches it.
    @test_throws ArgumentError Polyhedron([corners(cube); [SVector(0.1, 0.05, 0.0)]])
    # A point outside is a different error: it breaks convexity rather than going unused.
    @test_throws ArgumentError Polyhedron([corners(cube); [SVector(9.0, 0.0, 0.0)]], faces(cube))
end
