using Roly: BindingSite, shift_vertices, shift_color, isaligned, istouching, isincontact, color,
            standard_offset, contact_pairing
using Roly: BindingSite, twistfreedom, bondperiod, nphases, phase,
            standard_offset, isaligned
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

    # A bond keeps the symmetry the two sites have in common. 
    # Joining a k1-fold site to a k2-fold one must leave gcd(k1, k2) turns, since a
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

    ### contact_pairing under a phase
    # Vertex a of site 1 sits at angle 2πa/k1 about the bond axis; vertex b of site 2 sits at
    # 2πr/L - 2πb/k2, its cyclic order reversed by the gluing and its frame turned by the
    # phase. The pairing is exactly the coincidences, so check it against them directly.
    function coincidences(k1, k2, L, r)
        M = lcm(lcm(k1, k2), L)
        return sort([(a, b) for a in 0:(k1 - 1), b in 0:(k2 - 1)
                     if M * a ÷ k1 == mod(-(M * b ÷ k2) + M * r ÷ L, M)][:])
    end
    pairs0(k1, k2, r, L) = sort([(a - 1, b - k1 - 1)
                                 for (a, b) in contact_pairing(1:k1, (k1 + 1):(k1 + k2), r, L)])

    divisors(n) = [d for d in 1:n if n % d == 0]
    for k1 in 1:8, k2 in 1:8, q1 in divisors(k1), q2 in divisors(k2)
        L = lcm(q1, q2)
        # A twist freedom divides its site's vertex count, so L always divides lcm(k1, k2) and
        # the phase is expressible. This is what makes the offset an integer.
        @test lcm(k1, k2) % L == 0
        seen = Set()
        for r in 0:(L - 1)
            p = pairs0(k1, k2, r, L)
            @test p == coincidences(k1, k2, L, r)
            # The count never depends on the phase, so a bond keeps its residual
            # symmetry however it is turned
            @test length(p) == gcd(k1, k2)
            # Distinct phases give distinct pairings, so the graph records which
            # one a bond is in. Without that, geometrically different assemblies would share a
            # canonical form and be silently merged.
            @test p ∉ seen
            push!(seen, p)
        end
    end
    # r = 0 is the base phase and must reproduce the plain call exactly.
    for (k1, k2) in [(1, 1), (3, 3), (4, 4), (6, 3), (4, 6), (12, 8), (4, 5)]
        @test collect(contact_pairing(1:k1, 11:(10 + k2), 0, 4)) ==
              collect(contact_pairing(1:k1, 11:(10 + k2)))
    end
    # A phase a site's vertex count cannot express is refused rather than truncated.
    @test_throws ArgumentError collect(contact_pairing(1:1, 2:2, 1, 4))

    # phases
    site(gauge, stab, locking) =
        BindingSite(Pose{3,Float64,RotMatrix3{Float64}}(SVector(1.0, 0.0, 0.0),
                                                        one(RotMatrix3{Float64})),
                    1, 1:gauge, 1e-8, 1e-8, gauge, stab, locking)

    # A locking site's restricts down to its stabilizer; a rotation-free one restricts to its gauge
    # The two coincide when a face is no more symmetric than the body around it.
    @test twistfreedom(site(4, 2, true)) == 2
    @test twistfreedom(site(4, 2, false)) == 4
    @test twistfreedom(site(4, 4, true)) == twistfreedom(site(4, 4, false)) == 4

    # Phases come from symmetry: an unsymmetric particle has exactly one per bond,
    # however symmetric the face it bonds through.
    @test nphases(site(4, 1, true), site(4, 1, true)) == 1
    @test nphases(site(1, 1, true), site(1, 1, true)) == 1
    # A cube face onto a cube face, both 4-fold: still one, the four turns being symmetries.
    @test nphases(site(4, 4, true), site(4, 4, true)) == 1
    # A cube (4-fold) meeting a prism's square side face (2-fold): in-plane and out-of-plane.
    @test nphases(site(4, 4, true), site(4, 2, true)) == 2
    # ...and only one the other way round, since a lone cube turned about the bond axis maps
    # onto itself, so there is only one distinct dimer to find.
    @test nphases(site(4, 2, true), site(4, 4, true)) == 1
    # Freeing one side opens the bond up whatever the other says.
    @test nphases(site(4, 2, false), site(4, 2, true)) == 2
    @test nphases(site(4, 2, true), site(4, 2, true)) == 1

    # bondperiod is the lcm because the admissible set must be closed under both sites' turns.
    @test bondperiod(site(6, 6, true), site(4, 4, true)) == 12
    @test bondperiod(site(4, 2, true), site(6, 3, true)) == 6

    # A partner placed in phase r is recognised as being in phase r, and reading
    # the pair the other way round gives the same answer
    for L in (1, 2, 3, 4, 6)
        b1 = site(L, L, true)
        for r in 0:(L - 1)
            b2 = BindingSite(standard_offset(b1, r, L), 1, 1:L, 1e-8, 1e-8, L, L, true)
            @test phase(b1, b2) == r
            @test phase(b2, b1) == r
            @test isaligned(b1, b2)
        end
    end
    
    # An orientation in no phase at all is not a bond.
    b = site(2, 2, true)
    off = BindingSite(standard_offset(b, 1, 4), 1, 1:2, 1e-8, 1e-8, 2, 2, true)
    @test isnothing(phase(b, off))
    @test !isaligned(b, off)
end
