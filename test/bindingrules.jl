using Roly: nspecies, nbonds, nsites, dimension, species,
            interactionmatrix, bonded_colors, bonded_sites, bonded_species,
            siteloc2color, color2siteloc, color2species, isinert, compatible_sitelocs

@testset "bindingrules" begin
    sys = BindingRules([1 1 1 3; 1 2 1 4], UnitSquare)

    @test nspecies(sys) == 1
    @test nbonds(sys) == 2
    @test nsites(sys) == 4
    @test dimension(sys) == 2
    @test length(species(sys)) == 1
    @test species(sys, 1) === species(sys)[1]

    imat = interactionmatrix(sys)
    @test size(imat) == (4, 4)
    @test issymmetric(imat)

    bc = bonded_colors(sys)
    @test length(bc) == 2
    @test all(c isa NTuple{2,Int} for c in bc)

    bs = bonded_sites(sys)
    @test length(bs) == 2

    bsp = bonded_species(sys)
    @test length(bsp) == 2
    @test all(==((1, 1)), bsp)

    c1 = siteloc2color(sys, (1, 1))
    c3 = siteloc2color(sys, (1, 3))
    @test imat[c1, c3]
    @test color2siteloc(sys, c1) == [(1, 1)]
    @test color2species(sys, c1) == 1

    @test !isinert(sys, c1)
    @test !isinert(sys, (1, 1))

    @test compatible_sitelocs(sys, c1) == [(1, 3)]
    @test compatible_sitelocs(sys, c3) == [(1, 1)]

    sys1bond = BindingRules([1 1 1 3], UnitSquare)
    c2 = siteloc2color(sys1bond, (1, 2))
    @test isinert(sys1bond, c2)
    @test isinert(sys1bond, (1, 2))

    io = IOBuffer()
    show(io, sys)
    @test contains(String(take!(io)), "BindingRules")

    sys1 = BindingRules([1 1 2 1], UnitTriangle)
    sys2 = BindingRules([1 3 2 2], UnitTriangle)
    sys3 = BindingRules([1 3 3 2], UnitTriangle)

    g1 = graphrep(sys1)
    g2 = graphrep(sys2)
    g3 = graphrep(sys3)

    @test g1 ≃ g2
    @test !(g1 ≃ g3)
end

@testset "on-lattice detection" begin
    using Roly: PolygonParticleSpecies, PolyhedronParticleSpecies, PatchyDisk, Cube, Prism,
                Tetrahedron, Polyhedron, facenormal, nfaces
    using StaticArrays: SVector

    poly(n, a=1.0) = PolygonParticleSpecies(n, a; colors=fill(1, n))
    sides(p) = [i for i in 1:nfaces(p) if abs(facenormal(p, i)[3]) < 1e-8]
    faced(p, sticky; kw...) =
        PolyhedronParticleSpecies(p; colors=[i in sticky ? 1 : 2 for i in 1:nfaces(p)], kw...)

    # A system is on-lattice when every species tiles space at one size and no bond can leave
    # the tiling. Then two particles overlap only by occupying the same cell, and `overlap`
    # answers with a subtraction. Anything unrecognised just takes the real geometry.
    for (name, sys) in [
        ("squares",             BindingRules([1 1 1 1], poly(4))),
        ("triangles",           BindingRules([1 1 1 1], poly(3))),
        ("hexagons",            BindingRules([1 1 1 1], poly(6))),
        ("squares, 2 species",  BindingRules([1 1 2 1], [poly(4), poly(4)])),
        ("polycubes",           BindingRules([1 1 1 1], faced(Cube(), 1:6))),
        # Prism(4) at its default height *is* a cube, which a check on the constructor's name
        # would miss and a check on the geometry does not.
        ("cubes as Prism(4)",   BindingRules([1 1 1 1], faced(Prism(4), 1:6))),
        # A whole dart step turns a cube's site onto another of its own edges, so the partner
        # still arrives lattice-aligned.
        ("cubes, dart twist",   BindingRules([1 1 1 1], faced(Cube(), 1:6; twists=π/2))),
    ]
        @test sys._onlattice
    end

    for (name, sys) in [
        # Pentagons do not tile the plane.
        ("pentagons",           BindingRules([1 1 1 1], poly(5))),
        # Two tilings at different scales share no cells.
        ("squares, two sizes",  BindingRules([1 1 2 1], [poly(4), poly(4, 2.0)])),
        # Same size, but a square and a triangle bond into no one tiling.
        ("squares + triangles", BindingRules([1 1 2 1], [poly(4), poly(3)])),
        # A fraction of a dart step hands the partner over turned off the lattice, which is
        # exactly the case the site check exists to catch.
        ("cubes, 0.37 twist",   BindingRules([1 1 1 1], faced(Cube(), 1:6; twists=0.37))),
        # Prisms tile, but a bond between a cap and a side does not carry a cell to a cell, and
        # the species cannot tell which bonds the rules will allow.
        ("square prisms",       BindingRules([1 1 1 1], faced(Prism(4, 1.0; h=2.0),
                                                             sides(Prism(4, 1.0; h=2.0))))),
        ("tetrahedra",          BindingRules([1 1 1 1],
                                    PolyhedronParticleSpecies(Tetrahedron(); colors=fill(1, 4)))),
        ("patchy disks",        BindingRules([1 1 1 1],
                                    PatchyDisk([0.0, 2π/3, 4π/3]; colors=fill(1, 3)))),
    ]
        @test !sys._onlattice
    end

    # The point of it: the shortcut has to agree with the real geometry everywhere it is taken.
    # Enumerate each lattice both ways and require identical counts -- a single wrong overlap
    # answer would change them -- and require the OEIS sequence on top, so that agreeing on a
    # wrong answer is ruled out too. These are the *one-sided* counts: a 2D assembly cannot be
    # turned over, so mirror images stay distinct, unlike the 3D prisms in `test/enumeration.jl`
    # that tile the same lattices.
    for (n, want) in ((3, [1, 1, 1, 4, 6, 19, 43, 120]),    # https://oeis.org/A006534
                      (4, [1, 1, 2, 7, 18, 60, 196]),       # https://oeis.org/A000988
                      (6, [1, 1, 3, 10, 33, 147]))          # https://oeis.org/A006535
        sys = BindingRules([1 1 1 1], poly(n))
        @test sys._onlattice
        counts = [polyenum(sys; maxsize=i)[1] for i in eachindex(want)]
        @test counts == cumsum(want)
        slow = withoutlattice(sys)
        @test !slow._onlattice
        @test [polyenum(slow; maxsize=i)[1] for i in eachindex(want)] == counts
    end

    # And in 3D, where the shortcut replaces a separating-axis test rather than a cheap one.
    cubes = BindingRules([1 1 1 1], faced(Cube(), 1:6))
    polycubes = [1, 1, 2, 8, 29, 166]                        # https://oeis.org/A000162
    @test [polyenum(cubes; maxsize=i)[1] for i in eachindex(polycubes)] == cumsum(polycubes)
    @test [polyenum(withoutlattice(cubes); maxsize=i)[1] for i in eachindex(polycubes)] ==
          cumsum(polycubes)
end
