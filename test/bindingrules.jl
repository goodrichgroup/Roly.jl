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

    # Constructor from an interaction matrix (round-trip via intmat).
    ref = BindingRules([1 1 1 3; 1 2 1 4], UnitSquare)
    im  = Matrix{Bool}(interactionmatrix(ref))
    reconstructed = BindingRules(im, UnitSquare)
    @test interactionmatrix(reconstructed) == interactionmatrix(ref)
    @test nbonds(reconstructed) == nbonds(ref)
    @test bonded_colors(reconstructed) == bonded_colors(ref)

    # Multi-species intmat round-trip.
    ref2 = BindingRules([1 1 2 1], UnitTriangle)  # 2 species x 3 sites = 6 colors
    im2  = Matrix{Bool}(interactionmatrix(ref2))
    reconstructed2 = BindingRules(im2, [copy(UnitTriangle) for _ in 1:nspecies(ref2)])
    @test interactionmatrix(reconstructed2) == interactionmatrix(ref2)
    @test bonded_colors(reconstructed2) == bonded_colors(ref2)

    # Error cases.
    @test_throws ArgumentError BindingRules(falses(3, 4), UnitSquare)          # not square
    @test_throws ArgumentError BindingRules(falses(3, 3), UnitSquare)          # wrong ncolors
    bad_asymmetric = Bool[false true false false;
                          false false false false;
                          false false false false;
                          false false false false]
    @test_throws ArgumentError BindingRules(bad_asymmetric, UnitSquare)        # not symmetric
end
