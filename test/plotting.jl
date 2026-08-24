@testset "plotting" begin
    # every species must have a `plot_particlespecies!`
    # method reachable through `render`, and a polyform of them must draw as well.
    ext = Base.get_extension(Roly, :MakieExt)
    @test !isnothing(ext)

    tri = PolygonParticleSpecies(3)
    disk = PatchyDisk([0.0, 2π / 3, 4π / 3])
    # Distinct colors per face, so a rules table naming one bond leaves the other four inert
    # and the inert-site rendering below has something to grey out.
    cube = PolyhedronParticleSpecies(Cube())
    sphere = PatchySphere(Cube(), 1.0)

    for spcs in (tri, disk, cube, sphere, UnitTetrahedron, UnitDodecahedron, UnitPrism(5))
        @test render(spcs) isa Figure
    end

    # A 3D polyform, drawn particle by particle at its bonded pose.
    opposite(p, i) = findfirst(
        j -> isapprox(dot(facenormal(p, i), facenormal(p, j)), -1; atol=1e-8), 1:nfaces(p)
    )
    pairs = unique([minmax(i, opposite(Cube(), i)) for i in 1:6])
    rules = BindingRules(reduce(vcat, [[1 a 1 b] for (a, b) in pairs]), cube)
    polys = polygen(rules; maxsize=3)
    @test render(polys[end]) isa Figure

    # A 2D polyform, so the shared color resolution is covered in both dimensions.
    sys2d = BindingRules([1 1 1 1], tri)
    polys2d = polygen(sys2d; maxsize=3)
    @test render(polys2d[end]) isa Figure

    # An explicitly passed `speciesindex` picks the palette instead of being overridden by
    # the default of 1.
    _, _, colors1, bonding1 = ext._resolve_colors(cube, 1, nothing, nothing)
    _, _, colors2, _ = ext._resolve_colors(cube, 2, nothing, nothing)
    @test length(colors1) == length(colors2) == nsites(cube)
    @test colors1 != colors2
    # With no rules in scope every site counts as bonding, so nothing is tinted away.
    @test all(bonding1)

    # By default inert site are greyed out, which is what markers do.
    rules = BindingRules([1 1 1 6], cube)
    _, pal, colors, bonding = ext._resolve_colors(cube, 1, rules, nothing)
    @test count(!, bonding) == 4
    for i in 1:nsites(cube)
        @test (colors[i] == ext.INERT_COLOR) == !bonding[i]
    end
    # A polyhedron keeps the species hue throughout and separates by tint instead.
    _, pal, colors, _ = ext._resolve_colors(
        cube, 1, rules, nothing; bond_tint=ext.FACE_TINT, inert_color=nothing
    )
    for i in 1:nsites(cube)
        want = ext._tint(pal[i], bonding[i] ? ext.FACE_TINT : ext.BODY_TINT)
        @test colors[i] == want
    end
    # Tinting is skipped entirely when a `site_color` callback is given.
    _, _, colors, _ = ext._resolve_colors(cube, 1, rules, (_, i) -> ext.INERT_COLOR)
    @test all(==(ext.INERT_COLOR), colors)
end
