@testset "plotting" begin
    # Smoke tests for the Makie extension: every species must have a `plot_particlespecies!`
    # method reachable through `render`, and a polyform of them must draw as well. CairoMakie
    # is only a stand-in here — it has no depth buffer, so 3D output is not to be trusted
    # visually, but it exercises the same code path as GLMakie.
    tri = PolygonParticleSpecies(3)
    disk = PatchyDisk([0.0, 2π / 3, 4π / 3])
    cube = PolyhedronParticleSpecies(Cube(); labels=fill(1, 6))
    sphere = PatchySphere(Cube(), 1.0)

    for spcs in (tri, disk, cube, sphere, UnitTetrahedron, UnitDodecahedron, UnitPrism(5))
        @test render(spcs) isa Makie.Figure
    end

    # A 3D polyform, drawn particle by particle with its bonded poses.
    opposite(p, i) = findfirst(
        j -> isapprox(dot(Roly.facenormal(p, i), Roly.facenormal(p, j)), -1; atol=1e-8), 1:nfaces(p)
    )
    pairs = unique([minmax(i, opposite(Cube(), i)) for i in 1:6])
    sys = BindingRules(reduce(vcat, [[1 a 1 b] for (a, b) in pairs]), cube)
    polys, _ = polyenum(sys; max_size=4)
    @test render(polys[end]) isa Makie.Figure

    # A 2D polyform, so the shared color resolution is covered in both dimensions.
    sys2d = BindingRules([1 1 1 1], tri)
    polys2d, _ = polyenum(sys2d; max_size=3)
    @test render(polys2d[end]) isa Makie.Figure

    # `species_index` selects the palette; passing it explicitly must not be overridden.
    _, pal2, _ = Roly.MakieExt._resolve_colors(cube, 2, nothing, nothing)
    @test pal2 == Roly.MakieExt.species_palette(2, nsites(cube))
end
