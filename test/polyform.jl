using Roly: Polyform, nparticles, nsites, bindingrules, symmetrynumber, dimension,
            bonds, bondindex, composition, interior_edges, exterior_edges, tocanon, toorig,
            BindingRules, UnitSquare, nbonds, raise!, lower!, bindingsites, graphrep,
            collect_open_bindingsites, collect_compatible_pairs, particle_from_leadingvertex,
            PolygonParticleSpecies, species, polygen

@testset "polyform" begin
    sys = BindingRules([1 1 1 3; 1 2 1 4], UnitSquare)

    empty_poly = Polyform(sys)
    @test nparticles(empty_poly) == 0
    @test nsites(empty_poly) == 0

    mono = Polyform(sys, 1)
    @test nparticles(mono) == 1
    @test nsites(mono) == 4
    @test dimension(mono) == 2
    @test bindingrules(mono) === sys
    @test symmetrynumber(mono) == 1

    @test mono == Polyform(sys, 1)
    @test mono != empty_poly

    mono2 = copy(mono)
    @test mono == mono2
    copy!(empty_poly, mono)
    @test mono == empty_poly

    io = IOBuffer()
    show(io, mono)
    @test contains(String(take!(io)), "Polyform")

    @test !isempty(collect(interior_edges(mono)))
    @test isempty(collect(exterior_edges(mono)))
    @test isempty(collect(bonds(mono)))

    comp = composition(mono)
    @test comp[1] == 1
    @test sum(comp[2:end]) == 0

    @test length(collect(bindingsites(mono))) == nsites(mono)

    open_sites = collect_open_bindingsites(mono)
    @test length(open_sites) == 4

    for v in eachindex(mono.canon2orig)
        @test tocanon(mono, toorig(mono, v)) == v
    end

    di = copy(mono)
    site, siteloc = first(collect_compatible_pairs(di))
    @test !isnothing(raise!(di, site, siteloc))
    @test nparticles(di) == 2
    @test nsites(di) == 8

    @test !isempty(collect(exterior_edges(di)))
    @test length(collect(bonds(di))) == 1

    for e in exterior_edges(di)
        idx = bondindex(di, e.src, e.dst)
        @test !isnothing(idx)
        @test 1 <= idx <= nbonds(sys)
    end

    comp_di = composition(di)
    @test comp_di[1] == 2
    @test sum(comp_di[2:end]) == 1

    lower!(di)
    @test nparticles(di) == 1
    lower!(di)
    @test nparticles(di) == 0
    @test isnothing(lower!(di))

    ### Check raise / remove equality
    sys = BindingRules([1 1 1 4; 1 4 2 2; 1 1 3 2; 1 1 3 4; 2 1 3 3; 2 1 3 1], UnitSquare)
    poly = polygen(sys)[end]

    poly_raise = copy(poly)
    lower!(poly_raise)

    compatible_pairs = collect_compatible_pairs(poly_raise)
    site, siteloc = first(compatible_pairs)
    i = 1
    while ismissing(raise!(poly_raise, compatible_pairs[i]...))
        i += 1
        site_siteloc = compatible_pairs[i]
    end

    @test poly == poly_raise
    @test poly_raise.canon2orig == poly.canon2orig
    @test poly_raise.orig2canon == poly.orig2canon

    using Graphs: has_edge, nv
    using Roly: graphvertices, PolyhedronParticleSpecies, UnitCube, Cube

    # A monomer's canon2orig is the permutation `canonize!` applied to the species graph,
    # not the identity. Checking it via the edges: every edge of the species graph must
    # reappear between the corresponding canonical vertices.
    for spcs in (PolygonParticleSpecies(6; labels=[1, 2, 1, 2, 1, 2]),
                 PolygonParticleSpecies(4; labels=[1, 1, 1, 1]))
        s = BindingRules([1 1 1 3], spcs)
        mono = Polyform(s, 1)
        for e in edges(graphrep(species(s, 1)))
            @test has_edge(graphrep(mono), tocanon(mono, e.src), tocanon(mono, e.dst))
        end
        @test nv(graphrep(mono)) == nv(graphrep(species(s, 1)))
    end

    # Same for a 3D species, whose 24-vertex graph really does get permuted.
    s3 = BindingRules([1 1 1 2], PolyhedronParticleSpecies(Cube(); encoding=:dart))
    mono3 = Polyform(s3, 1)
    for e in edges(graphrep(species(s3, 1)))
        @test has_edge(graphrep(mono3), tocanon(mono3, e.src), tocanon(mono3, e.dst))
    end

    # Removing a particle compacts the graph's vertex numbering, so the surviving
    # particles' vertex blocks have to be compacted with it: they must still tile the
    # graph exactly, with nothing left hanging off the end.
    sys_pm = BindingRules([1 1 1 3; 1 2 1 4], PolygonParticleSpecies(4, 1.0; labels=[1, 1, 1, 1]))
    for p in polygen(sys_pm; maxsize=6)
        nparticles(p) < 2 && continue
        q = copy(p)
        lower!(q)
        blocks = sort(reduce(vcat, [collect(graphvertices(pt, sys_pm)) for pt in q.particles]))
        @test blocks == 1:nv(graphrep(q))
    end
end
