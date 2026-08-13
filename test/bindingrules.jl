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
        # 3D does not opt in at all. Cubes tile space and a cube version was written and
        # measured: 1.9% on polycubes, where nauty dominates, against the subtlest part of the
        # check to get right. Prisms tile too, and would additionally need the bond table, since
        # a cap-to-side bond does not carry a cell to a cell.
        ("polycubes",           BindingRules([1 1 1 1], faced(Cube(), 1:6))),
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
end

@testset "attachment orbit representatives" begin
    using Roly: compatible_sitelocs, attachment_reps, siteloc2color, collect_compatible_pairs,
                raise!, Polyform, graphrep, PolyhedronParticleSpecies, PolygonParticleSpecies,
                Cube, Prism, nfaces, color, bindingsites
    using NautyGraphs: NautyDiGraph

    # A cube with every face alike has one orbit of six sites, so six compatible mates collapse
    # to one representative. With every face distinct there is nothing to collapse.
    alike = BindingRules([1 1 1 1], PolyhedronParticleSpecies(Cube(); colors=fill(1, 6)))
    @test length(compatible_sitelocs(alike, 1)) == 6
    @test length(attachment_reps(alike, 1)) == 1

    distinct = BindingRules(reduce(vcat, [[1 i 1 j] for i in 1:6 for j in i:6]),
                            PolyhedronParticleSpecies(Cube()))
    for c in 1:Roly.ncolors(distinct)
        @test length(attachment_reps(distinct, c)) == length(compatible_sitelocs(distinct, c))
    end

    # The claim itself, checked rather than argued, on both sides of the bond: the children
    # reachable through *every* compatible mate site and from *every* open host site are exactly
    # those reachable from the representatives. Anything dropped is a repeat, so the two sets of
    # canonical forms must be equal.
    #
    # The host side is here because getting it wrong is silent and expensive. nauty's orbit array
    # is indexed by the numbering it was handed, and `canonize=true` then permutes the graph, so
    # reading the partition against the graph afterwards merges unrelated vertices — which merged
    # open sites and quietly lost most of the polycubes. See `orbitreps!`.
    function children(poly, sitelocs_of; allsites=true)
        sys = Roly.bindingrules(poly)
        out = Set{NautyDiGraph}()
        usedorbits = Set{Int}()
        for orig_v in Roly.canonical_vertices(poly)
            part = Roly.particle_from_leadingvertex(poly, orig_v)
            isnothing(part) && continue
            for k in 1:Roly.nsites(part, sys)
                site = bindingsites(part, sys, k)
                Roly._isbound_vertex(poly, part, first(site.vertices)) && continue
                Roly.isinert(sys, color(site)) && continue
                if !allsites
                    rep = poly.orbits[Roly.tocanon(poly, first(site.vertices))]
                    rep in usedorbits && continue
                    push!(usedorbits, rep)
                end
                for siteloc in sitelocs_of(sys, color(site))
                    mate = bindingsites(Roly.species(sys, siteloc[1]), siteloc[2])
                    for r in 0:(Roly.nregistrations(site, mate) - 1)
                        trial = copy(poly)
                        ismissing(raise!(trial, site, siteloc, r)) && continue
                        push!(out, copy(graphrep(trial)))
                    end
                end
            end
        end
        return out
    end

    sides(p) = [i for i in 1:nfaces(p) if abs(Roly.facenormal(p, i)[3]) < 1e-8]
    p3 = Prism(3, 1.0; h=2.0)
    systems = [
        ("cubes",     alike),
        ("prisms",    BindingRules([1 first(sides(p3)) 1 first(sides(p3))],
                          PolyhedronParticleSpecies(p3;
                              colors=[i in sides(p3) ? 1 : 2 for i in 1:nfaces(p3)]))),
        ("hexagons",  BindingRules([1 1 1 1], PolygonParticleSpecies(6, 1.0; colors=fill(1, 6)))),
        ("triangles", BindingRules([1 1 1 1], PolygonParticleSpecies(3, 1.0; colors=fill(1, 3)))),
    ]
    for (name, sys) in systems
        poly = Polyform(sys, 1)
        for _ in 1:3     # monomer, then grow, so the host is asymmetric in later rounds too
            every = children(poly, compatible_sitelocs)
            @test every == children(poly, attachment_reps)
            @test every == children(poly, attachment_reps; allsites=false)
            nxt = nothing
            for (site, loc, r) in collect_compatible_pairs(poly)
                trial = copy(poly)
                ismissing(raise!(trial, site, loc, r)) && continue
                nxt = trial; break
            end
            isnothing(nxt) && break
            poly = nxt
        end
    end
end
