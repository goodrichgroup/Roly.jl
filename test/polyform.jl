using Roly: Polyform, nparticles, nsites, assemblysystem, symmetrynumber, dimension,
            bonds, bondindex, composition, interior_edges, exterior_edges, tocanon, toorig,
            AssemblySystem, UnitSquare, nbonds, raise!, lower!, bindingsites,
            collect_open_bindingsites, collect_compatible_pairs

@testset "polyform" begin
    sys = AssemblySystem([1 1 1 3; 1 2 1 4], UnitSquare)

    empty_poly = Polyform(sys)
    @test nparticles(empty_poly) == 0
    @test nsites(empty_poly) == 0

    mono = Polyform(sys, 1)
    @test nparticles(mono) == 1
    @test nsites(mono) == 4
    @test dimension(mono) == 2
    @test assemblysystem(mono) === sys
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
end
