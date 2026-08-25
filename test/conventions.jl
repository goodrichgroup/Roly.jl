@testset "API conventions" begin
    using Roly: SpeciesSite, ParticleSite

    exported = [n for n in names(Roly) if Base.isexported(Roly, n)]
    # names Roly defines, as opposed to the handful it re-exports from Rotations and StaticArrays
    ours = filter(exported) do n
        v = getproperty(Roly, n)
        return (v isa Function || v isa Type) && parentmodule(v) === Roly
    end
    @test length(ours) > 60

    # A name reaching every user's namespace has to be one word. Underscores are for internals,
    # where the two halves are usually a qualifier and a noun that only make sense together.
    @test isempty(filter(n -> occursin('_', string(n)) && !endswith(string(n), "!"), ours))

    # Nothing exported takes or yields a bare graph vertex. Which numbering a vertex is in,
    # canonical or original, is an invariant callers inside the package keep; it must not be
    # something a user of the package has to know about.
    for n in (:bondindex, :interior_edges, :exterior_edges, :subpolyform)
        @test Base.ispublic(Roly, n)
        @test !Base.isexported(Roly, n)
    end

    # The functions that do take one demand to be told which numbering, with no default to fall
    # into. `UndefKeywordError` rather than a wrong answer.
    rules = BindingRules([1 1 1 3; 1 2 1 4], UnitSquare)
    dimer = first(p for p in polygen(rules; maxsize=2) if nparticles(p) == 2)
    part = dimer.particles[1]
    @test_throws UndefKeywordError Roly._same_particle(dimer, 1, 2)
    @test_throws UndefKeywordError Roly._isbound_vertex(dimer, part, 1)
    @test_throws UndefKeywordError Roly._vertex_to_particle_site(dimer, 1)

    # The two site addresses are different types, so neither can stand in for the other.
    @test SpeciesSite(1, 2) != ParticleSite(1, 2)
    @test !isa(SpeciesSite(1, 2), ParticleSite)
    @test_throws MethodError color(rules, ParticleSite(1, 1))
    @test_throws MethodError raise!(copy(dimer), bindingsite(dimer, first(opensites(dimer))), ParticleSite(1, 1))
    # but each still destructures like the pair it replaced
    @test (SpeciesSite(3, 4)...,) == (3, 4)
    @test (ParticleSite(3, 4)...,) == (3, 4)

    # each address indexes the thing it names
    @test bindingsite(dimer, ParticleSite(1, 1)) == bindingsite(dimer.particles[1], rules, 1)
    @test bindingsite(rules, SpeciesSite(1, 3)) == bindingsite(Roly.species(rules, 1), 3)
    @test_throws MethodError bindingsite(dimer, SpeciesSite(1, 1))
    @test_throws MethodError bindingsite(rules, ParticleSite(1, 1))

    # `bonds` speaks ParticleSite, the rules speak SpeciesSite
    @test eltype(collect(bonds(dimer))) == Pair{ParticleSite,ParticleSite}
    # the site accessors return addresses; the site itself is one index away
    @test eltype(opensites(dimer)) == ParticleSite
    @test eltype(exposedsites(dimer)) == ParticleSite
    @test opensites(dimer) ⊆ exposedsites(dimer)
    @test all(l -> bindingsite(dimer, l) isa BindingSite, exposedsites(dimer))
    @test eltype(Roly.possible_attachments(rules, 1)) == SpeciesSite
end
