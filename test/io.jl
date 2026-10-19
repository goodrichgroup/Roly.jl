@testset "io" begin
    samesite(a, b) = color(a) == color(b) && a.vertices == b.vertices && a.locking == b.locking &&
                     a.sitesym == b.sitesym && a.stab == b.stab && isapprox(a.pose.x, b.pose.x; atol=1e-10) &&
                     isapprox(Matrix(a.pose.psi), Matrix(b.pose.psi); atol=1e-10)
    samespecies(a, b) = typeof(a) == typeof(b) && graphrep(a) == graphrep(b) && nsites(a) == nsites(b) &&
                        all(i -> samesite(bindingsite(a, i), bindingsite(b, i)), 1:nsites(a))
    samerules(a, b) = interactionmatrix(a) == interactionmatrix(b) && nspecies(a) == nspecies(b) &&
                      all(samespecies(x, y) for (x, y) in zip(species(a), species(b)))
    sameparticles(a, b) = nparticles(a) == nparticles(b) && all(zip(a.particles, b.particles)) do (p, q)
        p.speciesindex == q.speciesindex && isapprox(p.pose.x, q.pose.x; atol=1e-9) &&
            isapprox(Matrix(p.pose.psi), Matrix(q.pose.psi); atol=1e-9)
    end
    samepoly(a, b) = graphrep(a) == graphrep(b) && sameparticles(a, b)

    sq = BindingRules([1 1 1 1], SymmetricUnitSquare)
    allrules = [
        "squares" => BindingRules([1 1 2 3; 2 1 3 3], UnitSquare),
        "rectangles" => BindingRules([1 1 1 3; 1 2 1 4], PolygonParticleSpecies([[1, 0.5], [-1, 0.5], [-1, -0.5], [1, -0.5]])),
        "cubes" => BindingRules([1 1 1 2], UnitCube),
        "prisms" => BindingRules([1 1 1 2], PolyhedronParticleSpecies(Prism(3); twists=0.3)),
        "patchy disks" => BindingRules([1 1 1 2], PatchyDisk([0, π / 2, π, 3π / 2])),
        "patchy spheres" => BindingRules([1 1 1 2], PatchySphere(Cube(); twists=0.2)),
        "meta" => BindingRules([1 1 1 1], MetaParticleSpecies(last(polygen(sq; maxsize=3)))),
    ]

    dir = mktempdir()
    for (name, rules) in allrules
        @testset "$name" begin
            polys = filter(p -> nparticles(p) > 0, polygen(rules; maxsize=3))

            # through the dictionaries, and through files
            @test samerules(Roly.fromdict(BindingRules, Roly.todict(rules)), rules)
            @test all(samepoly(Roly.fromdict(Polyform, Roly.todict(p); rules), p) for p in polys)
            @test samepoly(Roly.fromdict(Polyform, Roly.todict(last(polys); rules=true)), last(polys))

            @test samerules(fromfile(tofile(joinpath(dir, "rules.rly"), rules)), rules)
            p = fromfile(tofile(joinpath(dir, "poly.rly"), last(polys)))
            @test samepoly(p, last(polys))
            @test samerules(bindingrules(p), rules)
            back = fromfile(tofile(joinpath(dir, "polys.toml"), polys))
            @test !occursin("polyforms.rules", read(joinpath(dir, "polys.toml"), String))   # rules only once
            @test length(back) == length(polys)
            @test all(samepoly(a, b) for (a, b) in zip(back, polys))
            @test all(q -> bindingrules(q) === bindingrules(first(back)), back)

            # and through JLD2, whose files are plain HDF5
            @test samerules(fromfile(tofile(joinpath(dir, "rules.rlyb"), rules)), rules)
            @test samepoly(fromfile(tofile(joinpath(dir, "poly.jld2"), last(polys))), last(polys))
            back = fromfile(tofile(joinpath(dir, "polys.rlyb"), polys))
            @test length(back) == length(polys)
            @test all(samepoly(a, b) for (a, b) in zip(back, polys))
        end
    end

    # particles keep their order even when it cannot be attached in sequence
    chain = BindingRules([1 1 2 3; 2 1 3 3], UnitSquare)
    line = only(filter(p -> nparticles(p) == 3, polygen(chain; maxsize=3)))
    ends = [i for i in 1:3 if count(b -> i in (first(b).particle, last(b).particle), bonds(line)) == 1]
    middle = only(setdiff(1:3, ends))
    shuffled = Roly.subpolyform(line, [ends; middle])
    @test samepoly(Roly.fromdict(Polyform, Roly.todict(shuffled); rules=chain), shuffled)

    # the file header
    text = read(tofile(joinpath(dir, "chain.rly"), line), String)
    @test startswith(text, "contenttype = \"roly/polysingle\"\nformatversion = 1\nrolyversion = ")
    other = joinpath(dir, "other.toml")
    write(other, "a = 1\n")
    @test_throws ArgumentError fromfile(other)
    @test_throws ArgumentError tofile(joinpath(dir, "empty.rly"), Polyform[])
    @test_throws ArgumentError tofile(joinpath(dir, "x.rly"), chain; format=:yaml)
    @test_throws MethodError tofile(joinpath(dir, "x.rly"), 3)

    # a species written without a reader cannot be read back
    @test_throws ArgumentError Roly.fromdict(ParticleSpecies, Dict("type" => "Roly.NoSuchSpecies"))
    @test_throws ArgumentError Roly.fromdict(ParticleSpecies, Dict("type" => "NoSuchPackage.Species"))
    @test_throws ArgumentError Roly.fromdict(ParticleSpecies, Dict("type" => "Roly.BindingRules"))
    @test Roly.todict(UnitSquare)["type"] == "Roly.PolygonParticleSpecies"
end
