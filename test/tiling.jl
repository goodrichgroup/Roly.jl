@testset "tiling" begin
    # The copies of a cell at the integer combinations of its translations, each lattice point
    # once. The translations generate the lattice but need not be a basis of it, so two
    # combinations can land on one point.
    function patch(t; reps)
        cell = unitcell(t)
        vs = translations(t)
        seen = Set{NTuple{dimension(t),Float64}}()
        offsets = eltype(vs)[]
        for c in Iterators.product(ntuple(_ -> (-reps):reps, length(vs))...)
            shift = isempty(vs) ? zero(eltype(vs)) : sum(c[i] * vs[i] for i in eachindex(vs))
            key = Tuple(round.(shift; digits=7) .+ 0.0)
            key in seen && continue
            push!(seen, key)
            push!(offsets, shift)
        end
        center = cell.particles
        outer = [Roly.translate(p, s) for s in offsets for p in center if !all(iszero, s)]
        return center, outer
    end

    # the chain monomer is a unit cell of the infinite chain: one translation closes its 1-3 bond
    # onto its own copy
    chainlike = BindingRules([1 1 1 3], UnitSquare)
    chains = tilings(chainlike; maxsize=1)
    @test length(chains) == 1
    @test all(iscomplete(t) && bondtypes(t) == [1] && length(translations(t)) == 1 for t in chains)

    # the square-lattice monomer tiles the plane with two translations. Partial closures along a
    # single axis, the infinite chains, are enumerated alongside the complete tilings
    squarerules = BindingRules([1 2 1 4; 1 3 1 1], UnitSquare)
    ts = tilings(squarerules; maxsize=1)
    @test count(iscomplete, ts) > 0
    @test all(length(bondtypes(t)) == 2 && length(translations(t)) == 2 for t in ts if iscomplete(t))
    @test any(t -> !iscomplete(t) && length(bondtypes(t)) == 1, ts)

    # a closed dimer has no open sites left to close
    dimerrules = BindingRules([1 1 2 1], UnitSquare)
    @test isempty(tilings(dimerrules; maxsize=2))

    # a cell's sites are read off the cell. `BindingRules` shifts every species' colors into a
    # range of its own, so two square species whose open sites do not start at color 1:
    shifted = BindingRules([1 4 2 2; 2 1 1 3], UnitSquare)
    shiftedchain = polygen(shifted; maxsize=2)[end]
    @test [color(bindingsite(shiftedchain, l)) for l in opensitelocs(shiftedchain)] == [3, 5]
    @test (3, 5) in Roly.bonded_colors(shifted)
    @test length(tilings(shifted; maxsize=2)) == 1

    @testset "a Tiling answers for itself" begin
        # one entry per structure, not one per way the search reached it. A lattice does not care
        # which order its generators came in nor which way each points, and there are 2^d*d! ways
        # to walk the same one
        cube = BindingRules([1 1 1 1], PolyhedronParticleSpecies(Cube(); colors=fill(1, 6)))
        cts = tilings(cube; maxsize=1)
        @test count(iscomplete, cts) == 1
        @test !any(cts[i] == cts[j] for i in eachindex(cts) for j in 1:(i - 1))

        # the cell comes back as a polyform, which is what makes a tiling something to look at
        whole = first(filter(iscomplete, cts))
        @test unitcell(whole) isa Polyform
        @test nparticles(unitcell(whole)) == tilingorder(whole) == 1
        @test length(translations(whole)) == dimension(whole) == 3
        @test bondtypes(whole) == [1, 1, 1]
        @test bindingrules(unitcell(whole)) === cube
        @test nbonds(whole) == 3 && nbonds(unitcell(whole)) == 0
        @test symmetrynumber(whole) >= 1

        # site 1 binds site 2, a quarter turn, so no translation carries a single square onto a
        # bonded copy: the cell has to hold two of them, related by that turn
        turn = BindingRules([1 1 1 2; 1 3 1 4], UnitSquare)
        @test isempty(tilings(turn; maxsize=1))
        turned = tilings(turn; maxsize=2)
        @test !isempty(turned)
        @test all(tilingorder(t) == 2 for t in turned)
        pair = first(turned)
        @test nparticles(unitcell(pair)) == 2
        @test nbonds(unitcell(pair)) == 1
        @test bindingrules(unitcell(pair)) === turn
        # the bond joining the two copies belongs to the cell, so it is counted there
        @test !isempty(bondtypes(pair))

        # a partial closure runs along fewer directions than the space has
        strip = first(filter(t -> !iscomplete(t), cts))
        @test length(translations(strip)) < dimension(strip)

        # a cell that is several copies of a smaller one is never reported: raising the bound adds
        # no repeats of what is already there
        for k in 1:3
            @test length(tilings(chainlike; maxsize=k)) == 1
        end
    end

    @testset "3D" begin
        using Roly: nfaces, facenormal, twistfreedom

        # every face binds every face, so the cube monomer is a unit cell of space itself
        cubic = BindingRules([1 1 1 1], PolyhedronParticleSpecies(Cube(); colors=fill(1, 6)))
        cts = tilings(cubic; maxsize=1)
        # a complete tiling closes all three axes, each through the system's one bond type
        complete = filter(iscomplete, cts)
        @test !isempty(complete)
        @test all(t -> length(translations(t)) == 3 && bondtypes(t) == [1, 1, 1], complete)
        # the partial closures come back alongside: a column with one translation, a sheet with two
        @test sort(unique(length(translations(t)) for t in cts)) == [1, 2, 3]
        @test all(!iscomplete(t) for t in cts if length(translations(t)) < 3)
        @test all(tilingorder(t) == 1 for t in cts)

        # four inert sides leave one axis to close, and one translation closes it
        column = BindingRules([1 1 1 1], PolyhedronParticleSpecies(Cube(); colors=[1, 2, 2, 2, 2, 1]))
        @test length(opensitelocs(first(polygen(column; maxsize=1)))) == 2
        @test all(t -> iscomplete(t) && length(translations(t)) == 1, tilings(column; maxsize=1))

        # distinct colors do not stop a particle tiling. Two faces that face each other start at
        # corresponding corners, so a bond between them turns the neighbour not at all, and a
        # six-color cube tiles the cubic lattice with no twists asked for
        pure(rules) = all(polygen(rules; maxsize=2)) do p
            nparticles(p) < 2 && return true
            g = p.particles[2].pose * inv(p.particles[1].pose)
            return isapprox(rem(rotation_angle(g.psi), 2pi, RoundNearest), 0; atol=1e-8)
        end
        opposites = [1 1 1 6; 1 2 1 5; 1 3 1 4]
        keyedcube = BindingRules(opposites, PolyhedronParticleSpecies(Cube()))
        @test pure(keyedcube)
        kcube = filter(iscomplete, tilings(keyedcube; maxsize=1))
        @test !isempty(kcube)
        # the lattice is the three unit axes, and each bond type is spent once, where the
        # one-color cube spends its single type three times
        @test all(t -> sort(norm.(translations(t))) ≈ [1, 1, 1], kcube)
        @test all(t -> sort(bondtypes(t)) == [1, 2, 3], kcube)

        # a solid whose opposite faces are anti-aligned mates the other way: bonded octahedra meet
        # turned by 180 degrees, so no translation carries one onto the other and there is no cell
        # of one copy
        octopposite(i) = findfirst(j -> isapprox(dot(facenormal(Octahedron(), i),
                                                     facenormal(Octahedron(), j)), -1; atol=1e-8),
                                   1:nfaces(Octahedron()))
        octpairs = unique([minmax(i, octopposite(i)) for i in 1:nfaces(Octahedron())])
        octa = BindingRules(reduce(vcat, [[1 a 1 b] for (a, b) in octpairs]),
                            PolyhedronParticleSpecies(Octahedron()))
        @test all(polygen(octa; maxsize=2)) do p
            nparticles(p) < 2 && return true
            g = p.particles[2].pose * inv(p.particles[1].pose)
            return isapprox(rotation_angle(g.psi), pi; atol=1e-8)
        end
        @test isempty(tilings(octa; maxsize=1))
        # two of them do close along one direction, a column, but never the space
        @test !isempty(tilings(octa; maxsize=2))
        @test !any(iscomplete, tilings(octa; maxsize=2))

        # unkeyed faces reach a lattice without any of this. A non-locking site admits every twist
        # its own symmetry allows, four for a square face, so the translation is among them
        free = PolyhedronParticleSpecies(Cube(); locking=false)
        @test twistfreedom(bindingsite(free, 1), bindingsite(free, 6)) == 4
        freerules = BindingRules(opposites, free)
        @test count(p -> nparticles(p) == 2, polygen(freerules; maxsize=2)) == 12
        @test any(iscomplete, tilings(freerules; maxsize=1))

        # and patches say the same thing as faces, with the twist named directly rather than taken
        # from the corner order
        axes = [SVector(-0.5, 0.0, 0.0), SVector(0.0, -0.5, 0.0), SVector(0.0, 0.0, -0.5),
                SVector(0.0, 0.0, 0.5), SVector(0.0, 0.5, 0.0), SVector(0.5, 0.0, 0.0)]
        keyed = PatchyParticleSpecies(NautyDiGraph(6), 0.5, axes, [0.0, 0.0, 0.0, pi, 0.0, 0.0]; colors=1:6)
        keyedpatches = BindingRules(opposites, keyed)
        @test allunique(color(bindingsite(keyed, i)) for i in 1:6)
        @test pure(keyedpatches)
        kpatch = filter(iscomplete, tilings(keyedpatches; maxsize=1))
        @test !isempty(kpatch)
        @test all(t -> sort(norm.(translations(t))) ≈ [1, 1, 1], kpatch)
        @test all(t -> sort(bondtypes(t)) == [1, 2, 3], kpatch)

        # cells of two cubes come back alongside the one-cube cells
        @test sort(unique(tilingorder(t) for t in tilings(cubic; maxsize=2))) == [1, 2]
    end

    @testset "a tiling describes a structure that exists" begin
        # A tiling claims that copies of its cell at every lattice point form one valid structure.
        # Build that patch and check it against the geometry, which is the half the graph a tiling
        # is identified by cannot speak about.
        cubic = BindingRules([1 1 1 1], PolyhedronParticleSpecies(Cube(); colors=fill(1, 6)))
        for (rules, maxsize) in ((chainlike, 1), (squarerules, 2), (cubic, 1))
            for t in tilings(rules; maxsize)
                cell = unitcell(t)

                # nothing anywhere in the patch overlaps, rests on a face no bond can use, or
                # meets misaligned, all three of which `_overlap_and_contacts` refuses
                center, outer = patch(t; reps=1)
                everything = vcat(center, outer)
                @test all(eachindex(everything)) do i
                    return !first(Roly._overlap_and_contacts(view(everything, 1:(i - 1)), everything[i], rules))
                end

                # every bond the tiling counts is really formed, and is seen once from each side,
                # the cell having a neighbour in both directions along every translation
                _, far = patch(t; reps=2)
                met = sum(p -> length(last(Roly._overlap_and_contacts(far, p, rules))), center)
                @test met == 2 * (nbonds(t) - nbonds(cell))

                # and the sites add up: every bond spends two of them, and what is left over is
                # exactly what the tiling reports as still exposed
                @test 2 * nbonds(t) + length(exposedsitelocs(t)) == nsites(cell)
            end
        end
    end

    @testset "every structure once" begin
        # Two states carrying different graphs for one structure would both survive `unique`, so
        # the check has to come from the geometry. Lay a patch of each tiling down in its canonical
        # frame and record the local environment of every particle of the cell: two cuts of one
        # structure give the same environments, and different structures all but never do.
        function fingerprint(t; radius=2.6)
            center, outer = patch(t; reps=3)
            parts = vcat(center, outer)
            envs = Set{Any}()
            for p in center
                env = sort!([(Roly.speciesindex(q), round.(Tuple(q.pose.psi); digits=5),
                              round.(Tuple(q.pose.x - p.pose.x); digits=5) .+ 0.0)
                             for q in parts if norm(q.pose.x - p.pose.x) <= radius])
                push!(envs, (Roly.speciesindex(p), round.(Tuple(p.pose.psi); digits=5), env))
            end
            return envs
        end

        cases = ((BindingRules([1 1 1 3], UnitSquare), 3, 1),
                 (BindingRules([1 1 1 3; 1 2 1 4], UnitSquare), 3, 22),
                 (BindingRules([1 1 1 2; 1 3 1 4], UnitSquare), 2, 1),
                 (BindingRules([1 1 1 4; 1 2 1 5; 1 3 1 6], UnitHexagon), 2, 11),
                 (BindingRules([1 1 1 1], PolyhedronParticleSpecies(Cube(); colors=fill(1, 6))), 2, 14),
                 (BindingRules([1 1 1 6; 1 2 1 5; 1 3 1 4], UnitCube), 2, 32),
                 (BindingRules([1 4 2 2; 2 1 1 3], UnitSquare), 2, 1))
        for (rules, maxsize, count) in cases
            ts = tilings(rules; maxsize)
            @test length(ts) == count
            @test allunique(ts)
            @test allunique(fingerprint(t) for t in ts)
        end
    end
end
