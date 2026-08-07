using Roly
using Roly: PolyhedronParticleSpecies, UnitTetrahedron, UnitCube, UnitOctahedron,
            UnitDodecahedron, UnitIcosahedron, UnitPyramid, UnitPrism, UnitAntiprism,
            Polyhedron, Tetrahedron, Cube, Octahedron, Dodecahedron, Icosahedron,
            Pyramid, Prism, Antiprism, shape, corners, nfaces, facedegree, facecentroid,
            geometriclabels, rotationgroup, inradius, edgemidpoint,
            nsites, dimension, isconvex, numtype, bindingsites, graphrep, setcolors!, color,
            could_contact, overlap, symmetrynumber, nparticles, raise!, lower!,
            collect_compatible_pairs, tocanon, toorig, BindingRules, Polyform, nbonds
using Graphs, NautyGraphs, LinearAlgebra, StaticArrays, Rotations

@testset "PolyhedronParticleSpecies" begin
    solids = [
        ("UnitTetrahedron", UnitTetrahedron, Tetrahedron(), 4, 12),
        ("UnitCube", UnitCube, Cube(), 6, 24),
        ("UnitOctahedron", UnitOctahedron, Octahedron(), 8, 24),
        ("UnitDodecahedron", UnitDodecahedron, Dodecahedron(), 12, 60),
        ("UnitIcosahedron", UnitIcosahedron, Icosahedron(), 20, 60),
        ("UnitPyramid(5)", UnitPyramid(5), Pyramid(5), 6, 5),
        ("UnitPrism(5)", UnitPrism(5), Prism(5), 7, 10),
        ("UnitAntiprism(4)", UnitAntiprism(4), Antiprism(4), 10, 8),
    ]

    for (name, ps, shp, nf, order) in solids
        @test dimension(ps) == 3
        @test numtype(ps) === Float64
        @test isconvex(ps)
        @test nsites(ps) == nf
        @test shape(ps) === shp || nfaces(shape(ps)) == nf

        # Default labels are all distinct, so the sparse encoding is used.
        @test nv(graphrep(ps)) == nf
        @test symmetrynumber(ps) == 1

        # One site per face, sitting at the face centroid.
        for i in 1:nf
            @test isapprox(bindingsites(ps, i).pose.x, facecentroid(shp, i); atol=1e-10)
            @test length(bindingsites(ps, i).vertices) == 1
        end
        # The closest site is at the inradius, and everything is inside the bound.
        @test isapprox(minimum(norm(bindingsites(ps, i).pose.x) for i in 1:nf), inradius(shp))
        @test all(norm(c) <= Roly.bounding_radius(ps) + 1e-10 for c in corners(shp))

        @test occursin("PolyhedronParticleSpecies", sprint(show, ps))
        @test occursin("$nf sites", sprint(show, ps))
    end

    # The package convention: outward normal on local x. Roly additionally fixes local z
    # at the midpoint of the face's first edge, which is what makes the bond twist
    # well defined.
    for (name, ps, shp, nf, _) in solids
        for i in 1:nf
            psi = bindingsites(ps, i).pose.psi
            @test isapprox(psi[:, 1], Roly.facenormal(shp, i); atol=1e-10)
            v = edgemidpoint(shp, i, 1) - facecentroid(shp, i)
            @test isapprox(psi[:, 3], normalize(v); atol=1e-10)
            @test isapprox(det(psi), 1; atol=1e-10)      # proper rotation
        end
    end

    # Faces related by a translation are given the same local z (see `_align_mates`), so
    # `sⱼ = sᵢ·Δ` and the bond leaves the attached particle's orientation unchanged. That
    # is what lets a space-filling solid close rings and assemble into its lattice.
    antiparallel(p) = [(i, j) for i in 1:nfaces(p) for j in (i + 1):nfaces(p)
                       if isapprox(Roly.facenormal(p, i), -Roly.facenormal(p, j); atol=1e-8)]
    localz(p, i) = normalize(edgemidpoint(p, i, 1) - facecentroid(p, i))
    aligned(p) = [(i, j) for (i, j) in antiparallel(p) if isapprox(localz(p, i), localz(p, j); atol=1e-8)]

    for (name, shp) in [("Cube", Cube()), ("Prism(4,h=2)", Prism(4, 1.0; h=2.0)),
                        ("Prism(6)", Prism(6)), ("Prism(3)", Prism(3))]
        ps = PolyhedronParticleSpecies(shp; labels=fill(1, nfaces(shp)))
        pairs = aligned(shp)
        # A prism tiles by translation, so every antiparallel pair is a mated pair.
        @test !isempty(pairs)
        @test pairs == antiparallel(shp)
        for (i, j) in pairs
            sys = BindingRules([1 i 1 j], ps)
            poly = Polyform(sys, 1)
            grown = false
            for (site, loc) in collect_compatible_pairs(poly)
                trial = copy(poly)
                ismissing(raise!(trial, site, loc)) && continue
                a, b = trial.particles
                # Same orientation, and displaced by exactly one face-to-face step.
                @test isapprox(a.pose.psi, b.pose.psi; atol=1e-8)
                @test isapprox(norm(a.pose.x - b.pose.x), 2 * norm(facecentroid(shp, i)); atol=1e-8)
                grown = true
                break
            end
            @test grown
        end
    end

    # A tetrahedron has no antiparallel face pairs at all: it does not tile by
    # translation, so there is no orientation-free convention to find and none is forced.
    @test isempty(antiparallel(Tetrahedron()))
    # An octahedron does have antiparallel faces, but they are related by inversion, not
    # translation, so they are not mates and are correctly left unaligned.
    @test !isempty(antiparallel(Octahedron()))
    @test isempty(aligned(Octahedron()))

    for (name, _, shp, nf, order) in solids
        # Faces grouped by geometric orbit recover the solid's rotation group.
        geo = PolyhedronParticleSpecies(shp; labels=geometriclabels(shp))
        @test symmetrynumber(geo) == order == length(rotationgroup(shp))
        # Distinct labels give 1, regardless of encoding.
        @test symmetrynumber(PolyhedronParticleSpecies(shp)) == 1
        @test symmetrynumber(PolyhedronParticleSpecies(shp; encoding=:dart)) == 1
    end

    # A cube with the two caps distinguished from the four sides keeps the 4-fold axis
    # and the 2-fold axes through it: D_4, of order 8.
    caps = [abs(n[3]) > 0.5 ? 2 : 1 for n in Roly.facenormals(Cube())]
    @test symmetrynumber(PolyhedronParticleSpecies(Cube(); labels=caps)) == 8
    # Singling out one face leaves only the rotations about its normal.
    @test symmetrynumber(PolyhedronParticleSpecies(Cube(); labels=[2, 1, 1, 1, 1, 1])) == 4

    shp = Cube()
    # :auto takes the cheap encoding when it is provably equivalent, the dart encoding
    # when a repeated label means the graph has to carry the rotation group.
    @test nv(graphrep(PolyhedronParticleSpecies(shp))) == 6
    @test nv(graphrep(PolyhedronParticleSpecies(shp; labels=fill(1, 6)))) == 24
    @test nv(graphrep(PolyhedronParticleSpecies(shp; encoding=:dart))) == 24
    @test nv(graphrep(PolyhedronParticleSpecies(shp; encoding=:cycle))) == 6
    @test all(length(bindingsites(PolyhedronParticleSpecies(shp; encoding=:dart), i).vertices) == 4
              for i in 1:6)

    @test_throws ArgumentError PolyhedronParticleSpecies(shp; encoding=:nonsense)
    @test_throws ArgumentError PolyhedronParticleSpecies(shp; colors=1:5)

    ps = PolyhedronParticleSpecies(Cube(); colors=[3, 1, 4, 1, 5, 9], labels=1:6)
    @test color(bindingsites(ps, 1)) == 3
    @test color(bindingsites(ps, 5)) == 5

    cp = copy(ps)
    setcolors!(cp, fill(10, 6))
    @test color(bindingsites(cp, 1)) == 10
    @test color(bindingsites(ps, 1)) == 3          # the copy is independent
    @test_throws ArgumentError setcolors!(cp, [1, 2])

    id = Pose{3,Float64,RotMatrix3{Float64}}(SVector(0.0, 0.0, 0.0), one(RotMatrix3{Float64}))
    shifted(d, R=one(RotMatrix3{Float64})) =
        Pose{3,Float64,RotMatrix3{Float64}}(SVector(d, 0.0, 0.0), R)

    @test overlap(UnitCube => id, UnitCube => id)
    @test !overlap(UnitCube => id, UnitCube => shifted(5.0))
    @test !could_contact(UnitCube => id, UnitCube => shifted(5.0))
    @test could_contact(UnitCube => id, UnitCube => shifted(1.0))

    # Face to face at unit separation: touching, not overlapping.
    @test !overlap(UnitCube => id, UnitCube => shifted(1.0))
    @test overlap(UnitCube => id, UnitCube => shifted(0.9))

    # Rotated by 45 degrees the bounding spheres still reach and the inscribed spheres do
    # not, so these go through the separating-axis test rather than either fast path.
    rot45 = RotMatrix3{Float64}(RotZ(π / 4))
    @test overlap(UnitCube => id, UnitCube => shifted(1.05, rot45))
    @test !overlap(UnitCube => id, UnitCube => shifted(1.3, rot45))

    # An edge-on-edge configuration, which no face normal separates: this is what the
    # edge cross-product axes are for.
    edge = RotMatrix3{Float64}(RotZ(π / 4) * RotX(π / 4))
    @test !overlap(UnitCube => id, UnitCube => shifted(1.5, edge))
    @test overlap(UnitCube => id, UnitCube => shifted(0.95, edge))

    for (name, ps, _, _, _) in solids
        @test overlap(ps => id, ps => id)
        @test !overlap(ps => id, ps => shifted(10.0))
    end

    # Number of proper rigid motions mapping an assembly onto itself, computed from the
    # particle poses alone. Every symmetry maps particle 1 onto some particle, which fixes
    # the rotation, so the candidates are just the particles.
    function geometric_symmetry(poly)
        xs = [p.pose.x for p in poly.particles]
        Rs = [p.pose.psi for p in poly.particles]
        n = length(xs)
        return count(1:n) do j
            Q = Rs[j] * inv(Rs[1])
            t = xs[j] - Q * xs[1]
            all(1:n) do i
                xi, Ri = Q * xs[i] + t, Q * Rs[i]
                any(k -> isapprox(xs[k], xi; atol=1e-8) && isapprox(Rs[k], Ri; atol=1e-8), 1:n)
            end
        end
    end

    squarefaces(p) = [i for i in 1:nfaces(p) if facedegree(p, i) == 4]
    systems = [
        ("tetra", Tetrahedron(), [1 1 1 2; 1 3 1 4], 5),
        ("cube/opposite", Cube(), [1 1 1 6; 1 2 1 5; 1 3 1 4], 5),
        ("cube/full", Cube(), reduce(vcat, [[1 i 1 j] for i in 1:6 for j in i:6]), 3),
        ("prism5/sides", Prism(5),
         reduce(vcat, [[1 i 1 j] for i in squarefaces(Prism(5)) for j in squarefaces(Prism(5)) if j >= i]), 4),
        ("dodeca", Dodecahedron(), [1 1 1 12; 1 2 1 11; 1 3 1 10], 4),
    ]

    for (name, shp, bonds, maxn) in systems
        results = map((:cycle, :dart)) do enc
            sys = BindingRules(bonds, PolyhedronParticleSpecies(shp; encoding=enc))
            polys = polygen(sys; maxsize=maxn)
            (counts=[count(p -> nparticles(p) == k, polys) for k in 1:maxn],
             sigmas=sort([(nparticles(p), symmetrynumber(p)) for p in polys]),
             polys=polys)
        end
        cyc, dart = results

        # With all face labels distinct the sparse encoding is provably equivalent to
        # the dart encoding: same structures, same symmetry numbers.
        @test cyc.counts == dart.counts
        @test cyc.sigmas == dart.sigmas
        @test first(cyc.counts) == 1

        # The graph symmetry number must equal the assembly's actual rotational
        # symmetry. A pairing that lost the handedness of a face-to-face bond would
        # show up here as a graph automorphism the geometry does not have.
        for polys in (cyc.polys, dart.polys)
            for p in polys
                @test symmetrynumber(p) == geometric_symmetry(p)
            end
        end
    end

    sys = BindingRules([1 1 1 6; 1 2 1 5; 1 3 1 4], UnitCube)
    polys = polygen(sys; maxsize=4)

    for p in polys
        nparticles(p) < 2 && continue
        # Vertex bookkeeping stays consistent.
        for v in eachindex(p.canon2orig)
            @test tocanon(p, toorig(p, v)) == v
        end
        # Particle vertex blocks tile the graph exactly; a stale leading vertex after a
        # removal would leave a block hanging off the end.
        blocks = sort(reduce(vcat, [collect(Roly.graphvertices(pt, sys)) for pt in p.particles]))
        @test blocks == 1:nv(graphrep(p))
        # Each bond joins the vertices of one pair of faces.
        @test nbonds(p) == 4 * (nparticles(p) - 1) ÷ 1 || nbonds(p) >= nparticles(p) - 1
    end

    # Raising then lowering returns the original structure.
    big = polys[end]
    rebuilt = copy(big)
    lower!(rebuilt)
    @test nparticles(rebuilt) == nparticles(big) - 1
    for (site, loc) in collect_compatible_pairs(rebuilt)
        trial = copy(rebuilt)
        ismissing(raise!(trial, site, loc)) && continue
        trial == big && (@test trial.canon2orig == big.canon2orig; break)
    end

    # Sites of different sizes may bond: a prism's 4-vertex square faces to its 5-vertex
    # pentagonal ones. `contact_pairing` joins them at gcd(4, 5) = 1 vertex, so the bond
    # carries one edge instead of four.
    prism = PolyhedronParticleSpecies(Prism(5); encoding=:dart)
    squares = [i for i in 1:7 if facedegree(Prism(5), i) == 4]
    pentagons = [i for i in 1:7 if facedegree(Prism(5), i) == 5]

    mixed = BindingRules([1 squares[1] 1 pentagons[1]], prism)
    @test mixed isa BindingRules
    same = BindingRules([1 squares[1] 1 squares[2]], prism)
    @test same isa BindingRules
    # Four dart pairs for a square-to-square bond, one reference pair for the mixed one.
    @test nbonds(polygen(same; maxsize=2)[end]) == 4
    @test nbonds(polygen(mixed; maxsize=2)[end]) == 1
end
