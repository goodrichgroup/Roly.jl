# Geometric symmetry of assembled polyforms, checked against the graph's symmetry number.
#
# `sigma` is what identifies structures and what weights them, so if it ever disagreed with the
# geometry the enumeration would merge or split structures with nothing to say so. The species
# constructors check this for *particles* via `_check_encoding`; nothing checked it for the
# assemblies built out of them, which is the case that actually matters.
using Roly
using Roly: bindingsites, nsites, color, _siteturns, raise!, collect_compatible_pairs, Polyform
using LinearAlgebra, StaticArrays, Rotations

"""
Count the proper rotations about a polyform's centroid that map its oriented binding sites onto
themselves, matching position, color, and frame up to each receiving site's own gauge.

Two subtleties, both learned the hard way:

  * a rotation is pinned by where it sends site 1, *unless* two sites share a position -- which
    is exactly what a bond does. So candidates must be deduplicated, or every rotation is
    counted once per site coincident with site 1.
  * the induced site map must be a bijection. Coincident sites let a non-injective map satisfy
    the pointwise conditions.
"""
function site_rotations(poly)
    bs = [bindingsites(poly, i) for i in 1:nsites(poly)]
    isempty(bs) && return []
    c0 = sum(b.pose.x for b in bs) / length(bs)
    xs = [b.pose.x - c0 for b in bs]
    ps = [b.pose.psi for b in bs]
    gs = [b.gauge for b in bs]
    cs = [color(b) for b in bs]
    n = length(bs)
    atol = 1e-7 * max(1.0, maximum(norm, xs))
    turns(i) = collect(_siteturns(ps[i], gs[i]))

    found = typeof(ps[1] * inv(ps[1]))[]
    for a in 1:n
        cs[a] == cs[1] || continue
        for target in turns(a)
            Q = target * inv(ps[1])
            any(U -> isapprox(U, Q; atol=1e-7), found) && continue
            perm = map(1:n) do i
                findfirst(1:n) do j
                    cs[j] == cs[i] && isapprox(Q * xs[i], xs[j]; atol) &&
                        any(u -> isapprox(Q * ps[i], u; atol=1e-7), turns(j))
                end
            end
            any(isnothing, perm) && continue
            length(unique(perm)) == n || continue
            push!(found, Q)
        end
    end
    return found
end

geometric_symmetry(poly) = length(site_rotations(poly))

"""
The strongest reading: rotations mapping the union of the particles' actual corner sets onto
itself. `geometric_symmetry` compares site frames only up to each site's gauge, which is the
resolution the graph works at -- so on its own it could agree with `sigma` by sharing a blind
spot. Checking the solid closes that, for the species that have one.
"""
function body_symmetry(poly, sys, rotations)
    bs = [bindingsites(poly, i) for i in 1:nsites(poly)]
    c0 = sum(b.pose.x for b in bs) / length(bs)
    pts = [pt.pose * c - c0
           for pt in poly.particles
           for c in Roly.corners(Roly.shape(Roly.species(sys, pt.species_index)))]
    return count(Q -> all(p -> any(q -> isapprox(Q * p, q; atol=1e-7), pts), pts), rotations)
end

@testset "polyform symmetry matches the graph" begin
    sidefaces(p) = [i for i in 1:nfaces(p) if abs(Roly.facenormal(p, i)[3]) < 1e-8]
    faced(shp, sticky) = PolyhedronParticleSpecies(
        shp; colors=[i in sticky ? 1 : 2 for i in 1:nfaces(shp)])

    cases = [
        ("triangle",        BindingRules([1 1 1 1], PolygonParticleSpecies(3, 1.0; colors=fill(1, 3))), 4),
        ("square",          BindingRules([1 1 1 1], PolygonParticleSpecies(4, 1.0; colors=fill(1, 4))), 4),
        ("hexagon",         BindingRules([1 1 1 1], PolygonParticleSpecies(6, 1.0; colors=fill(1, 6))), 3),
        ("triangle, keyed", BindingRules([1 1 1 2], UnitTriangle), 4),
        ("patchy disk",     BindingRules([1 1 1 1], PatchyDisk([0.0, 2π/3, 4π/3]; colors=fill(1, 3))), 4),
        ("patchy sphere",   BindingRules([1 1 1 1], PatchySphere(Cube(), 1.0; colors=fill(1, 6))), 3),
        ("polycubes",       BindingRules([1 1 1 1], faced(Cube(), 1:6)), 4),
        ("polyominoes",     BindingRules([1 1 1 1], faced(Prism(4, 1.0; h=2.0), sidefaces(Prism(4, 1.0; h=2.0)))), 4),
        ("polyiamonds",     BindingRules([1 1 1 1], faced(Prism(3, 1.0; h=2.0), sidefaces(Prism(3, 1.0; h=2.0)))), 4),
        ("prism3, h=a",     BindingRules([1 1 1 1], faced(Prism(3), sidefaces(Prism(3)))), 4),
        ("cube, distinct",  BindingRules([1 1 1 2], PolyhedronParticleSpecies(Cube())), 4),
        ("cube, dart enc",  BindingRules([1 1 1 2], PolyhedronParticleSpecies(Cube(); encoding=:dart)), 4),
        ("twisted prism",   BindingRules([1 2 1 2], PolyhedronParticleSpecies(Prism(3);
                                colors=[2, 1, 1, 1, 2], twists=[0, 1, 0, 0, 0])), 4),
        ("free side faces", BindingRules([1 2 1 2], PolyhedronParticleSpecies(Prism(3);
                                colors=[2, 1, 1, 1, 2],
                                locking=[true, false, false, false, true])), 3),
    ]

    for (name, sys, maxsize) in cases
        polys = polygen(sys; maxsize)
        @test !isempty(polys)
        polyhedral = Roly.species(sys, 1) isa PolyhedronParticleSpecies
        for poly in polys
            rots = site_rotations(poly)
            @test symmetrynumber(poly) == length(rots)
            # ...and for a solid, every one of those really does map the body onto itself, so
            # the agreement is not two views sharing the graph's resolution.
            polyhedral && @test body_symmetry(poly, sys, rots) == length(rots)
        end
    end

    # Multi-species too, where the colors of two species must not be conflated.
    disks = [PatchyDisk([0.0, π]; colors=[1, 1]) for _ in 1:2]
    multi = BindingRules(reduce(vcat, [[s 1 t 1] for s in 1:2 for t in 1:2]), disks)
    for poly in polygen(multi; maxsize=4)
        @test symmetrynumber(poly) == geometric_symmetry(poly)
    end
end
