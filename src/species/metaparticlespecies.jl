"""
    MetaParticleSpecies{D,F,B,G,PF}

A `Polyform` wrapped as a `ParticleSpecies`, so that meta-assemblies can be built out of whole
polyforms.

Chosen unbound sites of the polyform become the species' sites. By default these are every one that is not
inert, and otherwise exactly those named, see [`exposablesites`](@ref). Geometry is delegated: two
meta-particles overlap exactly when any of their constituent particles do.

Meta-species describe *meta*-assembly systems, whose `BindingRules` say how meta-particles
attach to each other. A meta-particle cannot bond to a plain particle; wrap the plain species as
its own single-particle meta-species instead.

The encoding is the wrapped polyform's own graph, not a fresh graph over its open sites. A
polyform's shape is not determined by its binding sites — two open sites can sit in symmetric
poses while the polyforms behind them differ — so an encoding built from the sites alone would
have to either claim a symmetry it cannot justify or refuse to claim any. Carrying the polyform
graph instead makes the species' automorphisms *exactly* the polyform's own, so
[`symmetrynumber`](@ref) and the site equivalences are right by construction, and the graph of a
meta-assembly is the graph its particles would have formed directly. Bound sites stay in the
graph as interior vertices; only the open ones are sites, which is why `nv(graphrep(ps))` exceeds
`nsites(ps)`, as it does for a 3D dart encoding.

Open-site vertices are relabeled by color, above every interior label, so that automorphisms have
to preserve the colors a bond can see while interior structure keeps its own distinctions.
[`setcolors!`](@ref) maintains that and is not the generic method, which would re-derive the
labeling from the site poses; [`check_encoding`](@ref) is likewise not run, since it compares the
graph against those same poses.
"""
struct MetaParticleSpecies{D,F,B<:BindingSite,G<:AbstractNautyGraph,PF} <: ParticleSpecies{D,B}
    g::G
    sites::Vector{B}
    poly::PF
    rmax::F
end

"""
    exposablesites(poly::Polyform)

Every unbound binding site of `poly`, as rows of `(particle, site, color, inert)`.
"""
function exposablesites(poly::Polyform)
    rules = bindingrules(poly)
    index = Dict(leadingvertex(p) => i for (i, p) in enumerate(poly.particles))
    rows = @NamedTuple{particle::Int, site::Int, color::Int, inert::Bool}[]
    # iterate through binding sites and push them if unbound
    for orig_v in poly.canon2orig
        part = particle_from_leadingvertex(poly, orig_v)
        isnothing(part) && continue
        for k in 1:nsites(part, rules)
            s = bindingsite(part, rules, k)
            _isbound_vertex(poly, part, first(s.vertices); canonidxs=false) && continue
            push!(rows, (particle=index[leadingvertex(part)], site=k, color=color(s), inert=isinert(rules, color(s))))
        end
    end
    return rows
end

"""
    MetaParticleSpecies(poly::Polyform; colors=nothing)
    MetaParticleSpecies(poly::Polyform, sites; colors=nothing)

Wrap `poly` as a particle species whose sites are the ones named by `sites`, given as `(particle, site)`
pairs and exposed in that order. Without them every unbound, non-inert site is exposed, in
[`exposablesites`](@ref) order.

  - `colors`: one interaction color per exposed site; by default each keeps the color it has
    inside `poly`

Each site's `sitesym` and `locking` carry over from the polyform, while its `stab` is recomputed
against the polyform by [`_metasymmetries`](@ref).
"""
function MetaParticleSpecies(poly::Polyform; colors=nothing)
    MetaParticleSpecies(poly, [(r.particle, r.site) for r in exposablesites(poly) if !r.inert]; colors)
end

function MetaParticleSpecies(poly::Polyform{D}, sites; colors=nothing) where {D}
    rules = bindingrules(poly)
    picks = [(Int(p), Int(k)) for (p, k) in sites]  # materialize and narrow to Int
    n = length(picks)
    n > 0 || throw(ArgumentError("a meta-species needs at least one exposed site"))
    allunique(picks) || throw(ArgumentError("`sites` names the same site twice"))

    exposable = Set((r.particle, r.site) for r in exposablesites(poly))
    open = map(picks) do (p, k)
        1 <= p <= nparticles(poly) ||
            throw(ArgumentError("`poly` has $(nparticles(poly)) particles, so ($p, $k) is out of range"))
        part = poly.particles[p]
        1 <= k <= nsites(part, rules) ||
            throw(ArgumentError("particle $p has $(nsites(part, rules)) sites, so ($p, $k) is out of range"))
        (p, k) in exposable || throw(
            ArgumentError(
                "site ($p, $k) is bound inside the polyform; a bond already consumes it, so it cannot be exposed",
            ),
        )
        return bindingsite(part, rules, k)
    end

    cols = colors === nothing ? [color(s) for s in open] : collect(Int, colors)
    length(cols) == n || throw(DimensionMismatch("$n sites were exposed but $(length(cols)) colors were given"))

    # The polyform's graph in its original vertex order, where every particle owns a contiguous
    # block and each site keeps the vertex range it has inside the polyform.
    g = graphrep(poly)[poly.orig2canon]
    # `sitesym` and `locking` belong to the site and carry over. The orbits and the stabilizers
    # come from the polyform's own symmetry group.
    perms = _metasymmetries(poly, rotationgroup(poly), open, cols)
    orbits, stabs = _metaorbits(perms, n), _metastabs(perms, n)
    metasites = map(1:n) do i
        return BindingSite(
            open[i].pose,
            cols[i],
            open[i].vertices,
            open[i].touching_tolerance,
            open[i].alignment_tolerance,
            open[i].sitesym,
            stabs[i],
            open[i].locking,
        )
    end
    _labelsites!(g, metasites, orbits)

    F = numtype(poly)
    rmax = maximum(poly.particles; init=zero(F)) do part
        norm(part.pose.x) + bounding_radius(species(rules, speciesindex(part)))
    end
    return MetaParticleSpecies{D,F,eltype(metasites),typeof(g),typeof(poly)}(g, metasites, copy(poly), convert(F, rmax))
end

"""
    _metasymmetries(poly, group, sites, colors)

The site permutations induced by the rotations in `group` that carry every one of `sites` onto a
site of `sites` with the same color: the meta-particle's own symmetry group, acting on the sites
it offers.

The counterpart of [`_eachsitesymmetry`](@ref) for a meta-species, and the reason that one cannot
be used directly: it derives the group from the sites alone, but a polyform's shape is not
determined by its binding sites -- two open sites can sit in symmetric poses with different
polyforms behind them -- so it would claim symmetries the polyform does not have. The group has
to come from [`rotationgroup`](@ref rotationgroup(::Polyform))`(poly)` instead.

A rotation carrying an exposed site onto one that was *not* exposed is dropped. Exposing only one
of two equivalent sites is exactly what costs a meta-species the symmetry that related them.
"""
function _metasymmetries(poly::Polyform, group, sites, colors)
    n = length(sites)
    F = numtype(poly)
    centroid = sum(p.pose.x for p in poly.particles) / nparticles(poly)
    xs = [s.pose.x - centroid for s in sites]
    tol = sqrt(eps(F))
    atol = tol * max(one(F), maximum(norm, xs; init=one(F)))

    perms = Vector{Int}[]
    for Q in group
        perm = map(1:n) do i
            findfirst(1:n) do j
                colors[j] == colors[i] &&
                    isapprox(Q * xs[i], xs[j]; atol) &&
                    any(psi -> isapprox(Q * sites[i].pose.psi, psi; atol=tol),
                        _sitetwists(sites[j].pose.psi, sites[j].sitesym))
            end
        end
        any(isnothing, perm) || push!(perms, collect(Int, perm))
    end
    return perms
end

# Orbits of the exposed sites, as `siteorbits` reports them for a rigid species.
function _metaorbits(perms, n::Integer)
    orbit = collect(1:n)
    for perm in perms, i in 1:n
        lo, hi = minmax(orbit[i], orbit[perm[i]])
        hi == lo || replace!(orbit, hi => lo)
    end
    return _canonicalpartition(orbit)
end

# Stabilizer orders of the exposed sites, as `stabilizerorders` reports them for a rigid species.
# A rotation fixing a site fixes the centroid too, so it turns about the line through both, which
# is the site's normal; each such rotation realizes a distinct one of the site's `sitesym` turns.
_metastabs(perms, n::Integer) = [count(perm -> perm[i] == i, perms) for i in 1:n]

# Label every open-site vertex by its symmetry orbit, placed above every interior label.
#
# This is what `_check_labeling` demands of every other species -- a labeling has to be exactly
# the orbits -- with the orbits taken from the polyform's own rotation group rather than from the
# sites. Labeling by color alone would merge two orbits that happen to share one.
function _labelsites!(g, sites, orbits)
    sitevertices = Set(v for s in sites for v in s.vertices)
    base = maximum((label(g, v) for v in 1:nv(g) if v ∉ sitevertices); init=0)
    for (i, s) in enumerate(sites)
        for v in s.vertices
            setlabel!(g, v, base + orbits[i])
        end
    end
    return g
end

function Base.show(io::Core.IO, ps::MetaParticleSpecies)
    print(io, "$(dimension(ps))d MetaParticleSpecies[", nparticles(ps.poly), " particles, ", nsites(ps), " sites]")
end

Base.copy(ps::MetaParticleSpecies) = typeof(ps)(copy(ps.g), copy(ps.sites), copy(ps.poly), ps.rmax)

graphrep(ps::MetaParticleSpecies) = ps.g
nsites(ps::MetaParticleSpecies) = length(ps.sites)
bindingsite(ps::MetaParticleSpecies, i::Integer) = ps.sites[i]
bounding_radius(ps::MetaParticleSpecies) = ps.rmax

"""
    polyform(ps::MetaParticleSpecies)

The `Polyform` that `ps` wraps.
"""
polyform(ps::MetaParticleSpecies) = ps.poly

# The generic `setcolors!` re-derives labels and stabilizers from the site poses, which do not
# describe a polyform. Recoloring here keeps the polyform's interior labeling and rewrites only the
# open-site labels and stabilizers, both against the polyform's own symmetry group.
function setcolors!(ps::MetaParticleSpecies, colors::AbstractVector{<:Integer})
    n = nsites(ps)
    length(colors) == n || throw(ArgumentError("expected $n colors, got $(length(colors))"))
    perms = _metasymmetries(ps.poly, rotationgroup(ps.poly), ps.sites, colors)
    orbits, stabs = _metaorbits(perms, n), _metastabs(perms, n)
    for i in eachindex(ps.sites)
        ps.sites[i] = setstab(setcolor(ps.sites[i], colors[i]), stabs[i])
    end
    _labelsites!(ps.g, ps.sites, orbits)
    return nothing
end

"""
    overlap(p1::SpeciesAndPose{<:MetaParticleSpecies}, p2::SpeciesAndPose{<:MetaParticleSpecies})

Whether two posed meta-particles overlap, i.e. whether any of their constituent particles do.
"""
function overlap(p1::SpeciesAndPose{<:MetaParticleSpecies}, p2::SpeciesAndPose{<:MetaParticleSpecies}; kwargs...)
    (s1, pose1), (s2, pose2) = p1, p2
    rules1, rules2 = bindingrules(s1.poly), bindingrules(s2.poly)
    for a in s1.poly.particles
        pa = species(rules1, speciesindex(a)) => pose1 * a.pose
        for b in s2.poly.particles
            pb = species(rules2, speciesindex(b)) => pose2 * b.pose
            could_contact(pa, pb) || continue
            overlap(pa, pb; kwargs...) && return true
        end
    end
    return false
end

"""
    metarules(ps::MetaParticleSpecies)

Lift the interactions of `ps`' polyform to `ps` itself: the [`BindingRules`](@ref) under which two
copies of `ps` bond exactly where the sites they expose would have bonded as ordinary sites.

Only meaningful while `ps` still carries the colors it inherited from the polyform, which is the
default; a recoloring is a statement that the meta-assembly follows rules of its own, and those
have to be written out.
"""
function metarules(ps::MetaParticleSpecies)
    rules = bindingrules(polyform(ps))
    intmat = interactionmatrix(rules)
    cols = [color(bindingsite(ps, i)) for i in 1:nsites(ps)]
    all(c -> 1 <= c <= ncolors(rules), cols) || throw(
        ArgumentError(
            "`ps` carries colors the polyform's rules do not have, so its interactions cannot be " *
            "lifted from them. Write the meta rules out instead.",
        ),
    )
    rows = [[1, i, 1, j] for i in eachindex(cols) for j in i:length(cols) if intmat[cols[i], cols[j]]]
    bonds = isempty(rows) ? zeros(Int, 0, 4) : permutedims(reduce(hcat, rows))
    return BindingRules(bonds, [ps])
end

# Lay a meta-polyform's copies out as plain particles, renumbering each copy's vertices so that
# every site keeps a range of its own. Returns the particles, every contact between them (in
# `_overlap_and_contacts` form), and the sites each copy exposes, all in the new numbering.
function _unwrapparts(meta::Polyform)
    metarules = bindingrules(meta)
    first(species(metarules)) isa MetaParticleSpecies ||
        throw(ArgumentError("`meta` is not an assembly of meta-particles"))
    rules = bindingrules(polyform(species(metarules, 1)))

    P = eltype(polyform(species(metarules, 1)).particles)
    parts = P[]
    contacts = Tuple{UnitRange{Int},UnitRange{Int},Int,Int}[]
    sites = eltype(typeof(species(metarules, 1).sites))[]
    off = 0
    for mpart in meta.particles
        ps = species(metarules, speciesindex(mpart))
        cl = polyform(ps)
        for p in cl.particles
            sp = P(mpart.pose * p.pose, p.leadingvertex + off, p.speciesindex)
            ov, cts = _overlap_and_contacts(parts, sp, rules)
            # A valid meta-assembly never overlaps: its own `overlap` already ruled that out.
            ov && error("Internal error: a meta-assembly unwrapped to overlapping particles. Please file an issue.")
            append!(contacts, cts)
            push!(parts, sp)
        end
        append!(sites, (shift_vertices(mpart.pose * s, off) for s in ps.sites))
        off += nv(graphrep(ps))
    end
    return parts, contacts, sites
end

"""
    unwrap(meta::Polyform)

Read a meta-polyform as an ordinary [`Polyform`](@ref) over the species its polyforms are made of.

Every copy of every polyform becomes a particle in its own right, and every bond becomes a bond:
those the polyforms already carried and those the meta-assembly added. The result is the polyform
the particles would have formed had they been placed one at a time.
"""
function unwrap(meta::Polyform{D}) where {D}
    parts, contacts, _ = _unwrapparts(meta)
    rules = bindingrules(polyform(species(bindingrules(meta), 1)))

    g = NautyDiGraph(0)
    for part in parts
        blockdiag!(g, graphrep(species(rules, speciesindex(part))))
    end
    for (vs1, vs2, t, ntwists) in contacts
        for (v1, v2) in contact_pairing(vs1, vs2, t, ntwists)
            add_edge!(g, v1, v2)
            add_edge!(g, v2, v1)
        end
    end

    # `g` is built in original vertex order, so the canonical permutation is `canon2orig` itself.
    perm, autg = nauty(g; canonize=true)
    cvs = collect(Int, perm)
    P = eltype(parts)
    return Polyform{D,P,typeof(rules),typeof(g)}(g, convert(Int, autg.n), cvs, invperm(cvs), parts, rules)
end

"""
    metabonds(meta::Polyform)

The bonds a meta-polyform adds, as pairs of binding site vertex ranges in [`unwrap`](@ref)'s
numbering.

These are the bonds between copies. The ones inside a copy came with its polyform, so
`bonds(unwrap(meta))` is the two sets together.
"""
function metabonds(meta::Polyform)
    _, contacts, sites = _unwrapparts(meta)
    exposed = Set(s.vertices for s in sites)
    return [(c[1], c[2]) for c in contacts if c[1] in exposed && c[2] in exposed]
end

"""
    exposedsites(meta::Polyform)

Every site the copies of a meta-polyform expose, in [`unwrap`](@ref)'s numbering.

Includes the ones the meta-polyform's own bonds consume, since a cell built from it needs both
what it offers and what it has already spent; [`metabonds`](@ref) names the spent ones.
"""
exposedsites(meta::Polyform) = last(_unwrapparts(meta))
