"""
    MetaParticleSpecies{D,F,B,G,PF}

A `Polyform` wrapped as a `ParticleSpecies`, so that meta-assemblies can be built out of whole
polyforms.

Chosen unbound sites of the polyform become the species' sites. By default these are every one that is not
inert, and otherwise exactly those named, see [`exposedsites`](@ref). Geometry is delegated: two
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

Open-site vertices are relabeled by their symmetry orbit, above every interior label, so that
automorphisms have to preserve what a bond can see while interior structure keeps its own
distinctions. The orbits come from the polyform's own [`rotationgroup`](@ref
rotationgroup(::Polyform)), which is also what [`rotationgroup`](@ref
rotationgroup(::ParticleSpecies)) reports for the species: deriving a group from the sites would
claim symmetries the polyform does not have. [`setcolors!`](@ref) maintains the labeling and is
not the generic method, which would re-derive it from the site poses alone.
"""
struct MetaParticleSpecies{D,F,B<:BindingSite,G<:AbstractNautyGraph,PF} <: ParticleSpecies{D,B}
    g::G
    sites::Vector{B}
    poly::PF
    rmax::F
end

"""
    MetaParticleSpecies(poly::Polyform; colors=nothing)
    MetaParticleSpecies(poly::Polyform, sites; colors=nothing)

Wrap `poly` as a particle species whose sites are the ones named by `sites`, given as `(particle, site)`
pairs and exposed in that order. Without them every unbound, non-inert site is exposed, in
[`exposedsites`](@ref) order.

  - `colors`: one interaction color per exposed site; by default each keeps the color it has
    inside `poly`

Each site's `sitesym` and `locking` carry over from the polyform, while its `stab` is recomputed
against the polyform by [`siteorbits`](@ref) and [`stabilizerorders`](@ref).
"""
function MetaParticleSpecies(poly::Polyform; colors=nothing)
    MetaParticleSpecies(poly, opensites(poly); colors)
end

function MetaParticleSpecies(poly::Polyform{D}, sites; colors=nothing) where {D}
    rules = bindingrules(poly)
    picks = [(Int(p), Int(k)) for (p, k) in sites]  # materialize and narrow to Int
    n = length(picks)
    n > 0 || throw(ArgumentError("a meta-species needs at least one exposed site"))
    allunique(picks) || throw(ArgumentError("`sites` names the same site twice"))

    exposable = Set(exposedsites(poly))
    open = map(picks) do (p, k)
        1 <= p <= nparticles(poly) ||
            throw(ArgumentError("`poly` has $(nparticles(poly)) particles, so ($p, $k) is out of range"))
        part = poly.particles[p]
        1 <= k <= nsites(part, rules) ||
            throw(ArgumentError("particle $p has $(nsites(part, rules)) sites, so ($p, $k) is out of range"))
        ParticleSite(p, k) in exposable || throw(
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
    # are the ordinary ones, measured against the polyform's own symmetry group rather than one
    # derived from the sites, which cannot see the polyform behind them.
    group = rotationgroup(poly)
    c = rotationcenter(poly)
    poses, sitesyms = [s.pose + (-c) for s in open], [s.sitesym for s in open]
    orbits = siteorbits(poses, sitesyms, cols; group)
    stabs = stabilizerorders(poses, sitesyms, cols; group)
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
    return check_encoding(
        MetaParticleSpecies{D,F,eltype(metasites),typeof(g),typeof(poly)}(g, metasites, copy(poly), convert(F, rmax))
    )
end

# The two hooks `check_encoding` and its helpers reach through. A meta-species' symmetries are
# among its polyform's, turning about `rotationcenter` rather than about the species' pose origin,
# so the candidates and the frame the site poses are given in are overridden together.
#
# Candidates, not the group: exposing only some of the polyform's open sites, or coloring two of
# them apart, costs the meta-species symmetries the polyform still has. Every consumer filters,
# so `rotationgroup(ps)` is the subgroup that survives.
_rotationcandidates(ps::MetaParticleSpecies) = rotationgroup(polyform(ps))

function _sitegeometry(ps::MetaParticleSpecies)
    c = rotationcenter(polyform(ps))
    return ([s.pose + (-c) for s in ps.sites], [s.sitesym for s in ps.sites],
            [sitelabel(ps, i) for i in 1:nsites(ps)])
end

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
    poses, sitesyms, _ = _sitegeometry(ps)
    group = rotationgroup(ps)
    orbits = siteorbits(poses, sitesyms, colors; group)
    stabs = stabilizerorders(poses, sitesyms, colors; group)
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

function _metaspecies(meta::Polyform)
    spcs = species(bindingrules(meta))
    first(spcs) isa MetaParticleSpecies ||
        throw(ArgumentError("`meta` is not an assembly of meta-particles"))
    return spcs
end

# The color each site of `ps` carries inside the polyform it was taken from. A site keeps the
# vertex range it has there whatever it is recolored to, so any one of its vertices names it.
function _underlyingcolors(ps::MetaParticleSpecies)
    poly = polyform(ps)
    return map(1:nsites(ps)) do i
        v = first(bindingsite(ps, i).vertices)
        return color(bindingsite(poly, _vertex_to_particle_site(poly, v; canonidxs=false)))
    end
end

# Whether the rules `meta` was built under lift the rules its polyforms were built under: every
# bond they offer between copies has to stand for a bond the underlying particles can make.
#
# This asks about the rules, not about `meta`. A system that offers a bond its particles cannot
# make is not a system of those particles, and the polyforms that happen to avoid that bond are
# no more unwrappable for it -- they are structures of a different system that merely look
# familiar. Deciding it per polyform would make unwrapping a property of the sample rather than
# of the system.
function _checkmetalift(meta::Polyform)
    spcs = _metaspecies(meta)
    base = bindingrules(polyform(first(spcs)))
    all(ps -> bindingrules(polyform(ps)) === base, spcs) || throw(
        ArgumentError("these meta-species wrap polyforms built under different rules, so there are none to unwrap into"),
    )

    intmat = interactionmatrix(base)
    cols = [_underlyingcolors(ps) for ps in spcs]
    for (group1, group2) in bindingrules(meta)._bonded_sites, l1 in group1, l2 in group2
        intmat[cols[l1.species][l1.site], cols[l2.species][l2.site]] && continue
        throw(
            ArgumentError(
                "these meta rules bond site $(l1.site) of meta-species $(l1.species) to site " *
                "$(l2.site) of meta-species $(l2.species), but the sites those stand for do not " *
                "bond under the rules their polyforms were built under. A meta-system offering a " *
                "bond its particles cannot make is not a system of those particles, so nothing " *
                "built under it unwraps, this polyform included. `metarules` builds the rules " *
                "that lift.",
            ),
        )
    end
    return meta
end

# Lay a meta-polyform's copies out as plain particles, renumbering each copy's vertices so that
# every site keeps a range of its own. Returns the particles, every contact between them (in
# `_overlap_and_contacts` form), and the sites each copy exposes, all in the new numbering.
function _unwrapparts(meta::Polyform)
    metarules = bindingrules(meta)
    _metaspecies(meta)
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
            ov && _unwrapfailed(parts, sp, rules)
            append!(contacts, cts)
            push!(parts, sp)
        end
        append!(sites, (shift_vertices(mpart.pose * s, off) for s in ps.sites))
        off += nv(graphrep(ps))
    end
    return parts, contacts, sites
end

# `_overlap_and_contacts` refuses for three different reasons and reports all of them the same
# way, so ask it again to find out which. Only ever reached on the way to an error.
function _unwrapfailed(parts, part, rules)
    refuses(; kwargs...) = first(_overlap_and_contacts(parts, part, rules; kwargs...))
    # A meta-assembly's own `overlap` already ruled overlap out, so this one is ours.
    refuses(; allow_noninteracting=true, allow_misaligned=true) &&
        error("Internal error: a meta-assembly unwrapped to overlapping particles. Please file an issue.")
    why = if refuses(; allow_noninteracting=true)
        "meet at a twist the underlying rules do not allow"
    else
        "meet through a pair of sites the underlying rules leave inert"
    end
    return throw(
        ArgumentError(
            "this meta-polyform does not unwrap: two of its copies $why, so the particles could " *
            "never have assembled into this arrangement on their own. Meta rules bond by the " *
            "meta-species' own colors, which need not stand for a bond of the underlying rules; " *
            "`metarules` builds the ones that do.",
        ),
    )
end

"""
    unwrap(meta::Polyform)

Read a meta-polyform as an ordinary [`Polyform`](@ref) over the species its polyforms are made of.

Every copy of every polyform becomes a particle in its own right, and every bond becomes a bond:
those the polyforms already carried and those the meta-assembly added. The result is the polyform
the particles would have formed had they been placed one at a time.

Not every meta-polyform is one of those. Meta rules bond by the meta-species' own colors, and a
[`MetaParticleSpecies`](@ref) is free to recolor the sites it exposes, so a meta bond need not
stand for a bond the underlying rules allow. Throws an `ArgumentError` when any bond the rules
offer does not, whether or not `meta` itself uses that bond, and when two copies are brought into
contact through sites their species does not expose and the underlying rules do not bond.
"""
function unwrap(meta::Polyform{D}) where {D}
    _checkmetalift(meta)
    parts, contacts, _ = _unwrapparts(meta)
    original_rules = bindingrules(polyform(species(bindingrules(meta), 1)))

    g = NautyDiGraph(0)
    for part in parts
        blockdiag!(g, graphrep(species(original_rules, speciesindex(part))))
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
    return Polyform{D,P,typeof(original_rules),typeof(g)}(g, convert(Int, autg.n), cvs, invperm(cvs), parts, original_rules)
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
    unwrappedsites(meta::Polyform)

Every site the copies of a meta-polyform expose, in [`unwrap`](@ref)'s numbering.

Includes the ones the meta-polyform's own bonds consume, since a cell built from it needs both
what it offers and what it has already spent; [`metabonds`](@ref) names the spent ones.
"""
unwrappedsites(meta::Polyform) = last(_unwrapparts(meta))
