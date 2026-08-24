"""
    MetaParticleSpecies{D,F,B,G,PF}

A `Polyform` wrapped as a `ParticleSpecies`, so that assemblies can be built out of whole
assemblies.

Chosen unbound sites of the cluster become the species' sites — by default every one that is not
inert, and otherwise exactly those named, see [`exposablesites`](@ref). Geometry is delegated: two
meta-particles overlap exactly when any of their constituent particles do.

Meta-species describe *meta*-assembly systems, whose `BindingRules` say how blocks attach to
blocks. A meta-particle cannot bond to a plain particle; wrap the plain species as its own
single-particle meta-species instead.

The encoding is the wrapped cluster's own graph, not a fresh graph over its open sites. A
cluster's shape is not determined by its binding sites — two open sites can sit in symmetric
poses while the clusters behind them differ — so an encoding built from the sites alone would
have to either claim a symmetry it cannot justify or refuse to claim any. Carrying the cluster
graph instead makes the species' automorphisms *exactly* the cluster's own, so
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
    cluster::PF
    rmax::F
end

"""
    exposablesites(poly::Polyform)

Every *unbound* binding site of `poly`, as rows of `(particle, site, color, inert)`.

These are the sites a [`MetaParticleSpecies`](@ref) can expose, in the order it exposes them by
default. Bound sites are consumed by the bonds holding the cluster together and are never listed:
exposing one would let two blocks bond through a bond that already exists.

`inert` marks a site whose color takes part in no rule of the cluster's own system. Those are
skipped by default and are exactly the sites to name explicitly when a meta-assembly should use
them — selecting a site and giving it a live color is what activates it.
"""
function exposablesites(poly::Polyform)
    rules = bindingrules(poly)
    index = Dict(leadingvertex(p) => i for (i, p) in enumerate(poly.particles))
    rows = @NamedTuple{particle::Int, site::Int, color::Int, inert::Bool}[]
    for orig_v in poly.canon2orig
        part = particle_from_leadingvertex(poly, orig_v)
        isnothing(part) && continue
        for k in 1:nsites(part, rules)
            s = bindingsite(part, rules, k)
            _isbound_vertex(poly, part, first(s.vertices); canonidxs=false) && continue
            push!(rows, (particle=index[leadingvertex(part)], site=k, color=color(s),
                         inert=isinert(rules, color(s))))
        end
    end
    return rows
end

"""
    MetaParticleSpecies(poly::Polyform; colors=nothing)
    MetaParticleSpecies(poly::Polyform, sites; colors=nothing)

Wrap `poly` as a species whose sites are the ones named by `sites`, given as `(particle, site)`
pairs and exposed in that order. Without them every unbound, non-inert site is exposed, in
[`exposablesites`](@ref) order.

  - `colors`: one interaction color per exposed site; by default each keeps the color it has
    inside `poly`

Naming the sites is how a meta-assembly departs from the rules the cluster was built under: a
site that is inert inside `poly` becomes usable simply by being exposed and given a color that
some rule of the new system uses. Colors are the whole interface to those rules, so a fresh
coloring plus a `BindingRules` over it imposes whatever binding behaviour is wanted, and exposing
only one of two equivalent sites correctly costs the meta-species the symmetry that related them.

Each site's `sitesym` and `locking` carry over from the cluster, while its `stab` is recomputed
against the cluster by [`_sitestabilizer`](@ref): in 3D a face is free to twist about its own
normal, and how much of that freedom survives depends on the cluster behind it, not on the face.
"""
MetaParticleSpecies(poly::Polyform; colors=nothing) =
    MetaParticleSpecies(poly, [(r.particle, r.site) for r in exposablesites(poly) if !r.inert];
                        colors)

function MetaParticleSpecies(poly::Polyform{D}, sites; colors=nothing) where {D}
    rules = bindingrules(poly)
    picks = [(Int(p), Int(k)) for (p, k) in sites]
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
        (p, k) in exposable ||
            throw(ArgumentError("site ($p, $k) is bound inside the cluster; a bond already " *
                                "consumes it, so it cannot be exposed"))
        return bindingsite(part, rules, k)
    end

    cols = colors === nothing ? [color(s) for s in open] : collect(Int, colors)
    length(cols) == n ||
        throw(DimensionMismatch("$n sites were exposed but $(length(cols)) colors were given"))

    # The cluster's graph in its *original* vertex order, where every particle owns a contiguous
    # block and each site keeps the vertex range it has inside the cluster. Canonical order would
    # scatter a 3D site's darts, which have to stay contiguous and in their cyclic order for
    # `contact_pairing` to twist a bond correctly.
    g = first(induced_subgraph(graphrep(poly), poly.orig2canon))
    metasites = map(1:n) do i
        # `sitesym` and `locking` belong to the site and carry over; `stab` counts the turns about
        # this site's normal that carry the *cluster* onto itself, which only the cluster knows
        return BindingSite(open[i].pose, cols[i], open[i].vertices, open[i].touching_tolerance,
                           open[i].alignment_tolerance, open[i].sitesym,
                           _sitestabilizer(poly, open[i]), open[i].locking)
    end
    _labelsites!(g, metasites)

    F = numtype(poly)
    rmax = maximum(poly.particles; init=zero(F)) do part
        norm(part.pose.x) + bounding_radius(species(rules, speciesindex(part)))
    end
    return MetaParticleSpecies{D,F,eltype(metasites),typeof(g),typeof(poly)}(
        g, metasites, copy(poly), convert(F, rmax))
end

"""
    _sitestabilizer(poly, site)

How many of `site`'s `sitesym` turns about its own normal carry `poly` onto itself.

This is the cluster's version of what [`stabilizerorders`](@ref) reports for a rigid species, and
it cannot be read off the site: the turn has to move every particle of the cluster onto a
particle of the same species. `_ndistincttwists` divides by the attached site's stabilizer and takes the
`lcm` of the two twist freedoms, so overstating this drops attachment phases and understating it
only costs redundant enumeration.

In 2D a site has a single orientation, so the count is one by construction.
"""
function _sitestabilizer(poly::Polyform{D}, site) where {D}
    (D == 2 || site.sitesym <= 1) && return 1
    F = numtype(poly)
    rules = bindingrules(poly)
    tosite = inv(site.pose)

    # Every binding site of the cluster, in `site`'s frame, which puts the normal being turned
    # about on the local x axis. The test is over sites rather than particles because a
    # particle's pose is only defined up to its own symmetry -- a cube turned onto itself about a
    # face normal has a different pose but is the same particle -- and a colour identifies a site
    # within the rules, so matching the sites matches the cluster.
    frames = [(color(s), tosite * s.pose, s.sitesym)
              for p in poly.particles for s in bindingsites(p, rules)]
    tol = sqrt(eps(F))
    atol = tol * maximum(norm(f.x) for (_, f, _) in frames; init=one(F))

    return count(0:(site.sitesym - 1)) do m
        R = RotX(F(2π) * m / site.sitesym)
        all(frames) do (c, frame, _)
            any(frames) do (c2, frame2, gauge2)
                # a site's frame is only defined up to its own sitesym turns, exactly as
                # `_sitesymmetries` matches them
                c2 == c && isapprox(R * frame.x, frame2.x; atol) &&
                    any(psi -> isapprox(R * frame.psi, psi; atol=tol),
                        _sitetwists(frame2.psi, gauge2))
            end
        end
    end
end

# Give every open-site vertex a label determined by its color and placed above every interior
# label. Automorphisms then have to preserve what a bond can see (the colors) on top of the
# cluster's interior structure, and equal colors on genuinely equivalent sites stay equivalent.
function _labelsites!(g, sites)
    sitevertices = Set(v for s in sites for v in s.vertices)
    base = maximum((label(g, v) for v in 1:nv(g) if v ∉ sitevertices); init=0)
    palette = sort!(unique(color(s) for s in sites))
    for s in sites
        l = base + searchsortedfirst(palette, color(s))
        for v in s.vertices
            setlabel!(g, v, l)
        end
    end
    return g
end

Base.show(io::Core.IO, ps::MetaParticleSpecies) =
    print(io, "$(dimension(ps))d MetaParticleSpecies[", nparticles(ps.cluster), " particles, ",
          nsites(ps), " sites]")

Base.copy(ps::MetaParticleSpecies) =
    typeof(ps)(copy(ps.g), copy(ps.sites), copy(ps.cluster), ps.rmax)

graphrep(ps::MetaParticleSpecies) = ps.g
nsites(ps::MetaParticleSpecies) = length(ps.sites)
bindingsite(ps::MetaParticleSpecies, i::Integer) = ps.sites[i]
bounding_radius(ps::MetaParticleSpecies) = ps.rmax

"""
    cluster(ps::MetaParticleSpecies)

The `Polyform` that `ps` wraps.
"""
cluster(ps::MetaParticleSpecies) = ps.cluster

# The generic `setcolors!` re-derives labels and stabilizers from the site poses, which do not
# describe a cluster. Recoloring here keeps the cluster's interior labeling and rewrites only the
# open-site labels, so the automorphism group tracks the new colors.
function setcolors!(ps::MetaParticleSpecies, colors::AbstractVector{<:Integer})
    length(colors) == nsites(ps) ||
        throw(ArgumentError("expected $(nsites(ps)) colors, got $(length(colors))"))
    for i in eachindex(ps.sites)
        ps.sites[i] = setcolor(ps.sites[i], colors[i])
    end
    _labelsites!(ps.g, ps.sites)
    return nothing
end

"""
    overlap(p1::SpeciesAndPose{<:MetaParticleSpecies}, p2::SpeciesAndPose{<:MetaParticleSpecies})

Whether two posed meta-particles overlap, i.e. whether any of their constituent particles do.
"""
function overlap(p1::SpeciesAndPose{<:MetaParticleSpecies},
                 p2::SpeciesAndPose{<:MetaParticleSpecies}; kwargs...)
    (s1, pose1), (s2, pose2) = p1, p2
    rules1, rules2 = bindingrules(s1.cluster), bindingrules(s2.cluster)
    for a in s1.cluster.particles
        pa = species(rules1, speciesindex(a)) => pose1 * a.pose
        for b in s2.cluster.particles
            pb = species(rules2, speciesindex(b)) => pose2 * b.pose
            could_contact(pa, pb) || continue
            overlap(pa, pb; kwargs...) && return true
        end
    end
    return false
end

"""
    metarules(ps::MetaParticleSpecies)

Lift the interactions of `ps`' cluster to `ps` itself: the [`BindingRules`](@ref) under which two
copies of `ps` bond exactly where the sites they expose would have bonded as ordinary sites.

Only meaningful while `ps` still carries the colors it inherited from the cluster, which is the
default; a recoloring is a statement that the meta-assembly follows rules of its own, and those
have to be written out.
"""
function metarules(ps::MetaParticleSpecies)
    rules = bindingrules(cluster(ps))
    intmat = interactionmatrix(rules)
    cols = [color(bindingsite(ps, i)) for i in 1:nsites(ps)]
    all(c -> 1 <= c <= ncolors(rules), cols) || throw(
        ArgumentError(
            "`ps` carries colors the cluster's rules do not have, so its interactions cannot be " *
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
    rules = bindingrules(cluster(species(metarules, 1)))

    P = eltype(cluster(species(metarules, 1)).particles)
    parts = P[]
    contacts = Tuple{UnitRange{Int},UnitRange{Int},Int,Int}[]
    sites = eltype(typeof(species(metarules, 1).sites))[]
    off = 0
    for mpart in meta.particles
        ps = species(metarules, speciesindex(mpart))
        cl = cluster(ps)
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

Read a meta-polyform as an ordinary [`Polyform`](@ref) over the species its clusters are made of.

Every copy of every cluster becomes a particle in its own right, and every bond becomes a bond:
those the clusters already carried and those the meta-assembly added. The result is the polyform
the particles would have formed had they been placed one at a time.
"""
function unwrap(meta::Polyform{D}) where {D}
    parts, contacts, _ = _unwrapparts(meta)
    rules = bindingrules(cluster(species(bindingrules(meta), 1)))

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

These are the bonds between copies. The ones inside a copy came with its cluster, so
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
