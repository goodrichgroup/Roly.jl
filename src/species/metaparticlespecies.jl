"""
    MetaParticleSpecies{D,F,B,G,PF}

A `Polyform` wrapped as a `ParticleSpecies`, so that meta-assemblies can be built out of whole
polyforms.

Selected exposed sites of the polyform become the species' sites (by default are all non-inert sites are forwarded).
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
pairs and exposed in that order. Without them every open non-inert site is exposed, in
[`exposedsites`](@ref) order.

  - `colors`: one interaction color per exposed site; by default each keeps the color it has
    inside `poly`
  - `exposeinert`: whether to also expose previously inert sites to the new binding rules

Each site's `sitesym` and `locking` carry over from the polyform, while its `stab` is recomputed
against the polyform by [`siteorbits`](@ref) and [`stabilizerorders`](@ref).
"""
function MetaParticleSpecies(poly::Polyform; colors=nothing, exposeinert::Bool=false)
    return MetaParticleSpecies(poly, exposeinert ? exposedsites(poly) : opensites(poly); colors)
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
                "site ($p, $k) is bound inside the polyform; a bond already consumes it, so it cannot be exposed"
            ),
        )
        return bindingsite(part, rules, k)
    end

    cols = colors === nothing ? [color(s) for s in open] : collect(Int, colors)
    length(cols) == n || throw(DimensionMismatch("$n sites were exposed but $(length(cols)) colors were given"))

    # The polyform's graph in its original vertex order, where every particle owns a contiguous
    # block of vertices and each site keeps the vertex range it has inside the polyform.
    g = graphrep(poly)[poly.orig2canon]

    # `sitesym` and `locking` belong to the site and carry over. The orbits and the stabilizers
    # are the ordinary ones, measured against the polyform's own symmetry group rather than one
    # derived from just the site poses, which would be unaware of shape of the polyform.
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

# candidate rotation symmetries, which is the symmetry of the underlying polyform.
# Exposing only some of the polyform's open sites or recoloring some sites will lower the symmetry
_rotationcandidates(ps::MetaParticleSpecies) = rotationgroup(polyform(ps))

function _sitegeometry(ps::MetaParticleSpecies)
    c = rotationcenter(polyform(ps))
    return ([s.pose + (-c) for s in ps.sites], [s.sitesym for s in ps.sites], [sitelabel(ps, i) for i in 1:nsites(ps)])
end

# Label every open-site vertex by its symmetry orbit, placed above every interior label.
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
    BindingRules(ps::MetaParticleSpecies)

Lift the interactions of the wrapped polyform to the meta-species itself: the rules under which
two meta-particles bond exactly where the sites they expose would have bonded as ordinary sites.

Only meaningful while the meta-species still carries the colors it inherited from its polyform,
which is the default.
"""
BindingRules(ps::MetaParticleSpecies) = BindingRules([ps])

"""
    BindingRules(pss::AbstractVector{<:MetaParticleSpecies})

Lift the interactions of the wrapped polyforms to the meta-species themselves: the rules under
which two meta-particles bond exactly where the sites they expose would have bonded as ordinary
sites, across the species as well as within each.

Only works if the meta-species still carry the colors they inherited from their
polyforms, which is the default. All of them must wrap polyforms built under the same rules.
"""
function BindingRules(pss::AbstractVector{<:MetaParticleSpecies})
    isempty(pss) && throw(ArgumentError("The list of particle species is empty."))
    base = bindingrules(polyform(first(pss)))
    all(ps -> bindingrules(polyform(ps)) === base, pss) ||
        throw(ArgumentError("The meta particle species wrap polyforms built under different binding rules."))

    intmat = interactionmatrix(base)
    cols = [[color(bindingsite(ps, i)) for i in 1:nsites(ps)] for ps in pss]
    all(c -> 1 <= c <= ncolors(base), Iterators.flatten(cols)) ||
        throw(ArgumentError("These meta particle species carry colors the original polyforms' rules do not have."))

    # Each unordered pair of sites once: within a species from `i` on, across species all of them.
    rows = [
        [a, i, b, j] for a in eachindex(pss) for i in eachindex(cols[a]) for b in a:length(pss) for
        j in (a == b ? i : 1):length(cols[b]) if intmat[cols[a][i], cols[b][j]]
    ]
    bonds = isempty(rows) ? zeros(Int, 0, 4) : permutedims(reduce(hcat, rows))
    return BindingRules(bonds, pss)
end

"""
    MetaBindingRules{D}
    MetaPolyform{D,P}

Rules whose species are all [`MetaParticleSpecies`](@ref), and a `Polyform` built under such
rules -- the two things the meta level is made of, named so that dispatch can say so. Rules
mixing meta-species with plain ones are neither, since their element type is then the abstract
supertype, which is what makes a meta-particle bonded to a plain one a `MethodError`.
"""
const MetaBindingRules{D} = BindingRules{D,<:MetaParticleSpecies}
const MetaPolyform{D,P} = Polyform{D,P,<:MetaBindingRules{D}}

# The color each site of `ps` carries inside the polyform it was taken from. A site keeps the
# vertex range it has there whatever it is recolored to, so any one of its vertices names it.
function _underlyingcolors(ps::MetaParticleSpecies)
    poly = polyform(ps)
    return map(1:nsites(ps)) do i
        # look up the original polyform binding site through the graph vertex of the species' site
        v = first(bindingsite(ps, i).vertices)
        return color(bindingsite(poly, _vertex_to_particle_site(poly, v; canonidxs=false)))
    end
end

"""
    inducedrules(rules::BindingRules)

"Project" the meta-rules `rules` to down to the original binding rules of the underlying particles.
Returns the binding rules of the original particles with additional bonds that make all contacts in the meta-polyform
valid.

Meta rules bond by colors of the meta-species' own choosing, which need not stand for a bond the
underlying particles can make. Rather than refuse those, this adds these to the original binding rules as additional
bond types between the underlying particle species.
"""
function inducedrules(rules::MetaBindingRules)
    spcs = species(rules)
    origrules = bindingrules(polyform(first(spcs)))
    all(ps -> bindingrules(polyform(ps)) === origrules, spcs) ||
        throw(ArgumentError("The meta-species wrap polyforms built under different binding rules."))

    intmat = Matrix(interactionmatrix(origrules))
    cols = map(_underlyingcolors, spcs)
    # add bonds to the original binding rules
    for (sites1, sites2) in bonded_sites(rules)
        for l1 in sites1, l2 in sites2
            c1, c2 = cols[l1.species][l1.site], cols[l2.species][l2.site]
            intmat[c1, c2] = intmat[c2, c1] = true
        end
    end
    return BindingRules(intmat, species(origrules))
end

# Species `i` of `from` and species `i` of `rules` should be a "compatible" species. What that means exactly can be
# refined later. Right now, we check the vertex count and site poses.
# TODO: this is overly restrictive
function _checkspecies(from::BindingRules, rules::BindingRules, i::Integer)
    ok = i <= nspecies(rules)
    if ok
        a, b = species(from, i), species(rules, i)
        ok =
            nsites(a) == nsites(b) &&
            nv(graphrep(a)) == nv(graphrep(b)) &&
            all(k -> _siteoverlap(bindingsite(a, k), bindingsite(b, k)), 1:nsites(a))
    end
    ok || throw(
        ArgumentError(
            "species $i of `rules` is not the species that would stand in its place. `rules` has " *
            "to list the corresponding species in the same order.",
        ),
    )
    return nothing
end

# Check every species of `poly`, against the ones `rules` has in those places.
function _checkspecies(poly::Polyform, rules::BindingRules)
    foreach(q -> _checkspecies(bindingrules(poly), rules, speciesindex(q)), poly.particles)
    return poly
end

# check if the sites are at the same location and pose
function _siteoverlap(sa::BindingSite, sb::BindingSite)
    return isapprox(sa.pose.x, sb.pose.x; atol=sa.touching_tolerance + sb.touching_tolerance, rtol=0) &&
           isapprox(sa.pose.psi, sb.pose.psi; atol=sa.alignment_tolerance + sb.alignment_tolerance, rtol=0)
end

# Return the polyforms each species of `src` is replaced by, one per species, as a polyform of `rules`.
# Return `given` if present, return the underlying polyform of a meta species, or simply return the corresponding species.
function _resolvesubstitutions(src::BindingRules, rules::BindingRules, given)
    return map(collect(enumerate(species(src)))) do (i, ps)
        sub = get(given, i, nothing)
        # written by hand: a polyform of `rules` already, so its particles wear its own species
        if !isnothing(sub)
            bindingrules(sub) === rules || throw(ArgumentError("A substitution has to be a polyform of `rules`."))
            return sub
        end
        # a meta-species: a polyform of the rules it was lifted from, which `inducedrules` extends
        # rather than reuses, so its particles' species have to line up with `rules`
        ps isa MetaParticleSpecies && return _checkspecies(polyform(ps), rules)
        # standing for itself: `rules` has to have that same species in that same place
        _checkspecies(src, rules, i)
        return Polyform(rules, i)
    end
end

# Expand the meta-species; return all underyling particles of the underlying binding rules
function _underlying_particles(meta::MetaPolyform)
    subs = [polyform(ps) for ps in species(bindingrules(meta))]
    return _substitute_particles(meta, subs)
end

# Substitute every particle of `poly` of species `i` with the particles of polyform `subs[i]`.
function _substitute_particles(poly::Polyform, subs)
    P = particletype(first(subs))
    parts = P[]
    off = 0
    for part in poly.particles
        sub = subs[speciesindex(part)]
        subrules = bindingrules(sub)
        for q in sub.particles
            si = speciesindex(q)
            push!(parts, P(part.pose * q.pose, off + 1, si))
            off += nv(graphrep(species(subrules, si)))
        end
    end
    return parts
end

"""
    recast(poly::Polyform, rules::BindingRules; substitutions=Dict())

Recast `poly` as a [`Polyform`](@ref) of `rules`, substituting every particle with one or multiple particles
from `rules`.

`substitutions` maps a species index of `poly`'s own rules to the `Polyform` that species is
replaced by. It has to be a polyform of `rules`. A [`MetaParticleSpecies`](@ref) is substituted by the polyform it
wraps and needs no entry. Any other species type not in `substitutions`is relplaced by the species of `rules`
with the same species index.
"""
function recast(poly::Polyform{D}, rules::BindingRules; substitutions=Dict()) where {D}
    src = bindingrules(poly)
    subs = _resolvesubstitutions(src, rules, substitutions)

    parts = _substitute_particles(poly, subs)
    contacts = Contact[]
    g = NautyDiGraph(0)
    for (i, sp) in enumerate(parts)
        placed = view(parts, 1:(i - 1))
        overlap, cts = _overlap_and_contacts(placed, sp, rules)
        overlap && _recastfailed(placed, sp, rules)
        append!(contacts, cts)
        blockdiag!(g, graphrep(species(rules, speciesindex(sp))))
    end

    for contact in contacts
        for (v1, v2) in contact_pairing(contact)
            add_edge!(g, v1, v2)
            add_edge!(g, v2, v1)
        end
    end

    # `g` is built in original vertex order, so the canonical permutation is `canon2orig` itself.
    perm, autg = nauty(g; canonize=true)
    cvs = collect(Int, perm)
    return Polyform{D,particletype(rules),typeof(rules),typeof(g)}(
        g, convert(Int, autg.n), cvs, invperm(cvs), parts, rules
    )
end

# `_overlap_and_contacts` refuses for three different reasons and reports all of them the same
# way, so ask it again to find out which. Only ever reached on the way to an error.
function _recastfailed(parts, part, rules)
    refuses(; kwargs...) = first(_overlap_and_contacts(parts, part, rules; kwargs...))

    refuses(; allow_noninteracting=true, allow_misaligned=true) &&
        throw(ArgumentError("two of the particles overlap, the recast polyform is invalid."))
    why = if refuses(; allow_noninteracting=true)
        "at a twist `rules` does not allow"
    else
        "at a pair of sites `rules` leaves inert"
    end
    return throw(
        ArgumentError("Two of the particles touch $why, so recasting does not result in a polyform valid under `rules`")
    )
end

