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
    sys = bindingrules(poly)
    index = Dict(leading_vertex(p) => i for (i, p) in enumerate(poly.particles))
    rows = @NamedTuple{particle::Int, site::Int, color::Int, inert::Bool}[]
    for orig_v in poly.canon2orig
        part = particle_from_leadingvertex(poly, orig_v)
        isnothing(part) && continue
        for k in 1:nsites(part, sys)
            s = bindingsites(part, sys, k)
            _isbound_vertex(poly, part, first(s.vertices)) && continue
            push!(rows, (particle=index[leading_vertex(part)], site=k, color=color(s),
                         inert=isinert(sys, color(s))))
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

Each site's `gauge` and `locking` carry over from the cluster, while its `stab` is recomputed
against the cluster by [`_sitestabilizer`](@ref): in 3D a face is free to twist about its own
normal, and how much of that freedom survives depends on the cluster behind it, not on the face.
"""
MetaParticleSpecies(poly::Polyform; colors=nothing) =
    MetaParticleSpecies(poly, [(r.particle, r.site) for r in exposablesites(poly) if !r.inert];
                        colors)

function MetaParticleSpecies(poly::Polyform{D}, sites; colors=nothing) where {D}
    sys = bindingrules(poly)
    picks = [(Int(p), Int(k)) for (p, k) in sites]
    n = length(picks)
    n > 0 || throw(ArgumentError("a meta-species needs at least one exposed site"))
    allunique(picks) || throw(ArgumentError("`sites` names the same site twice"))

    exposable = Set((r.particle, r.site) for r in exposablesites(poly))
    open = map(picks) do (p, k)
        1 <= p <= nparticles(poly) ||
            throw(ArgumentError("`poly` has $(nparticles(poly)) particles, so ($p, $k) is out of range"))
        part = poly.particles[p]
        1 <= k <= nsites(part, sys) ||
            throw(ArgumentError("particle $p has $(nsites(part, sys)) sites, so ($p, $k) is out of range"))
        (p, k) in exposable ||
            throw(ArgumentError("site ($p, $k) is bound inside the cluster; a bond already " *
                                "consumes it, so it cannot be exposed"))
        return bindingsites(part, sys, k)
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
        # `gauge` and `locking` belong to the site and carry over; `stab` counts the turns about
        # this site's normal that carry the *cluster* onto itself, which only the cluster knows
        return BindingSite(open[i].pose, cols[i], open[i].vertices, open[i].touching_tolerance,
                           open[i].alignment_tolerance, open[i].gauge,
                           _sitestabilizer(poly, open[i]), open[i].locking)
    end
    _labelsites!(g, metasites)

    F = numtype(poly)
    rmax = maximum(poly.particles; init=zero(F)) do part
        norm(part.pose.x) + bounding_radius(species(sys, species_index(part)))
    end
    return MetaParticleSpecies{D,F,eltype(metasites),typeof(g),typeof(poly)}(
        g, metasites, copy(poly), convert(F, rmax))
end

"""
    _sitestabilizer(poly, site)

How many of `site`'s `gauge` turns about its own normal carry `poly` onto itself.

This is the cluster's version of what [`sitestabilizers`](@ref) reports for a rigid species, and
it cannot be read off the site: the turn has to move every particle of the cluster onto a
particle of the same species. `nphases` divides by the attached site's stabilizer and takes the
`lcm` of the two twist freedoms, so overstating this drops attachment phases and understating it
only costs redundant enumeration.

In 2D a site has a single orientation, so the count is one by construction.
"""
function _sitestabilizer(poly::Polyform{D}, site) where {D}
    (D == 2 || site.gauge <= 1) && return 1
    F = numtype(poly)
    sys = bindingrules(poly)
    tosite = inv(site.pose)

    # Every binding site of the cluster, in `site`'s frame, which puts the normal being turned
    # about on the local x axis. The test is over sites rather than particles because a
    # particle's pose is only defined up to its own symmetry -- a cube turned onto itself about a
    # face normal has a different pose but is the same particle -- and a colour identifies a site
    # within the rules, so matching the sites matches the cluster.
    frames = [(color(s), tosite * s.pose, s.gauge)
              for p in poly.particles for s in bindingsites(p, sys)]
    tol = sqrt(eps(F))
    atol = tol * maximum(norm(f.x) for (_, f, _) in frames; init=one(F))

    return count(0:(site.gauge - 1)) do m
        R = RotX(F(2π) * m / site.gauge)
        all(frames) do (c, frame, _)
            any(frames) do (c2, frame2, gauge2)
                # a site's frame is only defined up to its own gauge turns, exactly as
                # `_site_symmetries` matches them
                c2 == c && isapprox(R * frame.x, frame2.x; atol) &&
                    any(psi -> isapprox(R * frame.psi, psi; atol=tol),
                        _siteturns(frame2.psi, gauge2))
            end
        end
    end
end

# Give every open-site vertex a label determined by its color and placed above every interior
# label. Automorphisms then have to preserve what a bond can see (the colors) on top of the
# cluster's interior structure, and equal colors on genuinely equivalent sites stay equivalent.
function _labelsites!(g, sites)
    labs = labels(g)
    sitevertices = Set(v for s in sites for v in s.vertices)
    interior = (labs[v] for v in 1:nv(g) if v ∉ sitevertices)
    base = maximum(interior; init=0)
    palette = sort!(unique(color(s) for s in sites))
    for s in sites
        l = base + searchsortedfirst(palette, color(s))
        for v in s.vertices
            labs[v] = l
        end
    end
    setlabels!(g, labs)
    return g
end

Base.show(io::Core.IO, ps::MetaParticleSpecies) =
    print(io, "$(dimension(ps))d MetaParticleSpecies[", nparticles(ps.cluster), " particles, ",
          nsites(ps), " sites]")

Base.copy(ps::MetaParticleSpecies) =
    typeof(ps)(copy(ps.g), copy(ps.sites), copy(ps.cluster), ps.rmax)

graphrep(ps::MetaParticleSpecies) = ps.g
nsites(ps::MetaParticleSpecies) = length(ps.sites)
bindingsites(ps::MetaParticleSpecies, i::Integer) = ps.sites[i]
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
    sys1, sys2 = bindingrules(s1.cluster), bindingrules(s2.cluster)
    for a in s1.cluster.particles
        pa = species(sys1, species_index(a)) => pose1 * a.pose
        for b in s2.cluster.particles
            pb = species(sys2, species_index(b)) => pose2 * b.pose
            could_contact(pa, pb) || continue
            overlap(pa, pb; kwargs...) && return true
        end
    end
    return false
end
