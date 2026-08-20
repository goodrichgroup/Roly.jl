"""
    MetaParticleSpecies{D,F,B,G,PF}

A `Polyform` wrapped as a `ParticleSpecies`, so that assemblies can be built out of whole
assemblies.

The cluster's open binding sites become the species' sites, in the order
[`collect_open_bindingsites`](@ref) returns them. Geometry is delegated: two meta-particles
overlap exactly when any of their constituent particles do.

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
    MetaParticleSpecies(poly::Polyform; colors=nothing)

Wrap `poly` as a species whose sites are its open binding sites.

  - `colors`: one interaction color per open site; by default each site keeps the color it has
    inside `poly`

Only 2D clusters are supported. In 3D a site carries a twist freedom greater than one, which the
site's stabilizer within the cluster has to resolve: `nphases` takes the `lcm` of the two sites'
twist freedoms and divides by the attached site's stabilizer, so getting it wrong drops
attachment phases rather than merely costing work.
"""
function MetaParticleSpecies(poly::Polyform{D}; colors=nothing) where {D}
    D == 2 || throw(ArgumentError("`MetaParticleSpecies` supports 2D clusters only; in 3D a " *
                                  "site's twist freedom needs its stabilizer within the cluster, " *
                                  "and getting that wrong drops attachment phases"))
    open = collect_open_bindingsites(poly)
    n = length(open)
    n > 0 || throw(ArgumentError("`poly` has no open binding sites to expose"))

    cols = colors === nothing ? [color(s) for s in open] : collect(Int, colors)
    length(cols) == n ||
        throw(DimensionMismatch("`poly` has $n open sites but $(length(cols)) colors were given"))

    # the cluster's graph, with each open site's vertices carried over to canonical numbering
    g = copy(graphrep(poly))
    sites = map(1:n) do i
        vs = [tocanon(poly, v) for v in open[i].vertices]
        # a site occupies one vertex in 2D, so its canonical vertices stay contiguous
        r = minimum(vs):maximum(vs)
        length(r) == length(vs) ||
            throw(ArgumentError("open site $i does not occupy a contiguous vertex range"))
        # `gauge` and `locking` belong to the site; `stab` counts the symmetries of the whole
        # particle fixing it, which in 2D is 1 whenever the twist freedom is
        return BindingSite(open[i].pose, cols[i], r, open[i].touching_tolerance,
                           open[i].alignment_tolerance, open[i].gauge, 1, open[i].locking)
    end
    _labelsites!(g, sites)

    F = numtype(poly)
    sys = bindingrules(poly)
    rmax = maximum(poly.particles; init=zero(F)) do part
        norm(part.pose.x) + bounding_radius(species(sys, species_index(part)))
    end
    return MetaParticleSpecies{D,F,eltype(sites),typeof(g),typeof(poly)}(g, sites, copy(poly),
                                                                        convert(F, rmax))
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
