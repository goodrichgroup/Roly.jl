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

Unlike the built-in species, a meta-particle's shape is **not** determined by its binding sites:
two open sites can sit in symmetric poses while the clusters behind them differ. Everything the
built-ins derive from site poses is therefore declared conservatively here instead of computed:

  - graph labels are distinct by default, so no two sites are ever declared equivalent. `σ = 1`,
    and a symmetric cluster is enumerated once per orientation rather than once
  - `stab = 1` on every site: [`nphases`](@ref) divides by the *attached* site's stabilizer, so a
    larger value would skip attachment phases
  - [`check_encoding`](@ref) is deliberately not run, and [`setcolors!`](@ref) does not re-derive
    the labeling, because both read the symmetry off the site poses

Every one of those errs toward redundant enumeration and never toward a missing structure. Pass
`labels` to declare a symmetry you have established by other means.
"""
struct MetaParticleSpecies{D,F,B<:BindingSite,G<:AbstractNautyGraph,PF} <: ParticleSpecies{D,B}
    g::G
    sites::Vector{B}
    cluster::PF
    rmax::F
end

"""
    MetaParticleSpecies(poly::Polyform; colors=nothing, labels=nothing)

Wrap `poly` as a species whose sites are its open binding sites.

  - `colors`: one interaction color per open site; by default each site keeps the color it has
    inside `poly`
  - `labels`: graph labels, one per open site, declaring which sites are equivalent. Distinct by
    default — see the note on [`MetaParticleSpecies`](@ref) before passing anything else

Only 2D clusters are supported. In 3D a site carries a twist freedom greater than one, which the
cycle encoding cannot express and which cannot be declared away: `nphases` takes the `lcm` of the
two sites' twist freedoms, so understating one would drop attachment phases.
"""
function MetaParticleSpecies(poly::Polyform{D}; colors=nothing, labels=nothing) where {D}
    D == 2 || throw(ArgumentError("`MetaParticleSpecies` supports 2D clusters only; a 3D site's " *
                                  "twist freedom needs a dart encoding, and understating it would " *
                                  "drop attachment phases"))
    sys = bindingrules(poly)
    open = collect_open_bindingsites(poly)
    n = length(open)
    n > 0 || throw(ArgumentError("`poly` has no open binding sites to expose"))

    cols = colors === nothing ? [color(s) for s in open] : collect(Int, colors)
    length(cols) == n ||
        throw(DimensionMismatch("`poly` has $n open sites but $(length(cols)) colors were given"))
    labs = labels === nothing ? collect(1:n) : collect(Int, labels)
    length(labs) == n ||
        throw(DimensionMismatch("`poly` has $n open sites but $(length(labs)) labels were given"))

    g, ranges = cycleencoding(n; labels=labs)
    # `gauge` and `locking` are properties of the site itself and carry over; `stab` counts the
    # whole particle's symmetries fixing the site, which for a cluster is not knowable from the
    # site, so it is declared as 1 (see the note on `MetaParticleSpecies`).
    sites = [BindingSite(open[i].pose, cols[i], ranges[i], open[i].touching_tolerance,
                         open[i].alignment_tolerance, open[i].gauge, 1, open[i].locking)
             for i in 1:n]

    F = numtype(poly)
    rmax = maximum(poly.particles; init=zero(F)) do part
        norm(part.pose.x) + bounding_radius(species(sys, species_index(part)))
    end
    return MetaParticleSpecies{D,F,eltype(sites),typeof(g),typeof(poly)}(g, sites, copy(poly),
                                                                        convert(F, rmax))
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

# The generic `setcolors!` re-derives labels and stabilizers from the site poses and then checks
# the encoding against them. For a cluster the poses do not describe the particle, so recoloring
# here changes colors only and leaves the declared labeling alone.
function setcolors!(ps::MetaParticleSpecies, colors::AbstractVector{<:Integer})
    length(colors) == nsites(ps) ||
        throw(ArgumentError("expected $(nsites(ps)) colors, got $(length(colors))"))
    for i in eachindex(ps.sites)
        ps.sites[i] = setcolor(ps.sites[i], colors[i])
    end
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
