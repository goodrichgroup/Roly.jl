struct MetaParticleSpecies{D,F,P<:Pose{D,F},B<:BindingSite} <: ParticleSpecies{D,F,P}
    g::NautyDiGraph
    polyform::Polyform
    meta_sites::Vector{B}
    rmax::F
    active_site_indices::Vector{Int}  # canonical polyform site indices of the meta-sites
end

"""
    MetaParticleSpecies(poly::Polyform, site_indices)

Wrap `poly` as a `ParticleSpecies` for use in a meta-assembly system. Each element of
`site_indices` refers to a canonical binding site index of `poly` (as returned by
`bindingsites(poly, i)`) and becomes one binding site of the meta-species.

The geometry of the species (overlap, could_contact) is handled by checking all
pairwise particle overlaps across two posed copies of `poly`.
"""
function MetaParticleSpecies(poly::Polyform{D}, site_indices::AbstractVector{<:Integer}) where {D}
    sys = assemblysystem(poly)
    F = numtype(sys)
    P = posetype(sys)
    n = length(site_indices)

    raw_sites = [bindingsites(poly, i) for i in site_indices]
    B = eltype(raw_sites)
    meta_sites = B[BindingSite(s.pose, s.color, k:k, s.touching_tolerance, s.alignment_tolerance)
                   for (k, s) in enumerate(raw_sites)]

    g_base = n > 1 ? path_digraph(n) : SimpleDiGraph(max(n, 1))
    g = NautyDiGraph(g_base; vertex_labels=1:max(n, 1))

    rmax = maximum(
        norm(part.pose.x) + _bounding_radius(species(sys, part.species_index))
        for part in poly.particles; init=zero(F))

    return MetaParticleSpecies{D,F,P,B}(g, copy(poly), meta_sites, convert(F, rmax),
                                        collect(Int, site_indices))
end

"""
    MetaParticleSpecies(poly::Polyform)

Wrap `poly` as a `ParticleSpecies` using only its open (unbound, non-inert) binding sites.
"""
function MetaParticleSpecies(poly::Polyform)
    sys = assemblysystem(poly)
    inv_cv = invperm(canonical_vertices(poly))
    active_indices = Int[]
    k = 0
    for v in canonical_vertices(poly)
        part = particle(poly, v)
        isnothing(part) && continue
        for j in 1:nsites(part, sys)
            k += 1
            site = bindingsites(part, sys, j)
            isbound_vertex(poly, inv_cv[first(site.vertices)]) && continue
            isinert(sys, color(site)) && continue
            push!(active_indices, k)
        end
    end
    return MetaParticleSpecies(poly, active_indices)
end

_bounding_radius(ps::PolygonParticleSpecies) = ps.rmax
_bounding_radius(ps::SphereParticleSpecies) = ps.r

Base.show(io::Core.IO, ps::MetaParticleSpecies{D}) where {D} =
    print(io, "$(D)d MetaParticleSpecies with $(nsites(ps)) meta-sites")

Base.copy(ps::MetaParticleSpecies) =
    typeof(ps)(copy(ps.g), copy(ps.polyform), copy(ps.meta_sites), ps.rmax,
               copy(ps.active_site_indices))

dimension(::MetaParticleSpecies{D}) where {D} = D
graphrep(ps::MetaParticleSpecies) = ps.g
nsites(ps::MetaParticleSpecies) = length(ps.meta_sites)
bindingsites(ps::MetaParticleSpecies, i::Integer) = ps.meta_sites[i]

function setcolors!(ps::MetaParticleSpecies, colors::AbstractVector{<:Integer})
    length(colors) != nsites(ps) && throw(ArgumentError("incorrect number of colors"))
    for k in eachindex(ps.meta_sites)
        s = ps.meta_sites[k]
        ps.meta_sites[k] = BindingSite(s.pose, colors[k], s.vertices, s.touching_tolerance, s.alignment_tolerance)
    end
end

function could_contact(p1::SpeciesAndPose{<:MetaParticleSpecies},
                       p2::SpeciesAndPose{<:MetaParticleSpecies}; kwargs...)
    spcs1, pose1 = p1
    spcs2, pose2 = p2
    return norm(pose1.x - pose2.x) < spcs1.rmax + spcs2.rmax
end

function overlap(p1::SpeciesAndPose{<:MetaParticleSpecies},
                 p2::SpeciesAndPose{<:MetaParticleSpecies}; kwargs...)
    spcs1, pose1 = p1
    spcs2, pose2 = p2
    sys = assemblysystem(spcs1.polyform)
    for part1 in spcs1.polyform.particles, part2 in spcs2.polyform.particles
        world_pose1 = pose1 * part1.pose
        world_pose2 = pose2 * part2.pose
        ps1 = species(sys, part1.species_index)
        ps2 = species(sys, part2.species_index)
        overlap(ps1 => world_pose1, ps2 => world_pose2) && return true
    end
    return false
end

