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

    # TODO: n == 2 leads to a trivial cycle, which leads to issues down the road
    g_base = n > 1 ? cycle_digraph(max(n, 3)) : SimpleDiGraph(max(n, 1))
    g = NautyDiGraph(g_base; vertex_labels=n > 2 ? (1:n) : n == 2 ? (1:3) : 1:1)

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

"""
    AssemblySystem(polys::AbstractVector{<:Polyform})

Construct a meta-assembly system where each polyform in `polys` is a building block.
The open binding sites of each polyform become the binding sites of the corresponding
meta-species. Interactions are inherited from the original assembly system: two meta-sites
interact if and only if their original binding site colors interacted.
"""
function AssemblySystem(polys::AbstractVector{<:Polyform})
    isempty(polys) && throw(ArgumentError("polys must not be empty"))
    orig_sys = assemblysystem(first(polys))
    any(!=(orig_sys), assemblysystem(p) for p in polys) && throw(ArgumentError("polyforms do not come from the same assembly system"))
    orig_intmat = interactionmatrix(orig_sys)

    meta_species = [MetaParticleSpecies(poly) for poly in polys]

    # Capture original colors before AssemblySystem overwrites them via setcolors!
    orig_colors = [[color(ms.meta_sites[k]) for k in 1:nsites(ms)] for ms in meta_species]

    bonds = NTuple{4,Int}[]
    for i in eachindex(meta_species), k in 1:nsites(meta_species[i])
        c_ik = orig_colors[i][k]
        for j in i:length(meta_species)
            lstart = j == i ? k : 1
            for l in lstart:nsites(meta_species[j])
                orig_intmat[c_ik, orig_colors[j][l]] && push!(bonds, (i, k, j, l))
            end
        end
    end

    return AssemblySystem(bonds, meta_species)
end