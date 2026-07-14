struct MetaParticleSpecies{D,F,B<:BindingSite} <: ParticleSpecies{D,F,B}
    g::NautyDiGraph
    polyform::Polyform
    meta_sites::Vector{B}
    rmax::F
    active_site_indices::Vector{Int}
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
    n = length(site_indices)

    raw_sites = [bindingsites(poly, i) for i in site_indices]
    B = eltype(raw_sites)

    # Symmetry labels from the underlying polyform graph: two meta-sites with the same
    # label are in equivalent positions under the polyform's automorphism group.
    poly_graph_labels = labels(graphrep(poly))
    site_labels = Cint[poly_graph_labels[tocanon(poly, first(s.vertices))] for s in raw_sites]

    if n == 2
        # 2-vertex directed cycle is trivial (Z2 automorphism regardless of labels);
        # use a 4-vertex cycle with each site spanning 2 vertices to avoid this.
        g = NautyDiGraph(cycle_digraph(4); vertex_labels=Cint[site_labels[1], site_labels[1], site_labels[2], site_labels[2]])
        meta_sites = B[BindingSite(raw_sites[1].pose, raw_sites[1].color, 1:2, raw_sites[1].touching_tolerance, raw_sites[1].alignment_tolerance),
                       BindingSite(raw_sites[2].pose, raw_sites[2].color, 3:4, raw_sites[2].touching_tolerance, raw_sites[2].alignment_tolerance)]
    else
        g = if n <= 1
            NautyDiGraph(SimpleDiGraph(n); vertex_labels=site_labels)
        else
            NautyDiGraph(cycle_digraph(n); vertex_labels=site_labels)
        end
        meta_sites = B[BindingSite(s.pose, s.color, k:k, s.touching_tolerance, s.alignment_tolerance)
                       for (k, s) in enumerate(raw_sites)]
    end

    rmax = maximum(
        norm(part.pose.x) + _bounding_radius(species(sys, part.species_index))
        for part in poly.particles; init=zero(F))

    return MetaParticleSpecies{D,F,B}(g, copy(poly), meta_sites, convert(F, rmax),
                                      collect(Int, site_indices))
end

"""
    MetaParticleSpecies(poly::Polyform)

Wrap `poly` as a `ParticleSpecies` using only its open (unbound, non-inert) binding sites.
"""
function MetaParticleSpecies(poly::Polyform)
    sys = assemblysystem(poly)
    active_indices = Int[]
    k = 0
    for orig_v in canonical_vertices(poly)
        part = particle_from_leadingvertex(poly, orig_v)
        isnothing(part) && continue
        for j in 1:nsites(part, sys)
            k += 1
            site = bindingsites(part, sys, j)
            isbound_vertex(poly, tocanon(poly, first(site.vertices))) && continue
            isinert(sys, color(site)) && continue
            push!(active_indices, k)
        end
    end
    return MetaParticleSpecies(poly, active_indices)
end

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
    return nothing
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
    any(!=(orig_sys), assemblysystem(p) for p in polys) &&
        throw(ArgumentError("polyforms do not come from the same assembly system"))
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

function original_assemblysystem(sys::AssemblySystem{D,<:MetaParticleSpecies}) where {D}
    return assemblysystem(species(sys, 1).polyform)
end

function metacolor2origcolor(::AssemblySystem, c::Integer)
    return c
end
function metacolor2origcolor(sys::AssemblySystem{D,<:MetaParticleSpecies}, c::Integer) where {D}
    siteloc = first(color2siteloc(sys, c))
    spcs = species(sys, siteloc[1])
    orig_site_idx = spcs.active_site_indices[siteloc[2]]
    return color(bindingsites(spcs.polyform, orig_site_idx))
end
