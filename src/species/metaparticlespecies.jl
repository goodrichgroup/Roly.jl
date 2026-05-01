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

_bounding_radius(ps::MetaParticleSpecies) = ps.rmax
_bounding_radius(ps::PolygonParticleSpecies) = ps.rmax
_bounding_radius(ps::SphereParticleSpecies) = ps.r

function _original_color(sys::AssemblySystem{D,<:MetaParticleSpecies}, c::Integer) where {D}
    siteloc = color2siteloc(sys, c)[1]
    ms = species(sys, siteloc[1])
    return color(bindingsites(ms.polyform, ms.active_site_indices[siteloc[2]]))
end
_original_bonded_colors(sys::AssemblySystem{D,<:MetaParticleSpecies}) where {D} =
    bonded_colors(assemblysystem(species(sys, 1).polyform))

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

    active1 = spcs1.active_site_indices
    active2 = spcs2.active_site_indices
    n1 = nsites(spcs1.polyform)
    n2 = nsites(spcs2.polyform)
    for i in 1:n1
        is_active1 = i ∈ active1
        site1 = pose1 * bindingsites(spcs1.polyform, i)
        for j in 1:n2
            is_active1 && j ∈ active2 && continue  # valid meta-bond candidate, not a clash
            site2 = pose2 * bindingsites(spcs2.polyform, j)
            istouching(site1, site2) && return true
        end
    end
    return false
end

"""
    base_composition(p::Polyform)

Return the composition of a meta-polyform `p` as if it were assembled from the original
base-level particles. Sums the compositions of each constituent meta-particle's wrapped
polyform and adds one count per inter-meta-particle bond translated back to its original
bond type.
"""
function base_composition(p::Polyform)
    meta_sys = assemblysystem(p)
    isempty(p.particles) && throw(ArgumentError("polyform has no particles"))
    meta_ps1 = species(meta_sys, first(p.particles).species_index)
    meta_ps1 isa MetaParticleSpecies || throw(ArgumentError("polyform must be made of MetaParticleSpecies"))
    orig_sys = assemblysystem(meta_ps1.polyform)

    ns = nspecies(orig_sys)
    nb = nbonds(orig_sys)
    comp = zeros(Int, ns + nb)

    for part in p.particles
        comp .+= composition(species(meta_sys, part.species_index).polyform)
    end

    g = graphrep(p)
    labs = labels(g)
    orig_bondlist = bonded_colors(orig_sys)

    for (; src, dst) in edges(g)
        src >= dst && continue
        has_edge(g, dst, src) || continue
        src_siteloc = color2siteloc(meta_sys, label2color(meta_sys, labs[src]))[1]
        dst_siteloc = color2siteloc(meta_sys, label2color(meta_sys, labs[dst]))[1]
        src_meta_ps = species(meta_sys, src_siteloc[1])
        dst_meta_ps = species(meta_sys, dst_siteloc[1])
        src_orig_color = color(bindingsites(src_meta_ps.polyform, src_meta_ps.active_site_indices[src_siteloc[2]]))
        dst_orig_color = color(bindingsites(dst_meta_ps.polyform, dst_meta_ps.active_site_indices[dst_siteloc[2]]))
        bond = minmax(src_orig_color, dst_orig_color)
        i = findfirst(==(bond), orig_bondlist)
        isnothing(i) || (comp[ns + i] += 1)
    end

    return comp
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

"""
    AssemblySystem(polys::AbstractVector{<:Polyform}, bondlists::AbstractVector{<:AbstractVector{<:Integer}})

Construct a meta-assembly system where each polyform's allowed bonds are specified explicitly.

`bondlists[i]` is a vector of bond indices into `bonded_colors(orig_sys)` that polyform `i`
is permitted to participate in. An empty bond list makes polyform `i` passive: it exposes all
open binding sites but can only bind to active (non-empty bond list) polyforms. A non-empty
bond list makes polyform `i` active but restricts it to the specified bond types only; only
the corresponding binding sites are exposed as meta-sites.
"""
function AssemblySystem(polys::AbstractVector{<:Polyform}, bondlists::AbstractVector{<:AbstractVector{<:Integer}})
    length(polys) == length(bondlists) || throw(ArgumentError("polys and bondlists must have the same length"))
    isempty(polys) && throw(ArgumentError("polys must not be empty"))
    orig_sys = assemblysystem(first(polys))
    any(!=(orig_sys), assemblysystem(p) for p in polys) && throw(ArgumentError("polyforms do not come from the same assembly system"))
    orig_intmat = interactionmatrix(orig_sys)
    orig_bondcolors = bonded_colors(orig_sys)

    meta_species = map(polys, bondlists) do poly, bl
        if isempty(bl)
            MetaParticleSpecies(poly)
        else
            allowed_colors = Set{Int}(c for b in bl for c in orig_bondcolors[b])
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
                    color(site) ∈ allowed_colors || continue
                    push!(active_indices, k)
                end
            end
            MetaParticleSpecies(poly, active_indices)
        end
    end

    orig_colors = [[color(ms.meta_sites[k]) for k in 1:nsites(ms)] for ms in meta_species]

    bonds = NTuple{4,Int}[]
    for i in eachindex(meta_species), k in 1:nsites(meta_species[i])
        c_ik = orig_colors[i][k]
        for j in i:length(meta_species)
            isempty(bondlists[i]) && isempty(bondlists[j]) && continue
            lstart = j == i ? k : 1
            for l in lstart:nsites(meta_species[j])
                orig_intmat[c_ik, orig_colors[j][l]] || continue
                bidx = findfirst(==(minmax(c_ik, orig_colors[j][l])), orig_bondcolors)
                (isempty(bondlists[i]) || bidx ∈ bondlists[i]) || continue
                (isempty(bondlists[j]) || bidx ∈ bondlists[j]) || continue
                push!(bonds, (i, k, j, l))
            end
        end
    end

    return AssemblySystem(bonds, meta_species)
end

"""
    AssemblySystem(active_polys::AbstractVector{<:Polyform}, passive_polys::AbstractVector{<:Polyform})

Construct a meta-assembly system where active and passive polyforms are building blocks.

Active building blocks inherit interactions from the original assembly system and may bind
to both active and passive building blocks. Passive building blocks carry no direct
interactions and can only bind to active building blocks.
"""
function AssemblySystem(active_polys::AbstractVector{<:Polyform}, passive_polys::AbstractVector{<:Polyform})
    isempty(active_polys) && throw(ArgumentError("active_polys must not be empty"))
    orig_sys = assemblysystem(first(active_polys))
    for p in Iterators.flatten((active_polys, passive_polys))
        assemblysystem(p) !== orig_sys && throw(ArgumentError("polyforms do not come from the same assembly system"))
    end
    orig_intmat = interactionmatrix(orig_sys)

    n_active = length(active_polys)
    meta_species = [MetaParticleSpecies(poly) for poly in Iterators.flatten((active_polys, passive_polys))]

    orig_colors = [[color(ms.meta_sites[k]) for k in 1:nsites(ms)] for ms in meta_species]

    bonds = NTuple{4,Int}[]
    for i in eachindex(meta_species), k in 1:nsites(meta_species[i])
        c_ik = orig_colors[i][k]
        for j in i:length(meta_species)
            i > n_active && j > n_active && continue  # no passive-passive bonds
            lstart = j == i ? k : 1
            for l in lstart:nsites(meta_species[j])
                orig_intmat[c_ik, orig_colors[j][l]] && push!(bonds, (i, k, j, l))
            end
        end
    end

    return AssemblySystem(bonds, meta_species)
end

function make_bondslot_counter(metasys, bondslots)
    orig_sys = assemblysystem(Roly.species(metasys, 1).polyform)
    ns_meta = nspecies(metasys)
    nb_orig = nbonds(orig_sys)
    ns_orig = nspecies(orig_sys)

    particle_basebonds = [Roly.base_composition(Polyform(metasys, i))[ns_orig+1:end] for i in 1:ns_meta]
    function f(s, args...)
        metacomp = Roly.composition(s)
        basecomp = Roly.base_composition(s)
        slots = zeros(Int, nb_orig)
        bonds = basecomp[ns_orig+1:end]
        for (i, n) in pairs(metacomp[1:ns_meta])
            slots[bondslots[i]] .+= n
            bonds .-= n*particle_basebonds[i]
        end
        any(<(0), slots - bonds) && return false
        return true
    end

    return f
end