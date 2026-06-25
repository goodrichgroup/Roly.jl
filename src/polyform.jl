"""
    Polyform

`Polyform` is a struct containing all information required to instantiate a
polyform, which here is defined as an aggregate of building blocks that are connected
(bound) together at various binding sites. In a self-assembly context, a polyform is
often referred to a "self-assembled structure".
"""
mutable struct Polyform{D,P<:Particle,S<:AssemblySystem,G<:AbstractNautyGraph}
    graphrep::G
    sigma::Int
    canon2orig::Vector{Int}
    orig2canon::Vector{Int}
    particles::Vector{P}
    assemblysystem::S
end
function Polyform(sys::AssemblySystem{D}) where {D}
    P = posetype(sys)
    g = NautyDiGraph(0)
    return Polyform{D,Particle{P},typeof(sys),typeof(g)}(g, 1, Int[], Int[], Particle{P}[], sys)
end
function Polyform(sys::AssemblySystem{D}, i::Integer) where {D}
    P = posetype(sys)
    ps = species(sys, i)
    g = copy(graphrep(ps))
    canonize!(g)
    part = Particle(sys, i; leading_vertex=1)
    cvs = collect(vertices(g))
    return Polyform{D,Particle{P},typeof(sys),typeof(g)}(
        g, symmetrynumber(ps), cvs, invperm(cvs), [part], sys)
end

Base.copy(p::Polyform) = typeof(p)(copy(p.graphrep), p.sigma, copy(p.canon2orig),
                                   copy(p.orig2canon), copy(p.particles), p.assemblysystem)
function Base.copy!(dst::Polyform, src::Polyform)
    copy!(dst.graphrep, src.graphrep)
    dst.sigma = src.sigma
    copy!(dst.canon2orig, src.canon2orig)
    copy!(dst.orig2canon, src.orig2canon)
    copy!(dst.particles, src.particles)
    dst.assemblysystem = src.assemblysystem
    return dst
end

function Base.show(io::Core.IO, p::Polyform{D}) where {D}
    print(io, "Polyform{$D}[n=$(nparticles(p))]")
end
Base.show(io::Core.IO, ::Type{Polyform{D}}) where {D} = print(io, "Polyform{$D}")
Base.:(==)(p1::Polyform, p2::Polyform) = assemblysystem(p1) === assemblysystem(p2) && graphrep(p1) == graphrep(p2)

@inline function particle(p::Polyform, v::Integer)
    i = findfirst(pt -> pt.leading_vertex == v, p.particles)
    return isnothing(i) ? nothing : p.particles[i]
end
@inline nparticles(p::Polyform) = length(p.particles)
@inline nsites(p::Polyform) = sum(prt -> nsites(prt, p.assemblysystem), p.particles; init=0)
@inline graphrep(p::Polyform) = p.graphrep
@inline symmetrynumber(p::Polyform) = p.sigma
@inline assemblysystem(p::Polyform) = p.assemblysystem
@inline canonical_vertices(p::Polyform) = p.canon2orig
@inline is_leadingvertex(p::Polyform, v::Integer) = any(pt -> pt.leading_vertex == v, p.particles)


@inline dimension(::Polyform{D}) where {D} = D
@inline posetype(::Polyform{D,<:Particle{<:P}}) where {D,P} = P
@inline posetype(::Type{Polyform{D,<:Particle{<:P}}}) where {D,P} = P
@inline numtype(p::Polyform) = eltype(posetype(p))
@inline numtype(p::Type{Polyform}) = eltype(posetype(p))

# Convert between original (stable particle/site identifiers) and canonical (graph vertex) space.
@inline tocanon(p::Polyform, v::Integer) = p.orig2canon[v]
@inline toorig(p::Polyform, v::Integer) = p.canon2orig[v]

function _sync_orig2canonical!(poly::Polyform)
    resize!(poly.orig2canon, length(poly.canon2orig))
    for (i, v) in enumerate(poly.canon2orig)
        poly.orig2canon[v] = i
    end
    return
end

@inline interior_edges(p) = _filter_edges(p, Val(false))
@inline exterior_edges(p) = (e for e in _filter_edges(p, Val(true)) if e.src < e.dst)
function _filter_edges(p::Polyform, v::Val{exterior}) where exterior
    return _filter_edges(graphrep(p), v)
end
function _filter_edges(p::AbstractGraph, ::Val{exterior}) where exterior
    filtered_edges = Iterators.filter(edges(p)) do (; src, dst)
        isdouble = has_edge(p, dst, src)
        return exterior ? isdouble : !isdouble
    end
    return filtered_edges
end

function _vertex_to_particle_site(p::Polyform, v::Integer)
    sys = assemblysystem(p)
    v = toorig(p, v)
    for (i, part) in enumerate(p.particles)
        v in graphvertices(part, sys) || continue
        for j in 1:nsites(part, sys)
            v in bindingsites(part, sys, j).vertices && return (i, j)
        end
    end
    return nothing
end

"""
    bonds(p::Polyform)

Return a lazy iterator of named tuples `(p1, s1, p2, s2)` for each bond in `p`.
`p1`, `p2` are indices into `p.particles`; `s1`, `s2` are site indices within each particle.
"""
function bonds(p::Polyform)
    return (
        let ((i, j), (k, l)) = minmax(_vertex_to_particle_site(p, e.src), _vertex_to_particle_site(p, e.dst))
            (particle=i, site=j) => (particle=k, site=l)
        end
        for e in exterior_edges(p)
    )
end

"""
    bindingsites(p::Polyform, i::Integer)

Return the `i`th binding site of polyform `p`. The order of binding sites is
given by the isomorphism class of the graph encoding of `p` and is deterministic.
"""
function bindingsites(p::Polyform, i::Integer)
    sys = assemblysystem(p)
    k = 0
    for v in canonical_vertices(p)
        prtcl = particle(p, v)
        isnothing(prtcl) && continue

        for j in 1:nsites(prtcl, sys)
            k += 1
            if k == i
                return bindingsites(prtcl, sys, j)
            end
        end
    end
    return nothing
end
function bindingsites(p::Polyform)
    return (bindingsites(p, i) for i in 1:nsites(p))
end

# TODO This assumes that all particle species graphreps do not contain any double bonds.
# This is not true for 3d particles. However, we could enforce a special label for all "structural" vertices,
# so that we can additinally compare labels to make this work generally. This needs to be implemented before 3d
# particles work reliably.
function isbound_vertex(p::Polyform, v::Integer)
    outneighbor_bitvec = NautyGraphs.adjrow(graphrep(p), v)
    inneighbor_bitvec = NautyGraphs.adjcol(graphrep(p), v)
    return any(outneighbor_bitvec .* inneighbor_bitvec)

    # vspcs = label2species(assemblysystem(p), label(graphrep(p), v))
    # neighbor_bitvec = NautyGraphs.adjrow(graphrep(p), v)
    # return any(label2species(assemblysystem(p), label(graphrep(p), n)) != vspcs
    #            for n in eachindex(neighbor_bitvec) if neighbor_bitvec[n])
end

function _overlap_and_contacts(polyparticles::AbstractVector{<:Particle}, part::Particle, sys::AssemblySystem; allow_noninteracting=false, allow_misaligned=false, kwargs...)
    intmat = interactionmatrix(sys)

    contacts = NTuple{2,UnitRange{Int}}[]
    for polypart in polyparticles
        could_contact(polypart, part, sys) || continue
        overlap(polypart, part, sys) && return true, nothing

        for b1 in bindingsites(polypart, sys), b2 in bindingsites(part, sys)
            istouching(b1, b2) || continue

            interacting = intmat[color(b1), color(b2)]
            if !allow_noninteracting && !interacting
                return true, nothing
            end
            if !allow_misaligned && !isaligned(b1, b2)
                return true, nothing
            end
            push!(contacts, (b1.vertices, b2.vertices))
        end
    end
    return false, contacts
end
overlap_and_contacts(poly::Polyform, part::Particle; kwargs...) = _overlap_and_contacts(poly.particles, part, assemblysystem(poly); kwargs...)

function raise!(poly::Polyform, site::BindingSite, siteloc::BindingSiteLoc; kwargs...)
    sys = assemblysystem(poly)
    species_index, site_index = siteloc
    particle_species = species(sys, species_index)
    leading_vertex = nv(graphrep(poly)) + 1
    attached_particle = attach(particle_species, species_index, site_index, site; leading_vertex)

    has_overlap, contacting_vertices = overlap_and_contacts(poly, attached_particle; kwargs...)
    has_overlap && return missing

    push!(poly.particles, attached_particle)

    g_attach = graphrep(particle_species)
    add_vertices!(graphrep(poly); vertex_labels=labels(g_attach))
    for (; src, dst) in edges(g_attach)
        add_edge!(graphrep(poly), src + leading_vertex - 1, dst + leading_vertex - 1)
    end

    for (vs1, vs2) in contacting_vertices
        for (v1, v2) in zip(vs1, vs2)
            add_edge!(graphrep(poly), tocanon(poly, v1), v2)
            add_edge!(graphrep(poly), v2, tocanon(poly, v1))
        end
    end

    append!(poly.canon2orig, graphvertices(attached_particle, sys))
    perm, autg = nauty(graphrep(poly); canonize=true)
    poly.canon2orig .= @view poly.canon2orig[perm]
    _sync_orig2canonical!(poly)
    poly.sigma = convert(Int, autg.n)

    return poly
end

function lower!(poly::Polyform)
    n = nparticles(poly)
    if n == 0
        return nothing
    elseif n == 1
        pop!(poly.particles)
        resize!(poly.canon2orig, 0)
        resize!(poly.orig2canon, 0)
        rem_vertices!(graphrep(poly), vertices(graphrep(poly)))
        poly.sigma = 1
        return poly
    end

    sys = assemblysystem(poly)
    nv_g = nv(graphrep(poly))
    is_target = zeros(Bool, nv_g)
    forbidden = zeros(Bool, nv_g)
    explored = zeros(Bool, nv_g)
    queue = zeros(Cint, nv_g)

    for v in Iterators.reverse(canonical_vertices(poly))
        for k in v:-1:1
            is_leadingvertex(poly, k) || continue
            v = k
            break
        end

        part = particle(poly, v)
        vs0 = graphvertices(part, sys)
        vs = sort!(poly.orig2canon[vs0])

        are_cutvertices!(graphrep(poly), vs, is_target, forbidden, explored, queue) && continue

        lv = leading_vertex(part)
        idx = findfirst(pt -> pt.leading_vertex == lv, poly.particles)
        last_idx = lastindex(poly.particles)
        if idx < last_idx
            poly.particles[idx] = poly.particles[last_idx]
        end
        pop!(poly.particles)

        rem_vertices!(graphrep(poly), vs)
        deleteat!(poly.canon2orig, vs)

        for i in eachindex(poly.canon2orig)
            poly.canon2orig[i] -= searchsortedlast(vs0, poly.canon2orig[i])
        end

        perm, autg = nauty(graphrep(poly); canonize=true)
        poly.canon2orig .= @view perm[poly.canon2orig]
        _sync_orig2canonical!(poly)
        poly.sigma = convert(Int, autg.n)
        return poly
    end
    error("invariant violated: no removable particle found in connected polyform")
end

function collect_open_bindingsites!(sites, poly::Polyform)
    sys = assemblysystem(poly)
    empty!(sites)

    for orig_v in canonical_vertices(poly)
        part = particle(poly, orig_v)
        isnothing(part) && continue

        for k in 1:nsites(part, sys)
            site = bindingsites(part, sys, k)
            isbound_vertex(poly, tocanon(poly, first(site.vertices))) && continue
            isinert(sys, color(site)) && continue

            push!(sites, site)
        end
    end
    return sites
end
function collect_open_bindingsites(poly::Polyform)
    sys = assemblysystem(poly)
    sites = BindingSite{posetype(sys), numtype(sys)}[]
    return collect_open_bindingsites!(sites, poly)
end

function collect_compatible_pairs!(pairs, poly::Polyform)
    sys = assemblysystem(poly)
    empty!(pairs)

    for orig_v in canonical_vertices(poly)
        part = particle(poly, orig_v)
        isnothing(part) && continue

        for k in 1:nsites(part, sys)
            site = bindingsites(part, sys, k)
            isbound_vertex(poly, tocanon(poly, first(site.vertices))) && continue
            isinert(sys, color(site)) && continue

            for siteloc in compatible_sitelocs(sys, color(site))
                push!(pairs, (site, siteloc))
            end
        end
    end
    return pairs
end
function collect_compatible_pairs(poly::Polyform)
    sys = assemblysystem(poly)
    BS = BindingSite{posetype(sys), numtype(sys)}
    pairs = Tuple{BS,BindingSiteLoc}[]
    return collect_compatible_pairs!(pairs, poly)
end

function attach(ps::ParticleSpecies, species_index::Integer, site_index::Integer, at::BindingSite; leading_vertex)
    site = bindingsites(ps, site_index)
    particle_pose = standard_offset(at) * inv(site.pose)
    return Particle(particle_pose, leading_vertex, species_index)
end

"""
    composition(p::Polyform)

Return the composition vector of `p`: counts of each particle species followed by counts
of each bond type. Bond types are ordered as in `bonded_colors(assemblysystem(p))`, i.e.
by the sorted pair of colors of the connected binding sites.
"""
function composition(p::Polyform)
    sys = assemblysystem(p)
    ns = nspecies(sys)
    nb = nbonds(sys)
    comp = zeros(Int, ns + nb)

    for part in p.particles
        comp[part.species_index] += 1
    end

    g = graphrep(p)
    labs = labels(g)
    bondlist = bonded_colors(sys)

    for (; src, dst) in edges(g)
        src > dst && continue
        has_edge(g, dst, src) || continue
        bond = minmax(label2color(sys, labs[src]), label2color(sys, labs[dst]))
        i = findfirst(==(bond), bondlist)
        isnothing(i) || (comp[ns + i] += 1)
    end

    return comp
end
