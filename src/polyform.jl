"""
    Polyform

A `Polyform` is an aggregate of particles connected at binding
sites, represented by a directed graph.
"""
mutable struct Polyform{D,P<:Particle,S<:AssemblySystem,G<:AbstractNautyGraph}
    graphrep::G
    sigma::Int
    canon2orig::Vector{Int}
    orig2canon::Vector{Int}
    particles::Vector{P}
    assemblysystem::S
end

"""
    Polyform(sys::AssemblySystem{D}) where {D}

Create an empty polyform containing no particles.
"""
function Polyform(sys::AssemblySystem{D}) where {D}
    P = posetype(sys)
    g = NautyDiGraph(0)
    return Polyform{D,Particle{P},typeof(sys),typeof(g)}(g, 1, Int[], Int[], Particle{P}[], sys)
end

"""
    Polyform(sys::AssemblySystem{D}, i::Integer) where {D}

Create a single-particle polyform, consisting of species `i` of `sys`.
"""
function Polyform(sys::AssemblySystem{D}, i::Integer) where {D}
    P = posetype(sys)
    ps = species(sys, i)
    g = copy(graphrep(ps))
    canonize!(g)
    part = Particle(sys, i; leading_vertex=1)
    cvs = collect(vertices(g))
    return Polyform{D,Particle{P},typeof(sys),typeof(g)}(g, symmetrynumber(ps), cvs, invperm(cvs), [part], sys)
end

function Base.copy(p::Polyform)
    return typeof(p)(
        copy(p.graphrep), p.sigma, copy(p.canon2orig), copy(p.orig2canon), copy(p.particles), p.assemblysystem
    )
end
function Base.copy!(dst::Polyform, src::Polyform)
    copy!(dst.graphrep, src.graphrep)
    dst.sigma = src.sigma
    copy!(dst.canon2orig, src.canon2orig)
    copy!(dst.orig2canon, src.orig2canon)
    copy!(dst.particles, src.particles)
    dst.assemblysystem = src.assemblysystem
    return dst
end

Base.show(io::Core.IO, p::Polyform{D}) where {D} = print(io, "Polyform{$D}[n=$(nparticles(p))]")
Base.show(io::Core.IO, ::Type{Polyform{D}}) where {D} = print(io, "Polyform{$D}")
Base.:(==)(p1::Polyform, p2::Polyform) = assemblysystem(p1) === assemblysystem(p2) && graphrep(p1) == graphrep(p2)

"""
    nparticles(p::Polyform)

Return the number of particles in `p`.
"""
@inline nparticles(p::Polyform) = length(p.particles)

"""
    nsites(p::Polyform)

Return the total number of binding sites across all particles in `p`, including bound sites.
"""
@inline nsites(p::Polyform) = sum(prt -> nsites(prt, p.assemblysystem), p.particles; init=0)

"""
    assemblysystem(p::Polyform)

Return the `AssemblySystem` that `p` belongs to.
"""
@inline assemblysystem(p::Polyform) = p.assemblysystem

"""
    symmetrynumber(p::Polyform)

Return the symmetry number of `p`, i.e. the size of its automorphism group.
"""
@inline symmetrynumber(p::Polyform) = p.sigma

@inline graphrep(p::Polyform) = p.graphrep
@inline dimension(::Polyform{D}) where {D} = D
@inline posetype(::Polyform{D,<:Particle{<:P}}) where {D,P} = P
@inline posetype(::Type{Polyform{D,<:Particle{<:P}}}) where {D,P} = P
@inline numtype(p::Polyform) = eltype(posetype(p))
@inline numtype(p::Type{Polyform}) = eltype(posetype(p))

"""
    tocanon(p::Polyform, v::Integer)

Convert an original vertex index `v` (as stored in `Particle.leading_vertex` or
`BindingSite.vertices`) to the corresponding canonical graph vertex index in `graphrep(p)`.
"""
@inline tocanon(p::Polyform, v::Integer) = p.orig2canon[v]

"""
    toorig(p::Polyform, v::Integer)

Convert a canonical graph vertex index `v` (as returned by iterating `graphrep(p)`)
to the corresponding stable original vertex index.
"""
@inline toorig(p::Polyform, v::Integer) = p.canon2orig[v]

function _apply_perm!(poly::Polyform, perm)
    poly.canon2orig .= @view poly.canon2orig[perm]
    resize!(poly.orig2canon, length(poly.canon2orig))
    for i in eachindex(poly.canon2orig)
        poly.orig2canon[poly.canon2orig[i]] = i
    end
    return nothing
end

@inline canonical_vertices(p::Polyform) = p.canon2orig
@inline is_leadingvertex(p::Polyform, v::Integer) = any(pt -> pt.leading_vertex == v, p.particles)

@inline particles(p::Polyform, i::Integer) = p.particles[i]

# Look up the particle whose leading_vertex equals the original vertex v (O(n) scan).
@inline function particle_from_leadingvertex(p::Polyform, v::Integer)
    i = findfirst(pt -> pt.leading_vertex == v, p.particles)
    return isnothing(i) ? nothing : p.particles[i]
end

"""
    interior_edges(p::Polyform)

Return a lazy iterator over the internal particle edges of `graphrep(p)`.
"""
@inline interior_edges(p::Polyform) = _filter_edges(graphrep(p), Val(false))

"""
    exterior_edges(p::Polyform)

Return a lazy iterator over the external edges of `graphrep(p)`, corresponding to bonds in `p`.
"""
@inline exterior_edges(p::Polyform) = (e for e in _filter_edges(graphrep(p), Val(true)) if e.src < e.dst)

# TODO: assumes no double edges within a single particle's graphrep (not true for 3D particles).
# A label-based check is needed before 3D particles work reliably.
function _filter_edges(g::AbstractGraph, ::Val{exterior}) where {exterior}
    return Iterators.filter(edges(g)) do (; src, dst)
        isdouble = has_edge(g, dst, src)
        return exterior ? isdouble : !isdouble
    end
end

function isbound_vertex(p::Polyform, v::Integer)
    return any(NautyGraphs.adjrow(graphrep(p), v) .* NautyGraphs.adjcol(graphrep(p), v))
end

"""
    bondindex(poly::Polyform, src::Integer, dst::Integer)

Return the index into `bonded_colors(assemblysystem(poly))` for the bond between
canonical graph vertices `src` and `dst`, or `nothing` if they don't form a valid
bond type.
"""
function bondindex(poly::Polyform, src::Integer, dst::Integer)
    sys = assemblysystem(poly)
    p1, b1 = _vertex_to_particle_site(poly, src)
    p2, b2 = _vertex_to_particle_site(poly, dst)
    c1 = color(bindingsites(particles(poly, p1), sys, b1))
    c2 = color(bindingsites(particles(poly, p2), sys, b2))
    return findfirst(==(minmax(c1, c2)), bonded_colors(sys))
end

# Map a canonical graph vertex back to (particle_index, site_index).
function _vertex_to_particle_site(p::Polyform, v::Integer)
    sys = assemblysystem(p)
    orig_v = toorig(p, v)
    for (i, part) in enumerate(p.particles)
        orig_v in graphvertices(part, sys) || continue
        for j in 1:nsites(part, sys)
            orig_v in bindingsites(part, sys, j).vertices && return (i, j)
        end
    end
    return nothing
end

"""
    bindingsites(p::Polyform, i::Integer)

Return the `i`-th binding site of `p`. The ordering is determined by the canonical
graph labeling and is therefore deterministic across equivalent polyforms.
"""
function bindingsites(p::Polyform, i::Integer)
    sys = assemblysystem(p)
    k = 0
    for v in canonical_vertices(p)
        prtcl = particle_from_leadingvertex(p, v)
        isnothing(prtcl) && continue
        for j in 1:nsites(prtcl, sys)
            k += 1
            k == i && return bindingsites(prtcl, sys, j)
        end
    end
    return nothing
end

"""
    bindingsites(p::Polyform)

Return a lazy iterator over all binding sites of `p` in canonical order.
"""
bindingsites(p::Polyform) = (bindingsites(p, i) for i in 1:nsites(p))

"""
    bonds(p::Polyform)

Return a lazy iterator of bonds in `p`. Each element has the form

    (particle=i, site=j) => (particle=k, site=l)

where `i`, `k` are indices into `p.particles` and `j`, `l` are site indices
within those particles.
"""
function bonds(p::Polyform)
    return (
        let (lhs, rhs) = minmax(_vertex_to_particle_site(p, e.src), _vertex_to_particle_site(p, e.dst)),
            (p1, s1) = lhs,
            (p2, s2) = rhs

            (particle=p1, site=s1) => (particle=p2, site=s2)
        end for e in exterior_edges(p)
    )
end

"""
    composition(p::Polyform)

Return the composition vector of `p`: counts of each particle species (indices
`1:nspecies(sys)`) followed by counts of each bond type (indices `nspecies+1:end`).
Bond types are ordered as in `bonded_colors(assemblysystem(p))`.
"""
function composition(p::Polyform)
    sys = assemblysystem(p)
    ns = nspecies(sys)
    nb = nbonds(sys)
    comp = zeros(Int, ns + nb)

    for part in p.particles
        comp[part.species_index] += 1
    end

    for (; src, dst) in exterior_edges(p)
        i = bondindex(p, src, dst)
        isnothing(i) || (comp[ns + i] += 1)
    end

    return comp
end


function _overlap_and_contacts(
    polyparticles::AbstractVector{<:Particle},
    part::Particle,
    sys::AssemblySystem;
    allow_noninteracting=false,
    allow_misaligned=false,
    kwargs...,
)
    intmat = interactionmatrix(sys)
    contacts = NTuple{2,UnitRange{Int}}[]

    for polypart in polyparticles
        could_contact(polypart, part, sys) || continue
        overlap(polypart, part, sys) && return true, nothing

        for b1 in bindingsites(polypart, sys), b2 in bindingsites(part, sys)
            istouching(b1, b2) || continue
            interacting = intmat[color(b1), color(b2)]
            !allow_noninteracting && !interacting && return true, nothing
            !allow_misaligned && !isaligned(b1, b2) && return true, nothing
            push!(contacts, (b1.vertices, b2.vertices))
        end
    end
    return false, contacts
end

function overlap_and_contacts(poly::Polyform, part::Particle; kwargs...)
    return _overlap_and_contacts(poly.particles, part, assemblysystem(poly); kwargs...)
end

"""
    raise!(poly::Polyform, site::BindingSite, siteloc::BindingSiteLoc)

Attach a new particle to `poly` at the open binding site `site`, with the species
and site index given by `siteloc`. Returns `poly` on success, or `missing` if the
attachment is geometrically forbidden (overlap or misaligned contact).
"""
function raise!(poly::Polyform, site::BindingSite, siteloc::BindingSiteLoc; kwargs...)
    sys = assemblysystem(poly)
    species_index, site_index = siteloc
    particle_species = species(sys, species_index)
    leading_vertex = nv(graphrep(poly)) + 1
    particle_pose = standard_offset(site) * inv(bindingsites(particle_species, site_index).pose)
    attached_particle = Particle(particle_pose, leading_vertex, species_index)

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
    _apply_perm!(poly, perm)
    poly.sigma = convert(Int, autg.n)

    return poly
end

"""
    lower!(poly::Polyform)

Remove one particle from `poly`, choosing the last removable particle in
reverse canonical order (i.e. the particle whose removal keeps the polyform
connected). Returns `poly` on success, `nothing` if `poly` is already empty.

Mutates `poly` in place and re-canonicalizes the graph.
"""
function lower!(poly::Polyform)
    n = nparticles(poly)
    n == 0 && return nothing
    if n == 1
        pop!(poly.particles)
        resize!(poly.canon2orig, 0)
        resize!(poly.orig2canon, 0)
        rem_vertices!(graphrep(poly), vertices(graphrep(poly)))
        poly.sigma = 1
        return poly
    end

    sys = assemblysystem(poly)
    nv_g = nv(graphrep(poly))
    target = zeros(Bool, nv_g)
    visited = zeros(Bool, nv_g)
    queue = zeros(Cint, nv_g)

    for v in Iterators.reverse(canonical_vertices(poly))
        # Walk backward from v to find the leading vertex of its particle.
        for k in v:-1:1
            is_leadingvertex(poly, k) || continue
            v = k
            break
        end

        part = particle_from_leadingvertex(poly, v)
        is_cutset(graphrep(poly), @view(poly.orig2canon[graphvertices(part, sys)]); target, visited, queue) && continue

        return _remove_particle!(poly, part)
    end
    return error("invariant violated: no removable particle found in connected polyform")
end

"""
    _remove_particle!(poly::Polyform, part::Particle)

Remove particle `part` from `poly`. Does not check for connectedness before.
"""
function _remove_particle!(poly::Polyform, part::Particle)
    sys = assemblysystem(poly)
    vs0 = graphvertices(part, sys)
    vs = sort!(poly.orig2canon[vs0])

    lv = leading_vertex(part)
    idx = findfirst(pt -> pt.leading_vertex == lv, poly.particles)
    last_idx = lastindex(poly.particles)
    idx < last_idx && (poly.particles[idx] = poly.particles[last_idx])
    pop!(poly.particles)

    rem_vertices!(graphrep(poly), vs)
    deleteat!(poly.canon2orig, vs)

    for i in eachindex(poly.canon2orig)
        poly.canon2orig[i] -= searchsortedlast(vs0, poly.canon2orig[i])
    end

    perm, autg = nauty(graphrep(poly); canonize=true)
    _apply_perm!(poly, perm)
    poly.sigma = convert(Int, autg.n)
    return poly
end

function collect_open_bindingsites!(sites, poly::Polyform)
    sys = assemblysystem(poly)
    empty!(sites)
    for orig_v in canonical_vertices(poly)
        part = particle_from_leadingvertex(poly, orig_v)
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
    sites = BindingSite{posetype(sys),numtype(sys)}[]
    return collect_open_bindingsites!(sites, poly)
end

function collect_compatible_pairs!(pairs, poly::Polyform)
    sys = assemblysystem(poly)
    empty!(pairs)
    for orig_v in canonical_vertices(poly)
        part = particle_from_leadingvertex(poly, orig_v)
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
    BS = BindingSite{posetype(sys),numtype(sys)}
    pairs = Tuple{BS,BindingSiteLoc}[]
    return collect_compatible_pairs!(pairs, poly)
end
