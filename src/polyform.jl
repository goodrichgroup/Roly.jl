"""
    Polyform

A `Polyform` is an aggregate of particles connected at binding
sites, represented by a directed graph.
"""
mutable struct Polyform{D,P<:Particle,S<:BindingRules,G<:AbstractNautyGraph}
    graphrep::G
    sigma::Int
    canon2orig::Vector{Int}
    orig2canon::Vector{Int}
    particles::Vector{P}
    bindingrules::S
end

"""
    Polyform(rules::BindingRules{D}) where {D}

Create an empty polyform containing no particles.
"""
function Polyform(rules::BindingRules{D}) where {D}
    P = posetype(rules)
    g = NautyDiGraph(0)
    return Polyform{D,Particle{P},typeof(rules),typeof(g)}(g, 1, Int[], Int[], Particle{P}[], rules)
end

"""
    Polyform(rules::BindingRules{D}, i::Integer) where {D}

Create a single-particle polyform, consisting of species `i` of `rules`.
"""
function Polyform(rules::BindingRules{D}, i::Integer) where {D}
    P = posetype(rules)
    ps = species(rules, i)
    g = copy(graphrep(ps))
    part = Particle(rules, i; leadingvertex=1)
    perm = first(nauty(g; canonize=true))
    cvs = convert(Vector{Int}, perm)
    return Polyform{D,Particle{P},typeof(rules),typeof(g)}(g, symmetrynumber(ps), cvs, invperm(cvs), [part], rules)
end

function Base.copy(p::Polyform)
    return typeof(p)(
        copy(p.graphrep), p.sigma, copy(p.canon2orig), copy(p.orig2canon), copy(p.particles), p.bindingrules
    )
end
function Base.copy!(dst::Polyform, src::Polyform)
    copy!(dst.graphrep, src.graphrep)
    dst.sigma = src.sigma
    copy!(dst.canon2orig, src.canon2orig)
    copy!(dst.orig2canon, src.orig2canon)
    copy!(dst.particles, src.particles)
    dst.bindingrules = src.bindingrules
    return dst
end

Base.show(io::Core.IO, p::Polyform{D}) where {D} = print(io, "Polyform{$D}[n=$(nparticles(p))]")
Base.show(io::Core.IO, ::Type{Polyform{D}}) where {D} = print(io, "Polyform{$D}")
Base.:(==)(p1::Polyform, p2::Polyform) = bindingrules(p1) === bindingrules(p2) && graphrep(p1) == graphrep(p2)

"""
    nparticles(p::Polyform)

Return the number of particles in `p`.
"""
@inline nparticles(p::Polyform) = length(p.particles)

"""
    nsites(p::Polyform)

Return the total number of binding sites across all particles in `p`, including bound sites.
"""
@inline nsites(p::Polyform) = sum(prt -> nsites(prt, p.bindingrules), p.particles; init=0)

"""
    bindingrules(p::Polyform)

Return the `BindingRules` that `p` belongs to.
"""
@inline bindingrules(p::Polyform) = p.bindingrules

"""
    symmetrynumber(p::Polyform)

Return the symmetry number of `p`, i.e. the size of its automorphism group.
"""
@inline symmetrynumber(p::Polyform) = p.sigma

@inline graphrep(p::Polyform) = p.graphrep
@inline dimension(::Polyform{D}) where {D} = D
@inline posetype(::Polyform{D,<:Particle{<:P}}) where {D,P} = P
@inline posetype(::Type{<:Polyform{D,<:Particle{<:P}}}) where {D,P} = P
@inline numtype(p::Polyform) = eltype(posetype(p))
@inline numtype(::Type{<:Polyform{D,<:Particle{<:P}}}) where {D,P} = eltype(P)

"""
    tocanon(p::Polyform, v::Integer)

Convert an original vertex index `v` (as stored in `Particle.leadingvertex` or
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
    # `orig2canon` is rebuilt from the result three lines below, 
    # so we borrow it as scratch buffer
    resize!(poly.orig2canon, length(poly.canon2orig))
    @inbounds for i in eachindex(poly.orig2canon, perm)
        poly.orig2canon[i] = poly.canon2orig[perm[i]]
    end
    @inbounds for i in eachindex(poly.canon2orig, poly.orig2canon)
        poly.canon2orig[i] = poly.orig2canon[i]
    end
    @inbounds for i in eachindex(poly.canon2orig)
        poly.orig2canon[poly.canon2orig[i]] = i
    end
    return nothing
end

@inline is_leadingvertex(p::Polyform, v::Integer) = any(pt -> pt.leadingvertex == v, p.particles)

@inline particles(p::Polyform, i::Integer) = p.particles[i]

# Look up the particle whose leadingvertex equals the original vertex v (O(n) scan).
@inline function particle_from_leadingvertex(p::Polyform, v::Integer)
    i = findfirst(pt -> pt.leadingvertex == v, p.particles)
    return isnothing(i) ? nothing : p.particles[i]
end

"""
    interior_edges(p::Polyform)

Return a lazy iterator over the internal particle edges of `graphrep(p)`.
"""
@inline interior_edges(p::Polyform) = _filter_edges(p, Val(false))

"""
    exterior_edges(p::Polyform)

Return a lazy iterator over the external edges of `graphrep(p)`, those joining two particles.

A bond contributes *several* of these, one per vertex pair `contact_pairing` makes, so use
[`bonds`](@ref) to iterate bonds.
"""
@inline exterior_edges(p::Polyform) = (e for e in _filter_edges(p, Val(true)) if e.src < e.dst)

"""
    _same_particle(p::Polyform, u::Integer, v::Integer; canonidxs=true)

Return `true` if the graph vertices `u` and `v` belong to the same particle.

`canonidxs` says which numbering `u` and `v` are in.
"""
@inline function _same_particle(p::Polyform, u::Integer, v::Integer; canonidxs::Bool=true)
    # Each particle owns a contiguous block of original vertices starting at its leading vertex, so
    # `u` and `v` are split apart exactly when some leading vertex falls between them.
    canonidxs && ((u, v) = (toorig(p, u), toorig(p, v)))
    lo, hi = minmax(u, v)
    return !any(pt -> lo < leadingvertex(pt) <= hi, p.particles)
end

# An edge is a bond exactly when its endpoints belong to different particles
function _filter_edges(p::Polyform, ::Val{exterior}) where {exterior}
    return Iterators.filter(edges(graphrep(p))) do (; src, dst)
        same = _same_particle(p, src, dst)
        return exterior ? !same : same
    end
end

"""
    _isbound_vertex(p::Polyform, part::Particle, v::Integer; canonidxs=true)

Return `true` if the graph vertex `v` of particle `part` is bonded to another particle, i.e. if
it has a neighbor outside `part`'s own block of vertices.

`canonidxs` says which numbering `v` is in, the same way it does for [`bondindex`](@ref) and
[`_same_particle`](@ref).
"""
function _isbound_vertex(p::Polyform, part::Particle, v::Integer; canonidxs::Bool=true)
    own = graphvertices(part, bindingrules(p))
    neighs = NautyGraphs.adjrow(graphrep(p), canonidxs ? v : tocanon(p, v))
    for w in eachindex(neighs)
        neighs[w] || continue
        toorig(p, w) in own || return true
    end
    return false
end

"""
    bondindex(poly::Polyform, src::Integer, dst::Integer; canonidxs=true)

Return the index into `bonded_colors(bindingrules(poly))` for the bond between the graph
vertices `src` and `dst`, or `nothing` if they don't form a valid bond type.

`canonidxs` says which numbering `src` and `dst` are in: `true`, the default, for the canonical
one that `graphrep(poly)` and its `edges` are in, and `false` for the stable original one that
`BindingSite.vertices` and `Particle.leadingvertex` are in.
"""
function bondindex(poly::Polyform, src::Integer, dst::Integer; canonidxs::Bool=true)
    rules = bindingrules(poly)
    p1, b1 = _vertex_to_particle_site(poly, src; canonidxs)
    p2, b2 = _vertex_to_particle_site(poly, dst; canonidxs)
    c1 = color(bindingsite(particles(poly, p1), rules, b1))
    c2 = color(bindingsite(particles(poly, p2), rules, b2))
    return findfirst(==(minmax(c1, c2)), bonded_colors(rules))
end

# Map a graph vertex back to (particleindex, siteindex).
function _vertex_to_particle_site(p::Polyform, v::Integer; canonidxs::Bool=true)
    rules = bindingrules(p)
    orig_v = canonidxs ? toorig(p, v) : v
    for (i, part) in enumerate(p.particles)
        orig_v in graphvertices(part, rules) || continue
        for j in 1:nsites(part, rules)
            orig_v in bindingsite(part, rules, j).vertices && return (i, j)
        end
    end
    return nothing
end

"""
    bindingsite(p::Polyform, i::Integer)

Return the `i`-th binding site of `p`, counting through `p`'s particles in the order they are
stored and through each particle's own sites, exactly as on a [`ParticleSpecies`](@ref).

The ordering depends on how `p` was assembled. Use [`canonbindingsite`](@ref) for iterating through
binding sites in canonical order.
"""
function bindingsite(p::Polyform, i::Integer)
    rules = bindingrules(p)
    k = 0
    for prtcl in p.particles
        for j in 1:nsites(prtcl, rules)
            k += 1
            k == i && return bindingsite(prtcl, rules, j)
        end
    end
    return nothing
end

"""
    bindingsites(p::Polyform)

Return a lazy iterator over all binding sites of `p`, counting through `p`'s particles in the order they are
stored and through each particle's own sites, exactly as on a [`ParticleSpecies`](@ref).

The ordering depends on how `p` was assembled. Use [`canonbindingsite`](@ref) for iterating through
binding sites in canonical order.
"""
bindingsites(p::Polyform) = (bindingsite(p, i) for i in 1:nsites(p))

"""
    canonbindingsite(p::Polyform, i::Integer)

Return the `i`-th binding site of `p` in canonical order, which follows the canonical graph
labeling and is therefore the same for any two isomorphic polyforms.
"""
function canonbindingsite(p::Polyform, i::Integer)
    rules = bindingrules(p)
    k = 0
    for v in p.canon2orig
        prtcl = particle_from_leadingvertex(p, v)
        isnothing(prtcl) && continue
        for j in 1:nsites(prtcl, rules)
            k += 1
            k == i && return bindingsite(prtcl, rules, j)
        end
    end
    return nothing
end

"""
    canonbindingsites(p::Polyform)

Return a lazy iterator over all binding sites of `p` in canonical order.
"""
canonbindingsites(p::Polyform) = (canonbindingsite(p, i) for i in 1:nsites(p))

"""
    rotationcenter(p::Polyform)

The point every rotational symmetry of `p` fixes: the centroid of its particle positions.

A symmetry permutes the particles, so it fixes their centroid, and [`rotationgroup`](@ref
rotationgroup(::Polyform)) turns about this point rather than about `p`'s pose origin.
"""
rotationcenter(p::Polyform) = sum(pt.pose.x for pt in p.particles) / nparticles(p)

"""
    permutationgroup(p::Polyform)

Return the rotational symmetry group of `p` as particle permutations: one permutation per
rotation about the centroid of `p`'s particle positions that carries every particle onto a
particle of the same species, matching position and orientation.

The counterpart of [`rotationgroup(::Polyform)`](@ref), which returns the same group as
rotations, and the assembled analogue of [`permutationgroup(::ParticleSpecies)`](@ref), which
permutes a single particle's sites.

`length(permutationgroup(p))` should equal [`symmetrynumber`](@ref)`(p)`, which reads the same
order off the graph encoding. Use it to check an encoding rather than to compute the symmetry
number, which nauty gives far more cheaply.

The list holds one entry per rotation, so it may repeat: the action on particles is not
faithful. A straight chain of cubes turned about the chain axis leaves every particle where it
was and still moves the graph. Call `unique` on the result for the quotient group.
"""
function permutationgroup(p::Polyform)
    perms = Vector{Int}[]
    _eachpolyformsymmetry((_, perm) -> push!(perms, perm), p)
    return perms
end

"""
    rotationgroup(p::Polyform)

Return the rotational symmetry group of `p` as rotations about the centroid of `p`'s particle
positions: those carrying every particle onto a particle of the same species, matching position
and orientation.

The counterpart of [`permutationgroup(::Polyform)`](@ref), which returns the same group as
particle permutations.
"""
function rotationgroup(p::Polyform)
    rots = _rotationtype(posetype(p))[]
    _eachpolyformsymmetry((Q, _) -> push!(rots, Q), p)
    return rots
end

"""
    _eachpolyformsymmetry(f, p::Polyform)

Call `f(Q, perm)` on every rotation `Q` about the centroid of `p`'s particle positions that
carries each particle onto one of the same species, matching position and orientation, together
with the particle permutation `perm` it induces. Return how many there were.

A particle's frame need only match up to its species' own
[`rotationgroup`](@ref rotationgroup(::ParticleSpecies)), so the candidates are `Q = Rₐ S R₁⁻¹`
for each particle `a` of particle 1's species and each `S` in that group. Every element of the
group appears exactly once.
"""
function _eachpolyformsymmetry(f, p::Polyform)
    F = numtype(p)
    n = nparticles(p)
    # An empty polyform constrains nothing; report the identity, so the count still matches the
    # symmetry number its graph carries.
    n == 0 && (f(one(_rotationtype(posetype(p))), Int[]); return 1)

    rules = bindingrules(p)
    poses = [pt.pose for pt in p.particles]
    spcs = [speciesindex(pt) for pt in p.particles]
    own = Dict(s => rotationgroup(species(rules, s)) for s in unique(spcs))

    xs = [q.x - rotationcenter(p) for q in poses]
    tol = sqrt(eps(F))
    atol = tol * max(one(F), maximum(norm, xs))

    function permutation(Q)
        perm = zeros(Int, n)
        for i in 1:n
            j = findfirst(1:n) do j
                spcs[j] == spcs[i] &&
                    isapprox(Q * xs[i], xs[j]; atol) &&
                    any(S -> isapprox(Q * poses[i].psi, poses[j].psi * S; atol=tol), own[spcs[j]])
            end
            isnothing(j) && return nothing
            perm[i] = j
        end
        return perm
    end

    nfound = 0
    for a in 1:n
        spcs[a] == spcs[1] || continue
        for S in own[spcs[1]]
            Q = poses[a].psi * S * inv(poses[1].psi)
            perm = permutation(Q)
            isnothing(perm) && continue
            # A candidate is built from frames alone, and two particles of one species often
            # carry the same frame, so the same `Q` is reached once per such particle. Keeping
            # only the `a` that `Q` really sends particle 1 to picks each element out once.
            perm[1] == a || continue
            f(Q, perm)
            nfound += 1
        end
    end
    return nfound
end

"""
    _bondedges(p::Polyform)

Return one exterior edge per bond of `p`.

A bond joins two binding sites, and [`contact_pairing`](@ref) pairs `gcd(k₁, k₂)` of their
vertices, so a bond between two dart-encoded faces reaches `graphrep(p)` as several edges. Two
sites are joined by at most one bond, so the pair of sites an edge lands on names the bond.
"""
function _bondedges(p::Polyform)
    seen = Set{NTuple{2,NTuple{2,Int}}}()
    return filter(collect(exterior_edges(p))) do (; src, dst)
        key = minmax(_vertex_to_particle_site(p, src), _vertex_to_particle_site(p, dst))
        key in seen && return false
        push!(seen, key)
        return true
    end
end

"""
    bonds(p::Polyform)

Return a lazy iterator of bonds in `p`, each one once. Each element has the form

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
        end for e in _bondedges(p)
    )
end

"""
    nbonds(p::Polyform)

Return the number of bonds in `p`. Note that `nbonds(::BindingRules)` instead counts how many
*kinds* of bond a set of rules allows.
"""
nbonds(p::Polyform) = length(_bondedges(p))

"""
    composition(p::Polyform)

Return the composition vector of `p`: counts of each particle species (indices
`1:nspecies(rules)`) followed by counts of each bond type (indices `nspecies+1:end`).
Bond types are ordered as in `bonded_colors(bindingrules(p))`.
"""
function composition(p::Polyform)
    rules = bindingrules(p)
    ns = nspecies(rules)
    nb = nbonds(rules)
    comp = zeros(Int, ns + nb)

    for part in p.particles
        comp[part.speciesindex] += 1
    end

    for (; src, dst) in _bondedges(p)
        i = bondindex(p, src, dst)
        isnothing(i) || (comp[ns + i] += 1)
    end

    return comp
end

"""
    raise!(poly::Polyform, site::BindingSite, siteloc::BindingSiteLoc, t=0)

Attach a new particle to `poly` at the open binding site `site`, with the species and site index given by `siteloc`,
in twist `t` of the bond (see [`standard_twist`](@ref)).

Returns `poly` on success, or `missing` if the attachment is geometrically forbidden (overlap or misaligned contact).
"""
function raise!(poly::Polyform, site::BindingSite, siteloc::BindingSiteLoc, t::Integer=0; kwargs...)
    rules = bindingrules(poly)
    speciesindex, siteindex = siteloc
    particle_species = species(rules, speciesindex)
    leadingvertex = nv(graphrep(poly)) + 1
    mate = bindingsite(particle_species, siteindex)
    particle_pose = standard_twist(site, t, twistfreedom(site, mate)) * inv(mate.pose)
    attached_particle = Particle(particle_pose, leadingvertex, speciesindex)

    has_overlap, contacting_vertices = overlap_and_contacts(poly, attached_particle; kwargs...)
    has_overlap && return missing

    push!(poly.particles, attached_particle)

    g_attach = graphrep(particle_species)
    add_vertices!(graphrep(poly); vertex_labels=labels(g_attach))
    for (; src, dst) in edges(g_attach)
        add_edge!(graphrep(poly), src + leadingvertex - 1, dst + leadingvertex - 1)
    end

    # Every contact is paired in the twist it was found in, not the one this call
    for (vs1, vs2, twst, ntwists) in contacting_vertices
        for (v1, v2) in contact_pairing(vs1, vs2, twst, ntwists)
            add_edge!(graphrep(poly), tocanon(poly, v1), v2)
            add_edge!(graphrep(poly), v2, tocanon(poly, v1))
        end
    end

    append!(poly.canon2orig, graphvertices(attached_particle, rules))
    perm, autg = nauty(graphrep(poly); canonize=true)
    _apply_perm!(poly, perm)
    poly.sigma = convert(Int, autg.n)
    return poly
end

"""
    lower!(poly::Polyform)

Remove the last particle from `poly`, following canonical ordering.

Returns `poly` on success, `nothing` if `poly` is already empty.
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

    rules = bindingrules(poly)
    nv_g = nv(graphrep(poly))
    target = zeros(Bool, nv_g)
    dist = zeros(Int, nv_g)
    queue = zeros(Int, nv_g)

    # A particle owns >=1 vertex per site, so the scan will revisit vertices many times.
    # -> store for performance.
    tested = Int[]
    part = nothing
    for v in Iterators.reverse(poly.canon2orig)
        # Walk backward from v to find the leading vertex of its particle.
        for k in v:-1:1
            is_leadingvertex(poly, k) || continue
            v = k
            break
        end
        v in tested && continue
        push!(tested, v)

        part = particle_from_leadingvertex(poly, v)
        is_cutset(graphrep(poly), @view(poly.orig2canon[graphvertices(part, rules)]); target, dist, queue) || break
        part = nothing
    end
    isnothing(part) && error("Internal error: no removable particle found in connected polyform. Please file an issue.")
    return _remove_particle!(poly, part)
end

"""
    _remove_particle!(poly::Polyform, part::Particle)

Remove particle `part` from `poly`. Does not check for connectedness before.
"""
function _remove_particle!(poly::Polyform, part::Particle)
    rules = bindingrules(poly)
    vs0 = graphvertices(part, rules)
    vs = sort!(poly.orig2canon[vs0])

    lv = leadingvertex(part)
    idx = findfirst(pt -> pt.leadingvertex == lv, poly.particles)
    last_idx = lastindex(poly.particles)
    idx < last_idx && (poly.particles[idx] = poly.particles[last_idx])
    pop!(poly.particles)

    rem_vertices!(graphrep(poly), vs)
    deleteat!(poly.canon2orig, vs)

    # shift remaining indices and leading vertices
    for i in eachindex(poly.canon2orig)
        poly.canon2orig[i] > last(vs0) && (poly.canon2orig[i] -= length(vs0))
    end
    for i in eachindex(poly.particles)
        leadingvertex(poly.particles[i]) > last(vs0) || continue
        poly.particles[i] = shift_leadingvertex(poly.particles[i], -length(vs0))
    end

    perm, autg = nauty(graphrep(poly); canonize=true)
    _apply_perm!(poly, perm)
    poly.sigma = convert(Int, autg.n)
    return poly
end

"""
    subpolyform(poly::Polyform, particleids)

Return the sub-polyform induced by the particles `particleids`, renumbered in the given order.

The particle subset is expected to be connected; bonds to removed particles are dropped and their
sites become unbound.
"""
function subpolyform(poly::Polyform, particleids)
    rules = bindingrules(poly)
    newparticles = eltype(poly.particles)[]
    verts = Int[]
    lv = 1
    for p in particleids
        part = particles(poly, p)
        for ov in graphvertices(part, rules)
            push!(verts, tocanon(poly, ov))
        end
        push!(newparticles, typeof(part)(part.pose, lv, part.speciesindex))
        lv += nv(graphrep(species(rules, part.speciesindex)))
    end

    g = graphrep(poly)[verts]
    # `g` starts out in new-original vertex order, so the canon maps are the permutation itself.
    perm, autg = nauty(g; canonize=true)
    return typeof(poly)(g, convert(Int, autg.n), collect(Int, perm), invperm(perm), newparticles, rules)
end

function _overlap_and_contacts(
    polyparticles::AbstractVector{<:Particle},
    part::Particle,
    rules::BindingRules;
    allow_noninteracting=false,
    allow_misaligned=false,
    kwargs...,
)
    intmat = interactionmatrix(rules)
    contacts = Tuple{UnitRange{Int},UnitRange{Int},Int,Int}[]

    for polypart in polyparticles
        could_contact(polypart, part, rules) || continue
        overlap(polypart, part, rules) && return true, nothing

        for b1 in bindingsites(polypart, rules), b2 in bindingsites(part, rules)
            istouching(b1, b2) || continue
            interacting = intmat[color(b1), color(b2)]
            !allow_noninteracting && !interacting && return true, nothing
            twst = twist(b1, b2)
            !allow_misaligned && isnothing(twst) && return true, nothing
            push!(contacts, (b1.vertices, b2.vertices, something(twst, 0), twistfreedom(b1, b2)))
        end
    end
    return false, contacts
end

function overlap_and_contacts(poly::Polyform, part::Particle; kwargs...)
    return _overlap_and_contacts(poly.particles, part, bindingrules(poly); kwargs...)
end

# Walk the unbound binding sites of `poly` in canonical order, yielding each one's
# `(particle, site)` location together with the site itself. The four accessors below project it.
function _exposed(poly::Polyform)
    rules = bindingrules(poly)
    index = Dict(leadingvertex(p) => i for (i, p) in enumerate(poly.particles))
    out = Tuple{NTuple{2,Int},BindingSite{posetype(rules),numtype(rules)}}[]
    for orig_v in poly.canon2orig
        part = particle_from_leadingvertex(poly, orig_v)
        isnothing(part) && continue
        for k in 1:nsites(part, rules)
            s = bindingsite(part, rules, k)
            _isbound_vertex(poly, part, first(s.vertices); canonidxs=false) && continue
            push!(out, ((index[leadingvertex(part)], k), s))
        end
    end
    return out
end

"""
    exposedsitelocs(poly::Polyform)

The `(particle, site)` location of every *unbound* binding site of `poly`, in canonical order.

Bound sites are consumed by the bonds holding `poly` together and are never listed. Sites whose
color takes part in no rule are listed, even though nothing can attach through them as `poly`
stands; [`opensitelocs`](@ref) is this list without them.

These are the sites a [`MetaParticleSpecies`](@ref) may expose. It exposes the open ones by
default, and an inert one becomes usable simply by being named and given a live color.

A location here is `(particle, site)`, not the `(species, site)` of [`BindingSiteLoc`](@ref).
"""
exposedsitelocs(poly::Polyform) = [l for (l, _) in _exposed(poly)]

"""
    exposedsites(poly::Polyform)

Every unbound binding site of `poly`, in canonical order. See [`exposedsitelocs`](@ref).
"""
exposedsites(poly::Polyform) = [s for (_, s) in _exposed(poly)]

"""
    opensitelocs(poly::Polyform)

The `(particle, site)` location of every binding site of `poly` a partner can still attach
through: the unbound ones whose color some rule uses, in canonical order.

See [`exposedsitelocs`](@ref), which lists the inert ones too.
"""
function opensitelocs(poly::Polyform)
    rules = bindingrules(poly)
    return [l for (l, s) in _exposed(poly) if !isinert(rules, color(s))]
end

"""
    opensites(poly::Polyform)

Every binding site of `poly` a partner can still attach through, in canonical order. See
[`opensitelocs`](@ref).
"""
function opensites(poly::Polyform)
    rules = bindingrules(poly)
    return [s for (_, s) in _exposed(poly) if !isinert(rules, color(s))]
end

"""
    _deletable_species(poly; target, dist, queue)

Return `(top, runnerup, top_lv)`: the two largest species indices among the particles `lower!`
could delete from `poly`, and the top leading vertex of the particle achieving `top`.

`lower!` deletes from the highest label class holding a removable particle. Canonical position
runs with vertex label (nauty orders classes by color, `vertexlabels2labptn` by label), and
`_adjust_labels_and_colors` puts species `s`'s labels above species `s-1`'s, so species index
order is label order. Removability is read off `poly` rather than the child: attaching a
particle cannot disconnect what removing a different one leaves behind.

The runner-up lets the anchor particle exclude itself in O(1); see
[`collect_attachments!`](@ref).
"""
function _deletable_species(poly::Polyform; target, dist, queue)
    rules = bindingrules(poly)
    n = nparticles(poly)
    top, runnerup, top_lv = 0, 0, 0
    for part in poly.particles
        n > 1 &&
            is_cutset(graphrep(poly), @view(poly.orig2canon[graphvertices(part, rules)]); target, dist, queue) &&
            continue
        s = speciesindex(part)
        if s > top
            top, runnerup, top_lv = s, top, leadingvertex(part)
        elseif s > runnerup
            runnerup = s
        end
    end
    return top, runnerup, top_lv
end

"""
    collect_attachments!(attachments, poly::Polyform)

Fill `attachments` with every way of growing `poly` by one particle, as triples
`(site, siteloc, t)`: an open binding site of `poly`, the species and site index of the
particle to attach, and which of the bond's distinct twists to attach it in.

Two filters keep the list down to what reverse search can use. Only one partner site per
symmetry orbit is kept, since the others give the same child. And a candidate is dropped when
`lower!` would not undo it, i.e. when the child holds a removable particle from a species
index above the one being attached: the child's parent would then be a different polyform and
reverse search would reject the pair anyway.

The second filter is conservative about the anchor particle, the one carrying `site`. The
anchor is excluded from the removable set, because a new particle bonded to it alone leaves it
a cut vertex of the child. When the attachment also closes a ring the anchor does stay
removable, and excluding it then only lets a few extra candidates through, which reverse search
rejects.
"""
function collect_attachments!(attachments, poly::Polyform)
    rules = bindingrules(poly)
    empty!(attachments)
    nv_g = nv(graphrep(poly))
    # `dist` carries distances and the -1/-2 markers `is_cutset` needs, so it cannot be a Bool
    target, dist, queue = zeros(Bool, nv_g), zeros(Int, nv_g), zeros(Int, nv_g)

    top, runnerup, top_lv = _deletable_species(poly; target, dist, queue)
    for orig_v in poly.canon2orig
        anchor = particle_from_leadingvertex(poly, orig_v)
        isnothing(anchor) && continue

        # The anchor is excluded, so it takes the runner-up when it is itself the top scorer.
        deletable = leadingvertex(anchor) == top_lv ? runnerup : top
        for k in 1:nsites(anchor, rules)
            site = bindingsite(anchor, rules, k)
            _isbound_vertex(poly, anchor, first(site.vertices); canonidxs=false) && continue
            isinert(rules, color(site)) && continue

            for siteloc in distinct_attachments(rules, color(site))
                # The new particle is always removable, so `lower!` stops at its label class
                # unless a higher one survives.
                siteloc[1] < deletable && continue
                mate = bindingsite(species(rules, siteloc[1]), siteloc[2])
                for t in 0:(_ndistincttwists(site, mate) - 1)
                    push!(attachments, (site, siteloc, t))
                end
            end
        end
    end
    return attachments
end

"""
    collect_attachments(poly::Polyform)

Return every way of growing `poly` by one particle, allocating the vector.
See [`collect_attachments!`](@ref).
"""
function collect_attachments(poly::Polyform)
    rules = bindingrules(poly)
    BS = BindingSite{posetype(rules),numtype(rules)}
    return collect_attachments!(Tuple{BS,BindingSiteLoc,Int}[], poly)
end
