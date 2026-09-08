"""
An enumeration of tilings by reverse search, run alongside the one in `tiling.jl` until it is
trusted. Nothing here is wired into the exported API.

A state is a polyform together with a set of its open sites identified in pairs, which is the same
thing a [`Tiling`](@ref) is. Two moves reach every state: attaching a particle, which is `raise!`,
and identifying a pair of sites, which makes a periodic bond. Identifying comes after attaching,
because a particle cannot be removed while a bond refers to its sites.
"""
struct RSClose{V}
    a::Int          # first vertex of one site, in the cell's original numbering
    b::Int          # and of the other
    t::V            # the translation carrying the second onto the first
end

mutable struct RSTiling{P<:Polyform,G<:AbstractNautyGraph,V}
    cell::P
    closes::Vector{RSClose{V}}   # the identifications made, in the order they were made
    key::G                       # the cell's graph with a vertex per bond, canonized
    markerat::Vector{Int}        # where each close's vertex ended up in the canonical order
end

function RSTiling(cell::Polyform)
    V = SVector{dimension(bindingrules(cell)),numtype(bindingrules(cell))}
    s = RSTiling(cell, RSClose{V}[], copy(graphrep(cell)), Int[])
    return _rekey!(s)
end

Base.copy(s::RSTiling) = RSTiling(copy(s.cell), copy(s.closes), copy(s.key), copy(s.markerat))
function Base.copy!(dst::RSTiling, src::RSTiling)
    copy!(dst.cell, src.cell)
    copy!(dst.closes, src.closes)
    copy!(dst.key, src.key)
    copy!(dst.markerat, src.markerat)
    return dst
end

# Two states are the same when the structures they describe are, which is when their graphs are.
_samestate(a::RSTiling, b::RSTiling) = a.key == b.key

function Base.show(io::Core.IO, s::RSTiling)
    return print(io, "RSTiling[n=", nparticles(s.cell), ", closes=", length(s.closes), "]")
end

# Rebuild the key from the cell and the closes: the cell's own graph, plus a vertex per bond
# joined to the two sites it identifies, canonized. Records where each close's vertex landed, so
# that the parent can pick one of them in an order that does not depend on how the state was
# reached.
function _rekey!(s::RSTiling)
    rules = bindingrules(s.cell)
    g = NautyDiGraph(0)
    for part in s.cell.particles
        blockdiag!(g, graphrep(species(rules, speciesindex(part))))
    end

    # every bond is marked, the cell's own alongside the ones an identification makes. Marking only
    # the latter would record where the structure was cut, and two cuts of one tiling would then
    # look like two tilings
    marker = _markerlabel(rules)
    for e in exterior_edges(s.cell)
        _addmarker!(g, marker, toorig(s.cell, e.src), toorig(s.cell, e.dst))
    end
    first_of = Int[]
    for c in s.closes
        push!(first_of, nv(g) + 1)
        for (v1, v2) in contact_pairing(_contactof(s.cell, c))
            _addmarker!(g, marker, v1, v2)
        end
    end
    perm, _ = nauty(g; canonize=true)
    place = invperm(collect(Int, perm))
    s.key = g
    resize!(s.markerat, length(first_of))
    s.markerat .= (place[v] for v in first_of)
    return s
end

# The contact a close stands for: the two sites meet once one is carried onto the other, so the
# bond they make is the one they would make side by side.
function _contactof(cell::Polyform, c::RSClose)
    rules = bindingrules(cell)
    a = bindingsite(cell, _vertex_to_particle_site(cell, c.a; canonidxs=false))
    b = bindingsite(cell, _vertex_to_particle_site(cell, c.b; canonidxs=false))
    moved = translate(b, c.t)
    return Contact(a.vertices, b.vertices, twist(a, moved), twistfreedom(a, moved))
end

### the structure a state describes

# How wide the cell is, and how far a particle reaches past its own center.
function _rsreach(cell::Polyform)
    rules = bindingrules(cell)
    parts = cell.particles
    diameter = maximum(norm(p.pose.x - q.pose.x) for p in parts, q in parts)
    radius = maximum(bounding_radius(species(rules, speciesindex(p))) for p in parts)
    return diameter + 2radius, radius
end

_rsnear(parts, t, radius) =
    let d = norm(t)
        d > 0 || return false
        u = t / d
        lo, hi = extrema(dot(p.pose.x, u) for p in parts)
        return d <= (hi - lo) + 2radius + sqrt(eps(d))
    end

# A position as a dictionary key. Adding zero folds `-0.0` onto `0.0`, which are not `isequal` and
# so would otherwise land in different buckets, and a rounded position is exact enough to compare.
_rskey(v) = Tuple(round.(v; digits=7) .+ 0.0)

# The lattice points near enough to matter, walked outward from the origin along the translations
# the bonds name and their inverses. No basis is involved, which is the point: a set of bonds
# generates its lattice whether or not any subset of them is a basis of it.
function _rspoints(gens, span)
    V = eltype(gens)
    pts = V[zero(V)]
    isempty(gens) && return pts
    seen = Set{typeof(_rskey(zero(V)))}([_rskey(zero(V))])
    i = 1
    while i <= length(pts)
        p = pts[i]
        i += 1
        for g in gens, sgn in (1, -1)
            q = p + sgn * g
            norm(q) <= 2span || continue
            k = _rskey(q)
            k in seen && continue
            push!(seen, k)
            push!(pts, q)
        end
    end
    return pts
end

# Every bond the structure forms, or `false` if it is not a structure at all -- copies overlapping,
# resting on a face nothing can use, or claiming a site twice.
function _rsbonds(cell::Polyform, gens)
    rules = bindingrules(cell)
    parts = cell.particles
    span, radius = _rsreach(cell)
    pts = _rspoints(gens, span)

    partner = fill(-1, nv(graphrep(cell)))
    for l in opensitelocs(cell)
        partner[first(bindingsite(cell, l).vertices)] = 0
    end
    out = RSClose{eltype(gens)}[]
    for t in pts
        _rsnear(parts, t, radius) || continue
        for part in parts
            ov, cts = _overlap_and_contacts(parts, translate(part, t), rules)
            ov && return false, out, pts
            for c in cts
                v1, v2 = first(c.vs1), first(c.vs2)
                partner[v1] == v2 && partner[v2] == v1 && continue
                (partner[v1] == 0 && partner[v2] == 0) || return false, out, pts
                partner[v1], partner[v2] = v2, v1
                push!(out, RSClose(v1, v2, t))
            end
        end
    end
    return true, out, pts
end

# Whether the structure repeats under a translation its own lattice does not contain, in which case
# a smaller cell describes it and this one should not be reported. Every such translation carries
# the first particle onto some particle of the cell, so the candidates are the differences from it.
function _rsreducible(cell::Polyform, pts)
    parts = cell.particles
    length(pts) > 1 || return false
    occupied = Dict{typeof(_rskey(first(pts))),Int}()
    for q in pts, (i, p) in enumerate(parts)
        occupied[_rskey(p.pose.x + q)] = i
    end
    for j in 2:length(parts)
        _samepose(parts[1], parts[j]) || continue
        t = parts[j].pose.x - parts[1].pose.x
        any(q -> q ≈ t, pts) && continue
        all(parts) do p
            i = get(occupied, _rskey(p.pose.x + t), 0)
            return i != 0 && _samepose(p, parts[i])
        end && return true
    end
    return false
end

### the two moves

mutable struct RSAux{BS,G}
    attachments::Vector{Tuple{BS,SpeciesSiteLoc,Int}}
    pairs::Vector{NTuple{2,Int}}      # first vertices of the two sites a close would identify
    seen::Set{G}                      # children already offered, since two pairs can identify the
    maxsize::Int                      # same structure and each would then be walked into
end
Base.copy(a::RSAux) = typeof(a)(copy(a.attachments), copy(a.pairs), copy(a.seen), a.maxsize)

# Offer a child once. A second route to the same structure is not a second structure.
function _rsonce!(aux::RSAux, u::RSTiling)
    u.key in aux.seen && return missing
    push!(aux.seen, copy(u.key))
    return u
end

# The state with the closes recomputed from the translations `gens`: whatever bonds the structure
# they generate actually forms. `nothing` if it forms none, or is no structure at all.
function _rsfrom(cell::Polyform, gens)
    ok, bonds, pts = _rsbonds(cell, gens)
    ok || return nothing
    _rsreducible(cell, pts) && return nothing
    isempty(bonds) && return _rekey!(RSTiling(cell, bonds, copy(graphrep(cell)), Int[]))

    # A structure has as many cells as there are ways to cut it, and the walk needs one of them:
    # the parent has to be a function of the structure, not of the cut it was reached through, or
    # the same tiling is walked into once per cut. `Tiling` cuts canonically, so read the cut back
    # off one.
    t = Tiling(cell, [_contactof(cell, c) for c in bonds])
    cut = unitcell(t)
    # one entry per bond, not per marker: a bond between two dart-encoded faces is pinned by
    # several pairs of vertices and wears a marker for each
    kept = empty(bonds)
    seen = Set{NTuple{2,ParticleSiteLoc}}()
    for (_, (u, v)) in _markers(t)
        _isperiodic(t, u, v; canonidxs=false) || continue
        pair = minmax(_vertex_to_particle_site(t, u; canonidxs=false),
                      _vertex_to_particle_site(t, v; canonidxs=false))
        pair in seen && continue
        push!(seen, pair)
        push!(kept, RSClose(u, v, _translation(t, u, v; canonidxs=false)))
    end
    return _rekey!(RSTiling(cut, kept, copy(graphrep(cut)), Int[]))
end

# The parent: undo an identification if any were made, otherwise remove a particle.
#
# Bonds cannot be undone one at a time. Some are consequences of the lattice rather than generators
# of it -- a cell bonded at `v` and at `w` is bonded at `v + w` as well -- and dropping one of those
# leaves the same structure and so the same state. Worse, no single bond need be essential: of
# those three, any two generate what all three do, so removing any one changes nothing while
# removing two changes everything.
#
# So drop them in canonical order and stop as soon as the state changes, which is graph inequality
# and asks nothing of the lattice. Dropping all of them certainly changes it, so this terminates.
function _rsls!(k::RSTiling, s::RSTiling)
    if isempty(s.closes)
        copy!(k, s)
        lower!(k.cell)
        empty!(k.closes)
        return _rekey!(k)
    end
    keep = trues(length(s.closes))
    for i in sortperm(s.markerat; rev=true)
        keep[i] = false
        gens = [c.t for (j, c) in enumerate(s.closes) if keep[j]]
        cand = _rsfrom(s.cell, gens)
        (isnothing(cand) || _samestate(cand, s)) && continue
        return copy!(k, cand)
    end
    return error("Internal error: a tiling state survived dropping all of its bonds. Please file an issue.")
end

# The children: attach a particle while nothing is identified yet, then identify a pair of sites.
# Attaching after identifying is not offered, since the parent rule undoes identifications first
# and would never lead back here.
function _rsadj!(u::RSTiling, v::RSTiling, j::Integer, aux::RSAux)
    rules = bindingrules(v.cell)
    if nparticles(v.cell) == 0
        j > nspecies(rules) && return nothing
        copy!(u.cell, Polyform(rules, j))
        empty!(u.closes)
        return _rekey!(u)
    end

    if j == 1
        empty!(aux.seen)
        empty!(aux.attachments)
        isempty(v.closes) && nparticles(v.cell) < aux.maxsize &&
            collect_attachments!(aux.attachments, v.cell)
        _rspairs!(aux.pairs, v)
    end

    if j <= length(aux.attachments)
        site, loc, t = aux.attachments[j]
        copy!(u.cell, v.cell)
        empty!(u.closes)
        out = raise!(u.cell, site, loc, t)
        (ismissing(out) || isnothing(out)) && return out
        return _rsonce!(aux, _rekey!(u))
    end

    i = j - length(aux.attachments)
    i > length(aux.pairs) && return nothing
    a, b = aux.pairs[i]
    sa = bindingsite(v.cell, _vertex_to_particle_site(v.cell, a; canonidxs=false))
    sb = bindingsite(v.cell, _vertex_to_particle_site(v.cell, b; canonidxs=false))
    gens = push!([c.t for c in v.closes], sa.pose.x - sb.pose.x)
    child = _rsfrom(v.cell, gens)
    isnothing(child) && return missing
    length(child.closes) > length(v.closes) || return missing
    return _rsonce!(aux, copy!(u, child))
end

# The pairs of sites a close could identify: still unbonded, able to bond, and facing each other,
# since only antiparallel sites meet under a translation.
function _rspairs!(out, s::RSTiling)
    cell = s.cell
    rules = bindingrules(cell)
    intmat = interactionmatrix(rules)
    taken = Set(v for c in s.closes for v in (c.a, c.b))
    free = [bindingsite(cell, l) for l in opensitelocs(cell)]
    filter!(b -> !(first(b.vertices) in taken), free)

    empty!(out)
    for i in eachindex(free), k in (i + 1):length(free)
        a, b = free[i], free[k]
        intmat[color(a), color(b)] || continue
        isaligned(a, b) || continue
        istouching(a, b) && continue
        push!(out, (first(a.vertices), first(b.vertices)))
    end
    return out
end

### the entry point

"""
    rstilings(rules::BindingRules; maxsize)

Enumerate the tilings of `rules` whose cell holds at most `maxsize` particles, by reverse search.

Every polyform of `rules` is reached, and offered every way of identifying a pair of its open
sites; what comes back are the states that identify at least one, are geometrically sound, and do
not repeat under a translation their own lattice misses.
"""
function rstilings(rules::BindingRules; maxsize::Integer)
    v₀ = RSTiling(Polyform(rules))
    BS = sitetype(rules)
    aux = RSAux(Tuple{BS,SpeciesSiteLoc,Int}[], NTuple{2,Int}[],
                Set{typeof(v₀.key)}(), Int(maxsize))
    rsys = RSSystem(_rsls!, _rsadj!, v₀; compare=_samestate, aux)

    out = RSTiling[]
    reversesearch(rsys) do s, _
        isempty(s.closes) || push!(out, copy(s))
        return ACCEPT
    end
    return out
end
