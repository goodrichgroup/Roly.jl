"""
One bond a cell makes with a translate of itself: the contact the two sites make, and the
translation that brings them together.

A search state is a polyform together with a set of these, which is the same thing a
[`Tiling`](@ref) is. Two moves reach every state: attaching a particle, which is `raise!`, and
bonding a pair of open sites across the cut. Bonding comes after attaching, because a particle
cannot be removed while a bond refers to its sites.
"""
struct PeriodicContact{V}
    contact::Contact
    t::V            # the translation carrying the second site onto the first
end

mutable struct TilingState{P<:Polyform,G<:AbstractNautyGraph,V}
    cell::P
    periodic::Vector{PeriodicContact{V}}   # the bonds made across the cut, in the order made
    graphrep::G                            # the cell's graph with a vertex per bond, canonized
    canonorder::Vector{Int}                # each periodic bond's first marker, in that graph's order
end

function TilingState(cell::Polyform)
    V = SVector{dimension(bindingrules(cell)),numtype(bindingrules(cell))}
    s = TilingState(cell, PeriodicContact{V}[], copy(graphrep(cell)), Int[])
    return _recanonize!(s)
end

graphrep(s::TilingState) = s.graphrep

function Base.copy(s::TilingState)
    return TilingState(copy(s.cell), copy(s.periodic), copy(s.graphrep), copy(s.canonorder))
end
function Base.copy!(dst::TilingState, src::TilingState)
    copy!(dst.cell, src.cell)
    copy!(dst.periodic, src.periodic)
    copy!(dst.graphrep, src.graphrep)
    copy!(dst.canonorder, src.canonorder)
    return dst
end

# Two states are the same when the structures they describe are, which is when their graphs are.
_samestate(a::TilingState, b::TilingState) = a.graphrep == b.graphrep

function Base.show(io::Core.IO, s::TilingState)
    return print(io, "TilingState[n=", nparticles(s.cell), ", periodic=", length(s.periodic), "]")
end

# Rebuild the graph from the cell and the periodic bonds: the cell's own graph, plus a vertex per
# bond joined to the two sites it joins, canonized. Records where each periodic bond's first marker
# landed, so that the parent can pick one of them in an order that does not depend on how the state
# was reached.
function _recanonize!(s::TilingState)
    rules = bindingrules(s.cell)
    g = NautyDiGraph(0)
    for part in s.cell.particles
        blockdiag!(g, graphrep(species(rules, speciesindex(part))))
    end

    # every bond is marked, the cell's own alongside the periodic ones. Marking only the latter
    # would record where the structure was cut, and two cuts of one tiling would then look like
    # two tilings
    marker = _markerlabel(rules)
    for e in exterior_edges(s.cell)
        _addmarker!(g, marker, toorig(s.cell, e.src), toorig(s.cell, e.dst))
    end
    firstmarker = Int[]
    for c in s.periodic
        push!(firstmarker, nv(g) + 1)
        for (v1, v2) in contact_pairing(c.contact)
            _addmarker!(g, marker, v1, v2)
        end
    end
    perm, _ = nauty(g; canonize=true)
    place = invperm(collect(Int, perm))
    s.graphrep = g
    resize!(s.canonorder, length(firstmarker))
    s.canonorder .= (place[v] for v in firstmarker)
    return s
end

"""
    PeriodicContact(cell::Polyform, u::Integer, v::Integer, t)

The bond the sites at graph vertices `u` and `v` of `cell` make when `t` carries the second onto
the first, which is the bond they would make side by side.

`u` and `v` are in `cell`'s original numbering, and name their sites by any one of their vertices.
"""
function PeriodicContact(cell::Polyform, u::Integer, v::Integer, t)
    a = bindingsite(cell, _vertex_to_particle_site(cell, u; canonidxs=false))
    b = bindingsite(cell, _vertex_to_particle_site(cell, v; canonidxs=false))
    moved = translate(b, t)
    return PeriodicContact(Contact(a.vertices, b.vertices, twist(a, moved), twistfreedom(a, moved)), t)
end

### the structure a state describes

# How wide the cell is, and how far a particle reaches past its own center.
function _cellreach(cell::Polyform)
    rules = bindingrules(cell)
    parts = cell.particles
    diameter = maximum(norm(p.pose.x - q.pose.x) for p in parts, q in parts)
    radius = maximum(bounding_radius(species(rules, speciesindex(p))) for p in parts)
    return diameter + 2radius, radius
end

_maytouch(parts, t, radius) =
    let d = norm(t)
        d > 0 || return false
        u = t / d
        lo, hi = extrema(dot(p.pose.x, u) for p in parts)
        return d <= (hi - lo) + 2radius + sqrt(eps(d))
    end

# A position as a dictionary key. Adding zero folds `-0.0` onto `0.0`, which are not `isequal` and
# so would otherwise land in different buckets, and a rounded position is exact enough to compare.
_poskey(v) = Tuple(round.(v; digits=7) .+ 0.0)

# The lattice points near enough to matter, walked outward from the origin along the translations
# the bonds name and their inverses. No basis is involved, which is the point: a set of bonds
# generates its lattice whether or not any subset of them is a basis of it.
function _latticepoints(gens, span)
    V = eltype(gens)
    pts = V[zero(V)]
    isempty(gens) && return pts
    seen = Set{typeof(_poskey(zero(V)))}([_poskey(zero(V))])
    i = 1
    while i <= length(pts)
        p = pts[i]
        i += 1
        for g in gens, sgn in (1, -1)
            q = p + sgn * g
            norm(q) <= 2span || continue
            k = _poskey(q)
            k in seen && continue
            push!(seen, k)
            push!(pts, q)
        end
    end
    return pts
end

# Every bond the structure forms, or `false` if it is not a structure at all -- copies overlapping,
# resting on a face nothing can use, or claiming a site twice.
function _closurebonds(cell::Polyform, gens)
    rules = bindingrules(cell)
    parts = cell.particles
    span, radius = _cellreach(cell)
    pts = _latticepoints(gens, span)

    partner = fill(-1, nv(graphrep(cell)))
    for l in opensitelocs(cell)
        partner[first(bindingsite(cell, l).vertices)] = 0
    end
    out = PeriodicContact{eltype(gens)}[]
    for t in pts
        _maytouch(parts, t, radius) || continue
        for part in parts
            ov, cts = _overlap_and_contacts(parts, translate(part, t), rules)
            ov && return false, out, pts
            for c in cts
                v1, v2 = first(c.vs1), first(c.vs2)
                partner[v1] == v2 && partner[v2] == v1 && continue
                (partner[v1] == 0 && partner[v2] == 0) || return false, out, pts
                partner[v1], partner[v2] = v2, v1
                push!(out, PeriodicContact(c, t))
            end
        end
    end
    return true, out, pts
end

# Whether the structure repeats under a translation its own lattice does not contain, in which case
# a smaller cell describes it and this one should not be reported. Every such translation carries
# the first particle onto some particle of the cell, so the candidates are the differences from it.
function _isreducible(cell::Polyform, pts)
    parts = cell.particles
    length(pts) > 1 || return false
    occupied = Dict{typeof(_poskey(first(pts))),Int}()
    for q in pts, (i, p) in enumerate(parts)
        occupied[_poskey(p.pose.x + q)] = i
    end
    for j in 2:length(parts)
        _samepose(parts[1], parts[j]) || continue
        t = parts[j].pose.x - parts[1].pose.x
        any(q -> q ≈ t, pts) && continue
        all(parts) do p
            i = get(occupied, _poskey(p.pose.x + t), 0)
            return i != 0 && _samepose(p, parts[i])
        end && return true
    end
    return false
end

### the two moves

mutable struct TilingAux{BS,G}
    attachments::Vector{Tuple{BS,SpeciesSiteLoc,Int}}
    pairs::Vector{NTuple{2,Int}}      # first vertices of the two sites a periodic bond would join
    seen::Set{G}                      # children already offered, since two pairs can identify the
    maxsize::Int                      # same structure and each would then be walked into
end
Base.copy(a::TilingAux) = typeof(a)(copy(a.attachments), copy(a.pairs), copy(a.seen), a.maxsize)

# Offer a child once. A second route to the same structure is not a second structure.
function _once!(aux::TilingAux, u::TilingState)
    u.graphrep in aux.seen && return missing
    push!(aux.seen, copy(u.graphrep))
    return u
end

# The state with its periodic bonds recomputed from the translations `gens`: whatever bonds the
# structure they generate actually forms. `nothing` if it forms none, or is no structure at all.
function _statefrom(cell::Polyform, gens)
    ok, bonds, pts = _closurebonds(cell, gens)
    ok || return nothing
    _isreducible(cell, pts) && return nothing
    isempty(bonds) && return _recanonize!(TilingState(cell, bonds, copy(graphrep(cell)), Int[]))

    # A structure has as many cells as there are ways to cut it, and the walk needs one of them:
    # the parent has to be a function of the structure, not of the cut it was reached through, or
    # the same tiling is walked into once per cut. `Tiling` cuts canonically, so read the cut back
    # off one.
    t = Tiling(cell, [c.contact for c in bonds])
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
        push!(kept, PeriodicContact(cut, u, v, _translation(t, u, v; canonidxs=false)))
    end
    return _recanonize!(TilingState(cut, kept, copy(graphrep(cut)), Int[]))
end

# The parent: undo a periodic bond if any were made, otherwise remove a particle.
#
# Bonds cannot be undone one at a time. Some are consequences of the lattice rather than generators
# of it -- a cell bonded at `v` and at `w` is bonded at `v + w` as well -- and dropping one of those
# leaves the same structure and so the same state. Worse, no single bond need be essential: of
# those three, any two generate what all three do, so removing any one changes nothing while
# removing two changes everything.
#
# So drop them in canonical order and stop as soon as the state changes, which is graph inequality
# and asks nothing of the lattice. Dropping all of them certainly changes it, so this terminates.
function ls!(k::TilingState, s::TilingState)
    if isempty(s.periodic)
        copy!(k, s)
        lower!(k.cell)
        empty!(k.periodic)
        return _recanonize!(k)
    end
    keep = trues(length(s.periodic))
    for i in sortperm(s.canonorder; rev=true)
        keep[i] = false
        gens = [c.t for (j, c) in enumerate(s.periodic) if keep[j]]
        cand = _statefrom(s.cell, gens)
        (isnothing(cand) || _samestate(cand, s)) && continue
        return copy!(k, cand)
    end
    return error("Internal error: a tiling state survived dropping all of its bonds. Please file an issue.")
end

# The children: attach a particle while nothing is bonded across the cut yet, then bond a pair of
# sites. Attaching afterwards is not offered, since the parent rule undoes those bonds first and
# would never lead back here.
function adj!(u::TilingState, v::TilingState, j::Integer, aux::TilingAux)
    rules = bindingrules(v.cell)
    if nparticles(v.cell) == 0
        j > nspecies(rules) && return nothing
        copy!(u.cell, Polyform(rules, j))
        empty!(u.periodic)
        return _recanonize!(u)
    end

    if j == 1
        empty!(aux.seen)
        empty!(aux.attachments)
        isempty(v.periodic) && nparticles(v.cell) < aux.maxsize &&
            collect_attachments!(aux.attachments, v.cell)
        _freepairs!(aux.pairs, v)
    end

    if j <= length(aux.attachments)
        site, loc, t = aux.attachments[j]
        copy!(u.cell, v.cell)
        empty!(u.periodic)
        out = raise!(u.cell, site, loc, t)
        (ismissing(out) || isnothing(out)) && return out
        return _once!(aux, _recanonize!(u))
    end

    i = j - length(aux.attachments)
    i > length(aux.pairs) && return nothing
    a, b = aux.pairs[i]
    sa = bindingsite(v.cell, _vertex_to_particle_site(v.cell, a; canonidxs=false))
    sb = bindingsite(v.cell, _vertex_to_particle_site(v.cell, b; canonidxs=false))
    gens = push!([c.t for c in v.periodic], sa.pose.x - sb.pose.x)
    child = _statefrom(v.cell, gens)
    isnothing(child) && return missing
    length(child.periodic) > length(v.periodic) || return missing
    return _once!(aux, copy!(u, child))
end

# The pairs of sites a periodic bond could join: still unbonded, able to bond, and facing each
# other, since only antiparallel sites meet under a translation.
function _freepairs!(out, s::TilingState)
    cell = s.cell
    rules = bindingrules(cell)
    intmat = interactionmatrix(rules)
    taken = Set(v for c in s.periodic for v in (first(c.contact.vs1), first(c.contact.vs2)))
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

# The tiling a state stands for.
_astiling(s::TilingState) = Tiling(s.cell, [c.contact for c in s.periodic])

"""
    tilingenum(f, rules::BindingRules; maxsize)

Enumerate the periodic closures of `rules` whose cell holds at most `maxsize` particles, streaming
each to `f` as a [`Tiling`](@ref) together with the number of particles its cell holds.

`f(t, n)` returns `ACCEPT`, `REJECT` to leave that tiling unextended, or `BREAK` to stop.

Every polyform of `rules` is reached and offered every way of identifying a pair of its open
sites, so a cell is any structure the rules admit rather than copies of one chosen block. What
comes back identifies at least one pair, is geometrically sound, and does not repeat under a
translation its own lattice misses.
"""
function tilingenum(f::F, rules::BindingRules; maxsize::Integer) where {F}
    v₀ = TilingState(Polyform(rules))
    BS = sitetype(rules)
    aux = TilingAux(Tuple{BS,SpeciesSiteLoc,Int}[], NTuple{2,Int}[], Set{typeof(v₀.graphrep)}(), Int(maxsize))
    rsys = RSSystem(ls!, adj!, v₀; compare=_samestate, aux)

    reversesearch(rsys) do s, _
        isempty(s.periodic) && return ACCEPT      # a polyform on the way to one, not a tiling
        return f(_astiling(s), nparticles(s.cell))
    end
    return nothing
end

"""
    tilings(rules::BindingRules; maxsize)

Return the periodic closures of `rules` whose cell holds at most `maxsize` particles.

See [`tilingenum`](@ref) to take them as they are found, and to stop early.
"""
function tilings(rules::BindingRules; maxsize::Integer)
    out = Tiling{dimension(rules),particletype(rules),typeof(rules),
                 typeof(graphrep(Polyform(rules))),SVector{dimension(rules),numtype(rules)}}[]
    tilingenum(rules; maxsize) do t, _
        push!(out, t)
        return ACCEPT
    end
    return out
end

"""
    _cellclosures(f, cell::Polyform)

Stream every periodic closure whose cell is `cell` itself to `f`, as a [`Tiling`](@ref).

Asking whether one structure closes is not asking what a rule set tiles into, and differs in two
ways. The cell is fixed, so there are no attaching moves and no walk over polyforms. And a closure
that repeats under a translation its lattice misses is kept rather than rejected: `cell` is a cell
of the structure it describes even where a smaller one would do, which is what the question means.
"""
function _cellclosures(f::F, cell::Polyform) where {F}
    start = TilingState(cell)
    seen = Set([copy(start.graphrep)])
    stack = [start]
    while !isempty(stack)
        s = pop!(stack)
        for (a, b) in _freepairs!(NTuple{2,Int}[], s)
            sa = bindingsite(cell, _vertex_to_particle_site(cell, a; canonidxs=false))
            sb = bindingsite(cell, _vertex_to_particle_site(cell, b; canonidxs=false))
            ok, bonds, _ = _closurebonds(cell, push!([c.t for c in s.periodic], sa.pose.x - sb.pose.x))
            (ok && length(bonds) > length(s.periodic)) || continue
            child = _recanonize!(TilingState(cell, bonds, copy(graphrep(cell)), Int[]))
            child.graphrep in seen && continue
            push!(seen, copy(child.graphrep))

            signal = f(_astiling(child))
            signal == BREAK && return BREAK
            signal == REJECT || push!(stack, child)
        end
    end
    return ACCEPT
end
