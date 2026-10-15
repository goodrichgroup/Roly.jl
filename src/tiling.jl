"""
    PeriodicContact{V}

A bond a cell makes with a translate of itself: the [`Contact`](@ref) the two sites make, and the
translation `t` that carries the second site onto the first.
"""
struct PeriodicContact{V}
    contact::Contact
    t::V
end

"""
    PeriodicContact(cell::Polyform, u::Integer, v::Integer, t)

Return the periodic contact between the sites at graph vertices `u` and `v` of `cell` when the
translation `t` carries the second onto the first.

`u` and `v` are in `cell`'s original numbering, and name their sites by any one of their vertices.
"""
function PeriodicContact(cell::Polyform, u::Integer, v::Integer, t)
    a = bindingsite(cell, _vertex_to_particle_site(cell, u; canonidxs=false))
    b = bindingsite(cell, _vertex_to_particle_site(cell, v; canonidxs=false))
    moved = translate(b, t)
    return PeriodicContact(Contact(a.vertices, b.vertices, twist(a, moved), twistfreedom(a, moved)), t)
end

"""
    Tiling{D}

A periodic structure: one cell, a [`Polyform`](@ref), together with the bonds it makes with its
own translates. Returned by [`tilings`](@ref).

Two tilings are equal when they describe the same structure, whichever cell either was found
through. The structure is identified by its graph, in which every bond, inside the cell and across
to a translate alike, is a vertex of its own joined to the two sites it bonds.
"""
mutable struct Tiling{D,P<:Polyform,G<:AbstractNautyGraph,V}
    cell::P
    periodic::Vector{PeriodicContact{V}}   # the bonds across the cut, in the order made
    graphrep::G                            # the cell's graph with a vertex per bond, canonized
    canon2orig::Vector{Int}
    sigma::Int
    canonmarkers::Vector{Int}              # where each periodic bond's leading marker sits
end

_vectortype(p::Polyform) = SVector{dimension(bindingrules(p)),numtype(bindingrules(p))}

Tiling(cell::Polyform) = Tiling(cell, PeriodicContact{_vectortype(cell)}[])
function Tiling(cell::Polyform{D}, periodic::Vector{PeriodicContact{V}}) where {D,V}
    t = Tiling{D,typeof(cell),typeof(graphrep(cell)),V}(cell, periodic, copy(graphrep(cell)), Int[], 1, Int[])
    return _recanonize!(t)
end

function Base.copy(t::Tiling)
    return typeof(t)(copy(t.cell), copy(t.periodic), copy(t.graphrep), copy(t.canon2orig), t.sigma, copy(t.canonmarkers))
end
function Base.copy!(dst::Tiling, src::Tiling)
    copy!(dst.cell, src.cell)
    copy!(dst.periodic, src.periodic)
    copy!(dst.graphrep, src.graphrep)
    copy!(dst.canon2orig, src.canon2orig)
    dst.sigma = src.sigma
    copy!(dst.canonmarkers, src.canonmarkers)
    return dst
end

Base.:(==)(a::Tiling, b::Tiling) = bindingrules(a) === bindingrules(b) && a.graphrep == b.graphrep
Base.hash(t::Tiling, h::UInt) = hash(t.graphrep, h)

# Two tilings are the same when the structures they describe are, which is when their graphs are.
_samestate(a::Tiling, b::Tiling) = a.graphrep == b.graphrep

function Base.show(io::Core.IO, t::Tiling{D}) where {D}
    return print(io, "Tiling{$D}[n=$(nparticles(t)), d=$(length(translations(t)))]")
end

### the graph

"""
    _markerlabel(rules::BindingRules)

Return a graph vertex label no species of `rules` uses, worn by the vertices that stand in for
bonds.
"""
_markerlabel(rules::BindingRules) = maximum(maximum(labels(graphrep(species(rules, i)))) for i in 1:nspecies(rules)) + 1

# Record the bond between graph vertices `u` and `v` as a vertex of its own joined to both, rather
# than as the edge `u -- v`.
#
# A particle of a tiling bonds to its own translate, so a bond can join two sites of one particle,
# and the species may already carry an edge between exactly those two vertices, as a cycle
# encoding does for the cube's opposite faces. Writing the bond as an edge would then write
# nothing, and the tiling would be indistinguishable from one without the bond.
function _addmarker!(g, marker::Integer, u::Integer, v::Integer)
    add_vertices!(g, 1; vertex_labels=[marker])
    m = nv(g)
    add_edge!(g, u, m)
    add_edge!(g, m, u)
    add_edge!(g, v, m)
    add_edge!(g, m, v)
    return g
end

# Rebuild the graph from the cell and the periodic bonds: the cell's own graph, plus a vertex per
# bond joined to the two sites it joins, canonized. Every bond is marked, the cell's own alongside
# the periodic ones, so that the graph records that a bond exists and not where the structure was
# cut: two cuts of one tiling give one graph.
#
# Also record where each periodic bond's leading marker landed. One bond wears several markers when
# its sites are dart-encoded faces, so the leading one stands for it. These positions are what the
# parent orders the bonds by, which ties its choice to the tiling's own canonical graph rather than
# to the order the bonds were listed in.
function _recanonize!(t::Tiling)
    rules = bindingrules(t.cell)
    g = NautyDiGraph(0)
    for part in t.cell.particles
        blockdiag!(g, graphrep(species(rules, speciesindex(part))))
    end
    marker = _markerlabel(rules)
    for e in exterior_edges(t.cell)
        _addmarker!(g, marker, toorig(t.cell, e.src), toorig(t.cell, e.dst))
    end
    leading = Int[]
    for c in t.periodic
        push!(leading, nv(g) + 1)
        for (v1, v2) in contact_pairing(c.contact)
            _addmarker!(g, marker, v1, v2)
        end
    end
    perm, autg = nauty(g; canonize=true)
    t.graphrep = g
    t.canon2orig = collect(Int, perm)
    t.sigma = round(Int, autg.n)
    place = invperm(t.canon2orig)
    t.canonmarkers = [place[v] for v in leading]
    return t
end

# The marker vertices of `t`, in canonical order, each with the two site vertices it joins, in the
# cell's original numbering.
function _markers(t::Tiling)
    ml = _markerlabel(bindingrules(t))
    ls = labels(t.graphrep)
    return ((m, _markerends(t, m)) for m in vertices(t.graphrep) if ls[m] == ml)
end

function _markerends(t::Tiling, m::Integer)
    neighs = NautyGraphs.adjrow(t.graphrep, m)
    ends = (0, 0)
    for w in eachindex(neighs)
        neighs[w] || continue
        ends = first(ends) == 0 ? (t.canon2orig[w], 0) : (first(ends), t.canon2orig[w])
    end
    return ends
end

# Whether the bond between two of the cell's graph vertices crosses to a translate, which is to say
# the two sites it joins do not meet inside the cell. Position alone decides it: two sites that
# mate are antiparallel, so their poses never agree.
function _isperiodic(t::Tiling, u::Integer, v::Integer)
    a = bindingsite(t.cell, _vertex_to_particle_site(t.cell, u; canonidxs=false))
    b = bindingsite(t.cell, _vertex_to_particle_site(t.cell, v; canonidxs=false))
    return !isapprox(a.pose.x, b.pose.x; atol=a.touching_tolerance + b.touching_tolerance, rtol=0)
end

# The translation carrying the site at `v` onto the site at `u`.
function _translation(t::Tiling, u::Integer, v::Integer)
    a = bindingsite(t.cell, _vertex_to_particle_site(t.cell, u; canonidxs=false))
    b = bindingsite(t.cell, _vertex_to_particle_site(t.cell, v; canonidxs=false))
    return a.pose.x - b.pose.x
end

### the canonical cell

# The particle holding the first vertex in canonical order, which is where the cut starts.
function _rootparticle(t::Tiling)
    for v in vertices(t.graphrep)
        loc = _vertex_to_particle_site(t.cell, t.canon2orig[v]; canonidxs=false)
        isnothing(loc) || return loc.particle
    end
    return 1
end

# The first bond in canonical order that reaches an unplaced particle from a placed one, as that
# particle and the translation laying it against its partner. `nothing` once all are placed.
function _reachbond(t::Tiling, placed)
    for (_, (u, v)) in _markers(t)
        a = _vertex_to_particle_site(t.cell, u; canonidxs=false)
        b = _vertex_to_particle_site(t.cell, v; canonidxs=false)
        (isnothing(a) || isnothing(b)) && continue
        placed[a.particle] == placed[b.particle] && continue
        sa, sb = bindingsite(t.cell, a), bindingsite(t.cell, b)
        placed[a.particle] && return b.particle, sa.pose.x - sb.pose.x
        return a.particle, sb.pose.x - sa.pose.x
    end
    return nothing
end

# Re-choose which translate of each particle the cell holds, so that where the structure is cut
# follows from the canonical graph rather than from the cell the search happened to grow. Walking
# the canonical spanning tree and pulling each particle onto the bond that reaches it leaves every
# tree bond meeting inside the cell. Two equal tilings that arrived cut differently come out cut
# the same.
function _canonicalcut!(t::Tiling)
    parts = t.cell.particles
    placed = falses(length(parts))
    placed[_rootparticle(t)] = true
    while true
        step = _reachbond(t, placed)
        isnothing(step) && break
        j, delta = step
        parts[j] = translate(parts[j], delta)
        placed[j] = true
    end
    return t
end

# Place the cell in the frame of its root particle, so two equal tilings carry equal poses. A pose
# is a full frame on its own, so one particle is enough in any dimension.
function _canonicalframe!(t::Tiling)
    parts = t.cell.particles
    q = inv(parts[_rootparticle(t)].pose)
    for i in eachindex(parts)
        parts[i] = q * parts[i]
    end
    return t
end

# The cell `t`'s particles make as they now stand: a polyform of them, bonded wherever a bond's two
# sites meet inside the cell.
function _cutcell(t::Tiling{D}) where {D}
    rules = bindingrules(t)
    g = NautyDiGraph(0)
    for part in t.cell.particles
        blockdiag!(g, graphrep(species(rules, speciesindex(part))))
    end
    for (_, (u, v)) in _markers(t)
        _isperiodic(t, u, v) && continue
        add_edge!(g, u, v)
        add_edge!(g, v, u)
    end
    perm, autg = nauty(g; canonize=true)
    cvs = collect(Int, perm)
    return typeof(t.cell)(g, round(Int, autg.n), cvs, invperm(cvs), copy(t.cell.particles), rules)
end

# The tiling of `cell` closed by the periodic bonds `bonds`, cut and framed canonically. The cut
# re-chooses which bonds cross to a translate, so both the cell and the periodic bonds are read
# back off the moved particles.
function _canonical(cell::Polyform, bonds)
    t = Tiling(copy(cell), bonds)
    _canonicalcut!(t)
    _canonicalframe!(t)
    cut = _cutcell(t)
    # one entry per bond, not per marker: a bond between two dart-encoded faces is pinned by
    # several pairs of vertices and wears a marker for each
    periodic = empty(bonds)
    seen = Set{NTuple{2,ParticleSiteLoc}}()
    for (_, (u, v)) in _markers(t)
        _isperiodic(t, u, v) || continue
        pair = minmax(_vertex_to_particle_site(cut, u; canonidxs=false), _vertex_to_particle_site(cut, v; canonidxs=false))
        pair in seen && continue
        push!(seen, pair)
        push!(periodic, PeriodicContact(cut, u, v, _translation(t, u, v)))
    end
    return Tiling(cut, periodic)
end

### what a tiling is made of

"""
    unitcell(t::Tiling)

Return one cell of `t` as a [`Polyform`](@ref): its particles, keeping only the bonds that do not
cross to a translate.

Which bonds those are depends on where the structure is cut, and `t` is cut canonically, so equal
tilings give equal cells.
"""
unitcell(t::Tiling) = t.cell

"""
    tilingorder(t::Tiling)

Return how many particles one cell of `t` holds.
"""
tilingorder(t::Tiling) = nparticles(t)

"""
    translations(t::Tiling)

Return the translations the periodic bonds of `t` name, each once.

They generate the lattice `t` repeats under, but need not be a basis of it. Fewer than the
dimension is a partial closure, a column in the plane or a column or sheet in space. More arise
when the cell bonds to one translate through several pairs of sites, or along a sum of its other
translations.
"""
function translations(t::Tiling{D,P,G,V}) where {D,P,G,V}
    out = V[]
    for c in t.periodic
        v = _orient(c.t)
        any(w -> w ≈ v, out) || push!(out, v)
    end
    return out
end

_tol(::AbstractArray{F}) where {F} = sqrt(eps(F))

# A vector and its negative generate the same translations, so pick the one whose first
# significant component is positive.
function _orient(v::AbstractVector)
    i = findfirst(x -> abs(x) > _tol(v), v)
    return isnothing(i) || v[i] > 0 ? v : -v
end

@inline bindingrules(t::Tiling) = bindingrules(t.cell)
@inline nparticles(t::Tiling) = nparticles(t.cell)
@inline nsites(t::Tiling) = nsites(t.cell)
@inline dimension(::Tiling{D}) where {D} = D
@inline symmetrynumber(t::Tiling) = t.sigma
@inline graphrep(t::Tiling) = t.graphrep
@inline bindingsite(t::Tiling, loc::ParticleSiteLoc) = bindingsite(t.cell, loc)

# The two sites a periodic contact bonds, as addresses into the cell.
function _bondedsites(t::Tiling, c::PeriodicContact)
    a = _vertex_to_particle_site(t.cell, first(c.contact.vs1); canonidxs=false)
    b = _vertex_to_particle_site(t.cell, first(c.contact.vs2); canonidxs=false)
    return a, b
end

"""
    bonds(t::Tiling)

Return the bonds of one cell of `t`, each one once, as [`ParticleSiteLoc`](@ref) pairs. The bonds
a translate closes are listed alongside the ones inside the cell.
"""
function bonds(t::Tiling)
    out = Pair{ParticleSiteLoc,ParticleSiteLoc}[]
    for b in bonds(t.cell)
        push!(out, first(b) => last(b))
    end
    for c in t.periodic
        a, b = _bondedsites(t, c)
        push!(out, a => b)
    end
    return out
end

"""
    nbonds(t::Tiling)

Return the number of bonds of one cell of `t`, counting those a translate closes.
"""
nbonds(t::Tiling) = nbonds(t.cell) + length(t.periodic)

"""
    bondindex(t::Tiling, a::ParticleSiteLoc, b::ParticleSiteLoc)

Return the index into `bonded_colors(bindingrules(t))` for a bond between the sites `a` and `b`
of `t`, or `nothing` if their colors form no valid bond type.
"""
bondindex(t::Tiling, a::ParticleSiteLoc, b::ParticleSiteLoc) = bondindex(t.cell, a, b)

"""
    bondtypes(t::Tiling)

Return the bond type of every bond of one cell of `t`, each one once, indexing
`bonded_colors(bindingrules(t))`. The bonds a translate closes are counted alongside the ones
inside the cell.
"""
function bondtypes(t::Tiling)
    return map(bonds(t)) do (a, b)
        i = bondindex(t, a, b)
        isnothing(i) && error("Internal error: a bond has no bond type. Please file an issue.")
        return i
    end
end

# The sites of the cell that periodic bonds use up.
function _periodicsites(t::Tiling)
    out = Set{ParticleSiteLoc}()
    for c in t.periodic
        push!(out, _bondedsites(t, c)...)
    end
    return out
end

"""
    opensitelocs(t::Tiling)

Return the addresses of the sites of one cell of `t` that no bond uses, inside the cell or across
to a translate.
"""
function opensitelocs(t::Tiling)
    used = _periodicsites(t)
    return filter(l -> !(l in used), opensitelocs(t.cell))
end

"""
    exposedsitelocs(t::Tiling)

Return the addresses of the sites of one cell of `t` that no bond uses and that some bond could.
"""
function exposedsitelocs(t::Tiling)
    used = _periodicsites(t)
    return filter(l -> !(l in used), exposedsitelocs(t.cell))
end

"""
    iscomplete(t::Tiling)

Return `true` if `t` leaves no site open: every site of every particle is closed, by a neighbor in
the cell or by a translate. A tiling that is not complete closes along fewer directions than the
space has.
"""
iscomplete(t::Tiling) = isempty(opensitelocs(t))

### the structure a cell and its translations describe

# How wide the cell is, and how far a particle reaches past its own center.
function _cellreach(cell::Polyform)
    rules = bindingrules(cell)
    parts = cell.particles
    diameter = maximum(norm(p.pose.x - q.pose.x) for p in parts, q in parts)
    radius = maximum(bounding_radius(species(rules, speciesindex(p))) for p in parts)
    return diameter + 2radius, radius
end

# Whether a copy of the cell displaced by `t` can touch the cell: closer along `t` than the cell's
# own extent that way plus two particle radii.
function _maytouch(parts, t, radius)
    d = norm(t)
    d > 0 || return false
    u = t / d
    lo, hi = extrema(dot(p.pose.x, u) for p in parts)
    return d <= (hi - lo) + 2radius + sqrt(eps(d))
end

# A position as a dictionary key. Adding zero folds `-0.0` onto `0.0`, which are not `isequal` and
# so would otherwise land in different buckets, and a rounded position is exact enough to compare.
_poskey(v) = Tuple(round.(v; digits=7) .+ 0.0)

# The lattice points near enough to matter, walked outward from the origin along the translations
# the bonds name and their inverses. No basis is involved: a set of bonds generates its lattice
# whether or not any subset of them is a basis of it.
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

# Every bond the structure forms, or `false` if it is not a structure at all: copies overlapping,
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

# A translation carries a particle onto another only if the two wear the same species in the same
# orientation.
_samepose(a::Particle, b::Particle) = speciesindex(a) == speciesindex(b) && a.pose.psi ≈ b.pose.psi

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

# The tiling of `cell` closed by whatever bonds the structure generated by the translations `gens`
# forms, cut canonically. `nothing` if that structure is invalid, or a smaller cell describes it.
function _tilingfrom(cell::Polyform, gens)
    ok, bonds, pts = _closurebonds(cell, gens)
    ok || return nothing
    _isreducible(cell, pts) && return nothing
    isempty(bonds) && return Tiling(copy(cell), bonds)
    return _canonical(cell, bonds)
end

### the two moves of the search

mutable struct TilingAux{BS,G}
    attachments::Vector{Tuple{BS,SpeciesSiteLoc,Int}}
    pairs::Vector{NTuple{2,Int}}      # first vertices of the two sites a periodic bond would join
    seen::Set{G}                      # children already offered, since two pairs can close the
    maxsize::Int                      # same structure and each would then be walked into
end
Base.copy(a::TilingAux) = typeof(a)(copy(a.attachments), copy(a.pairs), copy(a.seen), a.maxsize)

# Offer a child once. A second route to the same structure is not a second structure.
function _once!(aux::TilingAux, u::Tiling)
    u.graphrep in aux.seen && return missing
    push!(aux.seen, copy(u.graphrep))
    return u
end

# The parent: undo a periodic bond if any were made, otherwise remove a particle.
#
# Bonds cannot be undone one at a time. Some are consequences of the lattice rather than generators
# of it, as a cell bonded at `v` and at `w` is bonded at `v + w` as well, and dropping one of those
# leaves the same structure. Nor need any single bond be essential: of those three, any two
# generate what all three do. So drop them in the order their markers take in the canonical graph
# and stop as soon as the structure changes, which is graph inequality. Dropping all of them
# certainly changes it, so this terminates.
function ls!(k::Tiling, s::Tiling)
    if isempty(s.periodic)
        copy!(k, s)
        lower!(k.cell)
        return _recanonize!(k)
    end
    keep = trues(length(s.periodic))
    for i in sortperm(s.canonmarkers; rev=true)
        keep[i] = false
        gens = [c.t for (j, c) in enumerate(s.periodic) if keep[j]]
        cand = _tilingfrom(s.cell, gens)
        (isnothing(cand) || _samestate(cand, s)) && continue
        return copy!(k, cand)
    end
    return error("Internal error: a tiling survived dropping all of its periodic bonds. Please file an issue.")
end

# The children: attach a particle while nothing is bonded across the cut yet, then bond a pair of
# sites. Attaching afterwards is not offered, since the parent rule undoes those bonds first and
# would never lead back here.
function adj!(u::Tiling, v::Tiling, j::Integer, aux::TilingAux)
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
    child = _tilingfrom(v.cell, gens)
    isnothing(child) && return missing
    length(child.periodic) > length(v.periodic) || return missing
    return _once!(aux, copy!(u, child))
end

# The pairs of sites a periodic bond could join: still unbonded, able to bond, and facing each
# other, since only antiparallel sites meet under a translation.
function _freepairs!(out, t::Tiling)
    cell = t.cell
    rules = bindingrules(cell)
    intmat = interactionmatrix(rules)
    taken = Set(v for c in t.periodic for v in (first(c.contact.vs1), first(c.contact.vs2)))
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

### the entry points

"""
    tilingenum(f, rules::BindingRules; maxsize)

Enumerate the periodic structures of `rules` whose cell holds at most `maxsize` particles,
calling `f(t, n)` on each as a [`Tiling`](@ref) `t` with `n` particles in its cell.

`f` returns `ACCEPT` to continue, `REJECT` to leave that tiling unextended, or `BREAK` to stop.

Every polyform of `rules` is reached and offered every way of bonding a pair of its open sites
across to a translate, so a cell is any structure the rules admit. What comes back closes at least
one such bond, is geometrically sound, and does not repeat under a translation its own lattice
misses. Each structure is reported once, whichever cells it could be cut into.
"""
function tilingenum(f::F, rules::BindingRules; maxsize::Integer) where {F}
    v₀ = Tiling(Polyform(rules))
    BS = sitetype(rules)
    aux = TilingAux(Tuple{BS,SpeciesSiteLoc,Int}[], NTuple{2,Int}[], Set{typeof(v₀.graphrep)}(), Int(maxsize))
    rsys = RSSystem(ls!, adj!, v₀; compare=_samestate, aux)

    reversesearch(rsys) do t, _
        isempty(t.periodic) && return ACCEPT      # a polyform on the way to a tiling, not one
        return f(copy(t), nparticles(t))
    end
    return nothing
end

"""
    tilings(rules::BindingRules; maxsize)

Return the periodic structures of `rules` whose cell holds at most `maxsize` particles, as
[`Tiling`](@ref)s.

See [`tilingenum`](@ref) to take them as they are found, and to stop early.
"""
function tilings(rules::BindingRules; maxsize::Integer)
    out = typeof(Tiling(Polyform(rules)))[]
    tilingenum(rules; maxsize) do t, _
        push!(out, t)
        return ACCEPT
    end
    return out
end
