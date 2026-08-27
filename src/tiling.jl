"""
    Tiling{V}

One periodic closure of a polyform: the lattice its copies form, and what one cell of that
lattice holds. Returned by [`tilings`](@ref).

  - `cell`: the unit cell, as a [`Polyform`](@ref) of the rules the tiled polyform was built under
  - `vectors`: the generating translations, at most `dimension` of them
  - `bondtypes`: the bond type of every bond one cell contributes, those inside the cell included
  - `complete`: whether the closure leaves no open site
  - `order`: how many copies of the tiled polyform the cell holds

A cell's composition is `order` copies of the polyform plus `bondtypes`. Reach the parts through
[`unitcell`](@ref), [`latticevectors`](@ref), [`bondtypes`](@ref), [`iscomplete`](@ref) and
[`tilingorder`](@ref) rather than the fields.
"""
struct Tiling{PF,V}
    cell::PF
    vectors::Vector{V}
    bondtypes::Vector{Int}
    complete::Bool
    order::Int
end

function Base.show(io::Core.IO, t::Tiling)
    return print(io, "Tiling[d=$(length(t.vectors)), order=$(t.order)", t.complete ? ", complete]" : "]")
end

function Base.:(==)(a::Tiling, b::Tiling)
    (iscomplete(a) == iscomplete(b) && tilingorder(a) == tilingorder(b)) || return false
    unitcell(a) == unitcell(b) || return false
    sort(bondtypes(a)) == sort(bondtypes(b)) || return false
    va, vb = latticevectors(a), latticevectors(b)
    length(va) == length(vb) || return false
    return all(v -> any(w -> _orient(v) ≈ _orient(w), vb), va)
end

"""
    unitcell(t::Tiling)

The cell `t` repeats, as a [`Polyform`](@ref): `order(t)` copies of the polyform that was tiled,
bonded to each other as they are inside the cell.
"""
unitcell(t::Tiling) = t.cell

"""
    latticevectors(t::Tiling)

The translations that generate `t`, at most one per dimension. Fewer means a partial closure: a
column in the plane, a column or a sheet in space.
"""
latticevectors(t::Tiling) = t.vectors

"""
    bondtypes(t::Tiling)

The bond type of every bond one cell of `t` contributes, indexing [`bonded_colors`](@ref). The
bonds inside the cell are counted alongside the ones its translates close.
"""
bondtypes(t::Tiling) = t.bondtypes

"""
    iscomplete(t::Tiling)

Whether `t` leaves no open site: every site of the cell is closed by one of its translates. A
tiling that is not complete closes along fewer directions than the space has, i.e. a 1d line in the 2d plane.
"""
iscomplete(t::Tiling) = t.complete

"""
    tilingorder(t::Tiling)

How many copies of the tiled polyform one cell of `t` holds.
"""
tilingorder(t::Tiling) = t.order

"""
    dimension(t::Tiling)

The dimension of the space `t` lives in, not the number of directions it closes -- that is
`length(latticevectors(t))`.
"""
dimension(t::Tiling) = dimension(unitcell(t))

"""
    tilings(poly::Polyform; nreps=2, maxorder=1)

Enumerate the periodic closures of `poly`: choices of up to `dimension` translation vectors
under which copies of a cell form valid bonds with each other and never overlap, each vector
contributing at least one bond.

  - `nreps`: neighbor shells placed and checked per vector.
  - `maxorder`: let a cell hold up to `maxorder` copies of `poly`, grown as a meta-polyform and
    then checked for translation tilings. Above 1 this reaches tilings in which `poly` appears in
    rotated configuraitons.
  - returns a vector of [`Tiling`](@ref)

Candidate vectors connect interacting, aligned pairs of open sites, so a closure that needs a
turn is found only through the cell it makes: growing cells with [`MetaParticleSpecies`](@ref)
puts the turn inside the cell, leaving a lattice that is a pure translation again.
"""
function tilings(poly::Polyform; nreps::Integer=2, maxorder::Integer=1)
    rules = bindingrules(poly)
    out = Tiling{typeof(poly),SVector{dimension(rules),numtype(rules)}}[]
    cells, metarules = _tilecells(poly, maxorder)
    for cell in cells
        # the cell is carried through the search as a meta-polyform and read back as a polyform
        # of `rules` once, here, rather than once per candidate
        _addtilings!(out, cell, metarules, rules, recast(cell.meta, rules), nreps)
    end
    return unique!(out)
end

function Base.hash(t::Tiling, h::UInt)
    h = hash(graphrep(unitcell(t)), h)
    h = hash(sort(bondtypes(t)), h)
    h = hash(iscomplete(t), h)
    h = hash(tilingorder(t), h)
    return hash(length(latticevectors(t)), h)
end

# spatial units are O(1), so use sqrt(eps) for absolute tolerance
_tol(::AbstractVector{F}) where {F} = sqrt(eps(F))

# A vector and its negative generate the same translations, so pick the one whose first
# significant component is positive.
function _orient(v::AbstractVector)
    i = findfirst(x -> abs(x) > _tol(v), v)
    return isnothing(i) || v[i] > 0 ? v : -v
end

"""
    TileCell{PF,S}

One candidate cell of a lattice, kept as the meta-polyform the search grows it from.

  - `meta`: the cell, `order` copies of the tiled polyform wrapped as one meta-polyform
  - `sites`: its meta sites, in [`bindingsites`](@ref) order
  - `spent`: which of those a translate can no longer use -- bound inside the cell, or inert
  - `basecolors`: the color each meta site carries in the *underlying* rules, for naming bond types
  - `order`: how many copies of the polyform the cell holds
"""
struct TileCell{PF,S}
    meta::PF
    sites::Vector{S}
    spent::BitVector
    basecolors::Vector{Int}
    order::Int
end

# Everything the shell search threads around: the cell being translated, what it is searched
# against, and the state of the current branch of the search.
struct _ShellSearch{PF,P,S,R,V}
    cell::TileCell{PF,S}
    rules::R                           # the *meta* rules: the cell is translated as meta-particles
    vectors::Vector{V}
    ofvertex::Vector{Int}              # first vertex of a site -> its index in `cell.sites`, else 0
    nreps::Int
    consumed::BitVector                # which of `cell.sites` a bond has already closed
    contacts::Vector{NTuple{2,Int}}    # bonds the placed shells have formed, as site indices
    placed::Vector{Vector{P}}          # one shell of copies per chosen vector
    chosen::Vector{Int}                # indices into `vectors`
end

# A contact comes back as the vertex range its site occupies, and a translated copy keeps the
# ranges of the cell it was translated from, so the first vertex identifies the site. An array
# rather than a `Dict`: the vertices are a contiguous range, so indexing is the whole lookup.
function _siteofvertex(cell::TileCell)
    ofvertex = zeros(Int, maximum(s -> last(s.vertices), cell.sites; init=0))
    for (i, s) in enumerate(cell.sites)
        ofvertex[first(s.vertices)] = i
    end
    return ofvertex
end
_lookup(s::_ShellSearch, vs) = first(vs) <= length(s.ofvertex) ? s.ofvertex[first(vs)] : 0

function _addtilings!(out, cell::TileCell, metarules, rules, cellpoly, nreps::Integer)
    free = [s for (i, s) in enumerate(cell.sites) if !cell.spent[i]]
    isempty(free) && return out
    vectors = _candidatelatticevectors(free, metarules)
    isempty(vectors) && return out

    search = _ShellSearch(
        cell,
        metarules,
        vectors,
        _siteofvertex(cell),
        Int(nreps),
        copy(cell.spent),
        NTuple{2,Int}[],
        Vector{eltype(cell.meta.particles)}[],
        Int[],
    )
    # the bonds inside the cell count towards what it spends, alongside the ones its translates
    # close. Both are named by the colors the sites carry in the underlying rules, not the meta
    # ones, since a bond type indexes that system's bonds
    inside = [(siteindex(cell.meta, a), siteindex(cell.meta, b)) for (a, b) in bonds(cell.meta)]
    spent = [_bondtype(rules, cell.basecolors, c) for c in inside]
    record!() = push!(
        out,
        Tiling(
            cellpoly,
            vectors[search.chosen],
            vcat(spent, [_bondtype(rules, cell.basecolors, c) for c in search.contacts]),
            all(search.consumed),
            cell.order,
        ),
    )
    _tilings!(record!, search, dimension(rules), 1)
    return out
end

# The cells to try: every meta-polyform of up to `maxorder` copies of `poly` that the assembly
# machinery can build. `maxorder == 1` needs no special case, since the meta-monomer is `poly`
# itself.
#
# `exposeinert=true` is not decoration. The search rejects a placement in which two cells merely
# *touch* at a pair of sites nothing can bond, and it can only see the sites the meta-species
# exposes. Left to its default the species would expose the bondable ones alone, and two cells
# resting against each other on inert faces would look like empty space -- `overlap` does not
# catch it either, since touching is not interpenetrating. Exposed, they are meta sites like any
# other, the lift leaves them inert, and the rejection fires.
function _tilecells(poly::Polyform, maxorder::Integer)
    mp = MetaParticleSpecies(poly; exposeinert=true)
    metarules = BindingRules(mp)
    S = sitetype(metarules)
    PF = typeof(Polyform(metarules))
    cells = TileCell{PF,S}[]
    isempty(opensites(poly)) && return cells, metarules

    ucolors = _underlyingcolors(mp)
    polyenum(metarules; maxsize=maxorder) do meta, _
        sites = collect(bindingsites(meta))
        # a translate can use a site that is unbound inside the cell and not inert; everything
        # else is spent before the search starts, and a `complete` closure is one that uses up
        # what is left
        spent = trues(length(sites))
        for l in opensites(meta)
            spent[siteindex(meta, l)] = false
        end
        basecolors = [ucolors[k] for _ in meta.particles for k in 1:nsites(mp)]
        push!(cells, TileCell(copy(meta), sites, spent, basecolors, nparticles(meta)))
        return ACCEPT
    end
    return cells, metarules
end

"""
    isunitcell(poly::Polyform; kwargs...)

Whether `poly` tiles space by translations with every open site closed — see [`tilings`](@ref).
"""
isunitcell(poly::Polyform; kwargs...) = any(t.complete for t in tilings(poly; kwargs...))

"""
    tilelatticevectors(poly::Polyform; kwargs...)

The lattice vectors of the first complete tiling, or `nothing` — see [`tilings`](@ref).
"""
function tilelatticevectors(poly::Polyform; kwargs...)
    i = findfirst(t -> t.complete, tilings(poly; kwargs...))
    return i === nothing ? nothing : tilings(poly; kwargs...)[i].vectors
end

"""
    cantile(rules::BindingRules; maxtilesize, kwargs...)

Search the structures of `rules` up to `maxtilesize` particles for a unit cell; returns the
first one found, or `nothing`.
"""
function cantile(rules::BindingRules; maxtilesize, kwargs...)
    hit = Ref{Any}(nothing)
    f(s, _) = isunitcell(s; kwargs...) ? (hit[]=copy(s); BREAK) : ACCEPT
    polyenum(f, rules; maxsize=maxtilesize)
    return hit[]
end

"""
    canchain(rules::BindingRules; maxlength=chainstatebound(rules) + 1, kwargs...)

Grow linear chains of `rules`, attaching only at their two ends, and return the first one that
closes periodically — or `nothing` if none does within `maxlength` particles.

  - `maxlength`: longest chain to try. The default is long enough that no periodic chain can hide
    below it, see [`chainstatebound`](@ref); lower it only to trade certainty for speed
  - other keyword arguments go to [`tilings`](@ref)

A chain that closes proves the rules admit arbitrarily large structures: truncating its periodic
continuation gives a valid polyform of every length. This is far cheaper than [`cantile`](@ref),
which sweeps every structure, because chains branch only at their ends, and a system whose
structures all close runs out of chains on its own well before the bound.

Restricting the search to chains costs nothing in reach. Any connected part of a valid structure
is valid here, and an infinite structure's bond graph is connected and locally finite, so it
contains an infinite path: branching growth always implies chain growth.

**In 2D** `nothing` at the default `maxlength` therefore means no infinite structure exists. A
chain long enough to repeat one of its finitely many states repeats the rigid motion between the
two occurrences forever, and in the plane that motion is a translation or a rotation — a rotation
keeps every particle on a circle, so its iterates must eventually collide. Unbounded growth
leaves only the translation, which is a periodic chain.

**In 3D** it does not: the motion can be a screw, whose translation along the axis never
collides. A screw by `2πp/q` closes into a translation over `q` copies and is found with
`maxorder ≥ q`, but an irrational one never closes at all, and no periodicity test can see it.
See [`isunbounded`](@ref).
"""
function canchain(rules::BindingRules;
    maxlength::Integer=chainstatebound(rules) + 1, kwargs...)
    frontier = [Polyform(rules, i) for i in 1:nspecies(rules)]
    seen = Set(hash(graphrep(p)) for p in frontier)
    while !isempty(frontier)
        poly = popfirst!(frontier)
        isempty(tilings(poly; kwargs...)) || return poly
        nparticles(poly) < maxlength || continue
        for child in _extendends(poly)
            hash(graphrep(child)) in seen && continue
            push!(seen, hash(graphrep(child)))
            push!(frontier, child)
        end
    end
    return nothing
end

"""
    isunbounded(rules::BindingRules; maxlength=chainstatebound(rules) + 1)

Whether `rules` admits arbitrarily large structures. See [`growthwitness`](@ref), which returns
the repeating motion behind a `true`.
"""
isunbounded(rules::BindingRules; kwargs...) = growthwitness(rules; kwargs...) !== nothing

"""
    growthwitness(rules::BindingRules; maxlength=chainstatebound(rules) + 1)

Find a rigid motion `g` and a polyform `C` such that `C, g(C), g²(C), …` is an infinite valid
structure, or `nothing` if none exists.

  - returns `(; generator, cell, period)`: the motion, the particles it repeats, and how many
    particles that is

Chains are grown one particle at a time, tracking each particle's state — its species, the site
its incoming bond uses, and that bond's phase. There are [`chainstatebound`](@ref) such states, so
a longer chain repeats one, and the motion `g` between the two occurrences repeats forever.
Restricting to chains loses nothing: any connected part of a valid structure is valid, and an
infinite structure's bond graph contains an infinite path.

By Chasles' theorem `g` is a screw — a rotation by `θ` about an axis together with a translation
`h` along it, degenerating to a pure rotation (`h = 0`) or a pure translation (`θ = 0`). The two
cases settle the question:

  - `h = 0`: every copy stays the same distance from the axis, so all of them lie in one bounded
    region. Infinitely many particles of positive volume cannot, so growth stops. In 2D this is
    every rotation, which is why the plane has no unbounded non-periodic growth
  - `h ≠ 0`: copies `m` periods apart sit `m·h` apart along the axis, so once `|m·h|` exceeds the
    cell's own axial extent they cannot touch. Only `m ≤ ⌈L/|h|⌉` placements can collide, and
    checking those settles the whole infinite chain

Both bounds are exact rather than cutoffs, so `false` is a proof and not a budget running out.
The shapes never enter: a particle's reach along the axis is bounded by its own radius, so
non-convex pieces can interlock only within the same window. What the guarantee does inherit is
the exactness of pairwise `overlap` for the species in play.
"""
function growthwitness(rules::BindingRules; maxlength::Integer=chainstatebound(rules) + 1)
    for i in 1:nspecies(rules)
        w = _walkchain(Polyform(rules, i), [(i, 0, 0)], maxlength)
        w === nothing || return w
    end
    return nothing
end

# Depth-first over directed chains, extending only at the particle last added, so the states seen
# so far are exactly the chain read from its start.
function _walkchain(poly::Polyform, states, maxlength)
    rules = bindingrules(poly)
    n = nparticles(poly)
    part = poly.particles[end]
    for k in 1:nsites(part, rules)
        site = bindingsite(part, rules, k)
        _isbound_vertex(poly, part, first(site.vertices); canonidxs=false) && continue
        isinert(rules, color(site)) && continue
        for loc in possible_attachments(rules, color(site))
            mate = bindingsite(rules, loc)
            for r in 0:(_ndistincttwists(site, mate)-1)
                child = copy(poly)
                ismissing(raise!(child, site, loc, r)) && continue
                # `raise!` forms every geometric contact, so growth may cross-link back onto the
                # chain. That is a valid structure and often the interesting one -- the walk only
                # insists that one particle was added, and always extends the newest
                nparticles(child) == n + 1 || continue

                state = (loc.species, loc.site, Int(r))
                j = findlast(==(state), states)
                if j !== nothing
                    w = _screwwitness(child, j, n + 1)
                    w === nothing || return w
                end
                if n + 1 < maxlength
                    w = _walkchain(child, push!(copy(states), state), maxlength)
                    w === nothing || return w
                end
            end
        end
    end
    return nothing
end

# The stretch between two occurrences of a state repeats under `g`. Whether that goes on forever
# is decided by the screw's pitch, and if it can, by finitely many placements.
function _screwwitness(poly::Polyform, i::Integer, j::Integer)
    rules = bindingrules(poly)
    cell = poly.particles[i:(j-1)]
    g = poly.particles[j].pose * inv(poly.particles[i].pose)

    axis, h = _screwpitch(g)
    h === nothing && return nothing               # a pure rotation cannot escape its own annulus

    us = (dot(p.pose.x, axis) for p in cell)
    rmax = maximum(bounding_radius(species(rules, speciesindex(p))) for p in cell)
    extent = maximum(us) - minimum(us) + 2 * rmax
    gm = g
    for _ in 1:ceil(Int, extent/abs(h))
        for p in cell
            moved = gm * p
            first(_overlap_and_contacts(cell, moved, rules)) === true && return nothing
        end
        gm = gm * g
    end
    return (; generator=g, cell=cell, period=length(cell))
end

# `(axis, pitch)` of a rigid motion, or `(nothing, nothing)` when it has no translation along its
# axis: a pure rotation in 3D, any rotation in 2D, or the identity.
function _screwpitch(g::Pose{D,F}) where {D,F}
    tol = sqrt(eps(F))
    # a composed motion accumulates its angle, so a full turn comes back as 2π rather than 0;
    # wrap to (-π, π] before asking whether it turns at all
    turning = abs(rem(rotation_angle(g.psi), 2 * F(π), RoundNearest)) > tol
    if !turning
        n = norm(g.x)
        n > tol || return nothing, nothing        # the identity: no repeat at all
        return g.x / n, n                         # a pure translation is a screw of zero angle
    end
    D == 3 || return nothing, nothing             # in the plane a turn has no axis to climb
    axis = rotation_axis(g.psi)
    pitch = dot(g.x, axis)
    return axis, abs(pitch) > tol ? pitch : nothing
end

"""
    chainstatebound(rules::BindingRules)

How long a chain of `rules` can get before it must repeat itself.

What a chain can do next depends only on the species at its end, which of that species' sites
carries the incoming bond, and with which phase — finitely many states. A longer chain visits one
twice, and the stretch between the two visits is a cell that repeats forever, so searching past
this length can find no periodic chain that a shorter one would have missed.
"""
function chainstatebound(rules::BindingRules)
    total = 0
    for i in 1:nspecies(rules)
        ps = species(rules, i)
        for k in 1:nsites(ps)
            site = bindingsite(ps, k)
            phases = 1
            for loc in possible_attachments(rules, color(site))
                mate = bindingsite(rules, loc)
                phases = max(phases, _ndistincttwists(mate, site))
            end
            total += phases
        end
    end
    return total
end

# How many of a particle's sites are bonded; a chain's ends are the particles with at most one.
function _bonddegree(poly::Polyform, part)
    rules = bindingrules(poly)
    return count(1:nsites(part, rules)) do k
        _isbound_vertex(poly, part, first(bindingsite(part, rules, k).vertices); canonidxs=false)
    end
end

# Every chain one particle longer, grown at an end. `raise!` forms all geometric contacts, so an
# attachment can branch the chain or close it into a ring; those are dropped, leaving paths.
function _extendends(poly::Polyform)
    rules = bindingrules(poly)
    out = typeof(poly)[]
    for part in poly.particles
        nparticles(poly) == 1 || _bonddegree(poly, part) <= 1 || continue
        for k in 1:nsites(part, rules)
            site = bindingsite(part, rules, k)
            _isbound_vertex(poly, part, first(site.vertices); canonidxs=false) && continue
            isinert(rules, color(site)) && continue
            for loc in possible_attachments(rules, color(site))
                mate = bindingsite(rules, loc)
                for r in 0:(_ndistincttwists(site, mate)-1)
                    child = copy(poly)
                    ismissing(raise!(child, site, loc, r)) && continue
                    _ischain(child) && push!(out, child)
                end
            end
        end
    end
    return out
end

# A path: every particle bonded to at most two others, and no cycle.
function _ischain(poly::Polyform)
    nbonds = sum(part -> _bonddegree(poly, part), poly.particles; init=0) ÷ 2
    return nbonds == nparticles(poly) - 1 &&
           all(part -> _bonddegree(poly, part) <= 2, poly.particles)
end

# Every recorded contact joins two sites that interact -- `_overlap_and_contacts` refuses the
# others before they reach here -- so their color pair is always in the bond list.
function _bondtype(rules, basecolors, (i1, i2))
    i = findfirst(==(minmax(basecolors[i1], basecolors[i2])), bonded_colors(rules))
    isnothing(i) && error("Internal error: a recorded contact has no bond type. Please file an issue.")
    return i
end

# Translations that lay an open site onto a compatible, already-aligned partner site. Both signs
# appear (the site pair swaps), so the shell search only needs non-negative coefficients.
function _candidatelatticevectors(sites, rules)
    intmat = interactionmatrix(rules)
    vecs = typeof(first(sites).pose.x)[]
    for s1 in sites, s2 in sites
        intmat[color(s1), color(s2)] || continue
        isaligned(s1, s2) || continue
        v = s1.pose.x - s2.pose.x
        norm(v) < _tol(v) && continue
        any(u -> u ≈ v, vecs) || push!(vecs, v)
    end
    return vecs
end

# Depth-first search over strictly growing vector index sets: place the newest vector's shells,
# record every consistent configuration, recurse, undo.
function _tilings!(record!::F, s::_ShellSearch, maxvecs::Integer, from::Integer) where {F}
    length(s.chosen) == maxvecs && return
    for idx in from:length(s.vectors)
        push!(s.chosen, idx)
        n0 = length(s.contacts)
        added = _placeshell!(s)
        # a shell without a single bond only builds disconnected unions; skip it
        if added !== nothing && length(s.contacts) > n0
            record!()
            _tilings!(record!, s, maxvecs, idx + 1)
        end
        if added !== nothing
            foreach(i -> s.consumed[i] = false, added)
            resize!(s.contacts, n0)
            pop!(s.placed)
        end
        pop!(s.chosen)
    end
    return
end

# Place every copy whose coefficient on the newest vector is positive (earlier-vector shells were
# placed by earlier stages). Returns the newly consumed site ranges, or `nothing` on an overlap,
# an invalid or bound-site contact, or a doubly-consumed site — with all bookkeeping undone.
# Returns the indices of the sites it consumed.
function _placeshell!(s::_ShellSearch)
    parts = s.cell.meta.particles
    k = length(s.chosen)
    n0 = length(s.contacts)
    added = Int[]
    shell = eltype(parts)[]
    function fail()
        foreach(i -> s.consumed[i] = false, added)
        resize!(s.contacts, n0)
        return nothing
    end

    for m in Iterators.product(ntuple(_ -> 0:(s.nreps), k - 1)..., 1:(s.nreps))
        t = sum(m[i] * s.vectors[s.chosen[i]] for i in 1:k)
        for part in parts
            sp = part + t
            for prior in s.placed
                first(_overlap_and_contacts(prior, sp, s.rules)) === true && return fail()
            end
            first(_overlap_and_contacts(shell, sp, s.rules)) === true && return fail()
            ov, cts = _overlap_and_contacts(parts, sp, s.rules)
            ov && return fail()
            for (; vs1, vs2) in cts
                # both endpoints must be sites of the cell, each closed at most once
                i1, i2 = _lookup(s, vs1), _lookup(s, vs2)
                (i1 == 0 || i2 == 0) && return fail()
                s.consumed[i1] && return fail()
                s.consumed[i1] = true
                push!(added, i1)
                if i2 != i1
                    s.consumed[i2] && return fail()
                    s.consumed[i2] = true
                    push!(added, i2)
                end
                push!(s.contacts, (i1, i2))
            end
            push!(shell, sp)
        end
    end
    push!(s.placed, shell)
    return added
end
