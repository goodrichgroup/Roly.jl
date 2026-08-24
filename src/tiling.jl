# Periodic witnesses: translated copies of a polyform closing open binding sites into valid
# bonds realize an infinite chain or crystal, whose per-cell composition no finite structure
# reaches. A tiling need not close every site — an infinite chain of squares closes one axis and
# leaves the other open — so partial closures are enumerated too, with completeness a property.

"""
    tilings(poly::Polyform; nreps=2, maxblock=1)

Enumerate the periodic closures of `poly`: choices of up to `dimension` translation vectors
under which copies of a cell form valid bonds with each other and never overlap, each vector
contributing at least one bond.

  - `nreps`: neighbor shells placed and checked per vector.
  - `maxblock`: construct meta-polyforms up to size `maxblock`, which are then checked for translation tilings.
    `maxblock > 1` makes it possible to find tilings where `poly` appears in rotated configurations.
  - returns a vector of `(; vectors, bondtypes, complete, order)`: the generating vectors, the
    bond type of every bond one cell contributes (those inside the block included), whether all
    open sites were closed, and how many copies of `poly` the cell holds

A cell's composition is `order` copies of `poly` plus `bondtypes`.

Candidate vectors connect interacting, aligned pairs of open sites, so a closure that needs a
turn is found only through the block it makes: growing blocks with [`MetaParticleSpecies`](@ref)
puts the turn inside the cell, leaving a lattice that is a pure translation again.
"""
function tilings(poly::Polyform; nreps::Integer=2, maxblock::Integer=1)
    rules = bindingrules(poly)
    V = fieldtype(posetype(rules), :x)
    T = NamedTuple{(:vectors, :bondtypes, :complete, :order),
        Tuple{Vector{V},Vector{Int},Bool,Int}}
    out = T[]
    for (parts, opensites, prebonds, order) in _tileblocks(poly, maxblock)
        _addtilings!(out, parts, opensites, prebonds, order, rules, nreps)
    end
    return out
end

# `sites` are every site the block exposes, `prebonds` those of them already joined to each other
# inside the block. Both are needed: the bonds count towards the cell's composition, and their
# endpoints are spent, so only the rest can carry the lattice.
function _addtilings!(out, parts, sites, prebonds, order, rules, nreps)
    consumed = Set{UnitRange{Int}}()
    for (r1, r2) in prebonds
        push!(consumed, r1)
        push!(consumed, r2)
    end
    free = [s for s in sites if s.vertices ∉ consumed]
    isempty(free) && return out
    vecs = _candidatelatticevectors(free, rules)
    isempty(vecs) && return out

    siteof = Dict(s.vertices => s for s in sites)
    pretypes = [_bondtype(rules, siteof, c) for c in prebonds]
    state = (consumed=consumed, contacts=NTuple{2,UnitRange{Int}}[],
        placed=Vector{eltype(parts)}[], chosen=Int[])
    record!() = push!(out, (vectors=vecs[state.chosen],
        bondtypes=vcat(pretypes,
            [_bondtype(rules, siteof, c) for c in state.contacts]),
        complete=length(state.consumed) == length(siteof),
        order=order))
    _tilings!(record!, state, parts, rules, vecs, siteof, nreps, dimension(rules), 1)
    return out
end

# The tiles to try: `poly` itself, then every block of up to `maxblock` copies that the assembly
# machinery can build out of it.
function _tileblocks(poly::Polyform, maxblock::Integer)
    rules = bindingrules(poly)
    inner = collect_open_bindingsites(poly)
    S = eltype(inner)
    P = eltype(poly.particles)
    blocks = Tuple{Vector{P},Vector{S},Vector{NTuple{2,UnitRange{Int}}},Int}[]
    push!(blocks, (poly.particles, inner, NTuple{2,UnitRange{Int}}[], 1))
    (maxblock <= 1 || isempty(inner)) && return blocks

    metasys = _blockrules(poly, inner, rules)
    metasys === nothing && return blocks
    for mb in polygen(metasys; maxsize=maxblock)
        nparticles(mb) == 1 && continue          # that is `poly` itself, already listed
        blk = _unwrapblock(mb, poly, inner, rules)
        blk === nothing || push!(blocks, blk)
    end
    return blocks
end

# Blocks bond to blocks exactly where the underlying sites do, so the meta-rules are read off the
# cluster's own interaction matrix.
function _blockrules(poly, inner, rules)
    mp = MetaParticleSpecies(poly)
    intmat = interactionmatrix(rules)
    rows = Vector{Int}[]
    for i in eachindex(inner), j in i:length(inner)
        intmat[color(inner[i]), color(inner[j])] || continue
        push!(rows, [1, i, 1, j])
    end
    isempty(rows) && return nothing
    return BindingRules(permutedims(reduce(hcat, rows)), mp)
end

# Lay a block's copies of `poly` out as plain particles, renumbering each copy's vertices so every
# site keeps a range of its own, and collect the bonds the copies make with each other.
function _unwrapblock(mb::Polyform, poly::Polyform, inner, rules)
    nvc = nv(graphrep(poly))
    P = eltype(poly.particles)
    parts = P[]
    sites = eltype(inner)[]
    prebonds = NTuple{2,UnitRange{Int}}[]
    for (b, mpart) in enumerate(mb.particles)
        off = (b - 1) * nvc
        pose = mpart.pose
        copyparts = [P(pose * p.pose, p.leadingvertex + off, p.speciesindex)
                     for p in poly.particles]
        for sp in copyparts
            ov, cts = _overlap_and_contacts(parts, sp, rules)
            ov && return nothing        # a valid meta-assembly should never land here
            # a contact also carries its phase and period, which the cell does not need
            append!(prebonds, ((c[1], c[2]) for c in cts))
        end
        append!(parts, copyparts)
        append!(sites, (shift_vertices(pose * s, off) for s in inner))
    end
    return parts, sites, prebonds, nparticles(mb)
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
continuation gives a valid cluster of every length. This is far cheaper than [`cantile`](@ref),
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
`maxblock ≥ q`, but an irrational one never closes at all, and no periodicity test can see it.
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

Find a rigid motion `g` and a cluster `C` such that `C, g(C), g²(C), …` is an infinite valid
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
        for siteloc in possible_attachments(rules, color(site))
            mate = bindingsite(species(rules, siteloc[1]), siteloc[2])
            for r in 0:(_ndistincttwists(site, mate)-1)
                child = copy(poly)
                ismissing(raise!(child, site, siteloc, r)) && continue
                # `raise!` forms every geometric contact, so growth may cross-link back onto the
                # chain. That is a valid structure and often the interesting one -- the walk only
                # insists that one particle was added, and always extends the newest
                nparticles(child) == n + 1 || continue

                state = (siteloc[1], siteloc[2], Int(r))
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
                mate = bindingsite(species(rules, loc[1]), loc[2])
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
            for siteloc in possible_attachments(rules, color(site))
                mate = bindingsite(species(rules, siteloc[1]), siteloc[2])
                for r in 0:(_ndistincttwists(site, mate)-1)
                    child = copy(poly)
                    ismissing(raise!(child, site, siteloc, r)) && continue
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

_bondtype(rules, siteof, (r1, r2)) =
    findfirst(==(minmax(color(siteof[r1]), color(siteof[r2]))), rules._bondlist)

# Translations that lay an open site onto a compatible, already-aligned partner site. Both signs
# appear (the site pair swaps), so the shell search only needs non-negative coefficients.
function _candidatelatticevectors(opensites, rules)
    intmat = interactionmatrix(rules)
    vecs = typeof(first(opensites).pose.x)[]
    for s1 in opensites, s2 in opensites
        intmat[color(s1), color(s2)] || continue
        isaligned(s1, s2) || continue
        v = s1.pose.x - s2.pose.x
        sum(abs2, v) < 1e-18 && continue
        any(u -> isapprox(u, v; atol=1e-9), vecs) || push!(vecs, v)
    end
    return vecs
end

# Depth-first search over strictly growing vector index sets: place the newest vector's shells,
# record every consistent configuration, recurse, undo.
function _tilings!(record!::F, state, parts, rules, vecs, siteof, nreps, maxvecs, from) where {F}
    length(state.chosen) == maxvecs && return
    for idx in from:length(vecs)
        push!(state.chosen, idx)
        n0 = length(state.contacts)
        added = _placeshell!(state, parts, rules, vecs, siteof, nreps)
        # a shell without a single bond only builds disconnected unions; skip it
        if added !== nothing && length(state.contacts) > n0
            record!()
            _tilings!(record!, state, parts, rules, vecs, siteof, nreps, maxvecs, idx + 1)
        end
        if added !== nothing
            foreach(r -> delete!(state.consumed, r), added)
            resize!(state.contacts, n0)
            pop!(state.placed)
        end
        pop!(state.chosen)
    end
    return
end

# Place every copy whose coefficient on the newest vector is positive (earlier-vector shells were
# placed by earlier stages). Returns the newly consumed site ranges, or `nothing` on an overlap,
# an invalid or bound-site contact, or a doubly-consumed site — with all bookkeeping undone.
function _placeshell!(state, parts, rules, vecs, siteof, nreps)
    k = length(state.chosen)
    n0 = length(state.contacts)
    added = UnitRange{Int}[]
    shell = eltype(parts)[]
    function fail()
        foreach(r -> delete!(state.consumed, r), added)
        resize!(state.contacts, n0)
        return nothing
    end

    for m in Iterators.product(ntuple(_ -> 0:nreps, k - 1)..., 1:nreps)
        t = sum(m[i] * vecs[state.chosen[i]] for i in 1:k)
        for part in parts
            sp = part + t
            for prior in state.placed
                first(_overlap_and_contacts(prior, sp, rules)) === true && return fail()
            end
            first(_overlap_and_contacts(shell, sp, rules)) === true && return fail()
            ov, cts = _overlap_and_contacts(parts, sp, rules)
            ov && return fail()
            for (r1, r2) in cts
                # both endpoints must be open sites of the cell, each closed at most once
                (haskey(siteof, r1) && haskey(siteof, r2)) || return fail()
                r1 in state.consumed && return fail()
                push!(state.consumed, r1)
                push!(added, r1)
                if r2 != r1
                    r2 in state.consumed && return fail()
                    push!(state.consumed, r2)
                    push!(added, r2)
                end
                push!(state.contacts, (r1, r2))
            end
            push!(shell, sp)
        end
    end
    push!(state.placed, shell)
    return added
end
