# Periodic witnesses: translated copies of a polyform closing open binding sites into valid
# bonds realize an infinite chain or crystal, whose per-cell composition no finite structure
# reaches. A tiling need not close every site — an infinite chain of squares closes one axis and
# leaves the other open — so partial closures are enumerated too, with completeness a property.

"""
    tilings(poly::Polyform; nreps=2, maxblock=1)

Enumerate the periodic closures of `poly`: choices of up to `dimension` translation vectors
under which copies of a cell form valid bonds with each other and never overlap, each vector
contributing at least one bond.

  - `nreps`: neighbor shells placed and checked per vector; a cutoff, like every witness search
  - `maxblock`: how many copies of `poly` a cell may contain. With `1` the cell is `poly` itself
    and only its translations are tried; larger values also try the blocks the assembly machinery
    can build out of `poly`, which is what admits a cell whose copies are *rotated* relative to
    each other
  - returns a vector of `(; vectors, bondtypes, complete, order)`: the generating vectors, the
    bond type of every bond one cell contributes — those inside the block included — whether all
    open sites were closed, and how many copies of `poly` the cell holds

A cell's composition is `order` copies of `poly` plus `bondtypes`.

Candidate vectors connect interacting, aligned pairs of open sites, so a closure that needs a
turn is found only through the block it makes: growing blocks with [`MetaParticleSpecies`](@ref)
puts the turn inside the cell, leaving a lattice that is a pure translation again.
"""
function tilings(poly::Polyform; nreps::Integer=2, maxblock::Integer=1)
    sys = bindingrules(poly)
    V = fieldtype(posetype(sys), :x)
    T = NamedTuple{(:vectors, :bondtypes, :complete, :order),
                   Tuple{Vector{V},Vector{Int},Bool,Int}}
    out = T[]
    for (parts, opensites, prebonds, order) in _tileblocks(poly, maxblock)
        _addtilings!(out, parts, opensites, prebonds, order, sys, nreps)
    end
    return out
end

# `sites` are every site the block exposes, `prebonds` those of them already joined to each other
# inside the block. Both are needed: the bonds count towards the cell's composition, and their
# endpoints are spent, so only the rest can carry the lattice.
function _addtilings!(out, parts, sites, prebonds, order, sys, nreps)
    consumed = Set{UnitRange{Int}}()
    for (r1, r2) in prebonds
        push!(consumed, r1)
        push!(consumed, r2)
    end
    free = [s for s in sites if s.vertices ∉ consumed]
    isempty(free) && return out
    vecs = _candidatelatticevectors(free, sys)
    isempty(vecs) && return out

    siteof = Dict(s.vertices => s for s in sites)
    pretypes = [_bondtype(sys, siteof, c) for c in prebonds]
    state = (consumed=consumed, contacts=NTuple{2,UnitRange{Int}}[],
             placed=Vector{eltype(parts)}[], chosen=Int[])
    record!() = push!(out, (vectors=vecs[state.chosen],
                            bondtypes=vcat(pretypes,
                                           [_bondtype(sys, siteof, c) for c in state.contacts]),
                            complete=length(state.consumed) == length(siteof),
                            order=order))
    _tilings!(record!, state, parts, sys, vecs, siteof, nreps, dimension(sys), 1)
    return out
end

# The tiles to try: `poly` itself, then every block of up to `maxblock` copies that the assembly
# machinery can build out of it. Growing blocks through `MetaParticleSpecies` is what lets a cell
# be a *rotated* pair of copies -- the translation search that follows only ever considers
# candidates that are already aligned, so a cell needing a turn is invisible to it otherwise.
function _tileblocks(poly::Polyform, maxblock::Integer)
    sys = bindingrules(poly)
    inner = collect_open_bindingsites(poly)
    S = eltype(inner)
    P = eltype(poly.particles)
    blocks = Tuple{Vector{P},Vector{S},Vector{NTuple{2,UnitRange{Int}}},Int}[]
    push!(blocks, (poly.particles, inner, NTuple{2,UnitRange{Int}}[], 1))
    (maxblock <= 1 || isempty(inner)) && return blocks

    metasys = _blockrules(poly, inner, sys)
    metasys === nothing && return blocks
    for mb in polygen(metasys; maxsize=maxblock)
        nparticles(mb) == 1 && continue          # that is `poly` itself, already listed
        blk = _unwrapblock(mb, poly, inner, sys)
        blk === nothing || push!(blocks, blk)
    end
    return blocks
end

# Blocks bond to blocks exactly where the underlying sites do, so the meta-rules are read off the
# cluster's own interaction matrix.
function _blockrules(poly, inner, sys)
    mp = MetaParticleSpecies(poly)
    intmat = interactionmatrix(sys)
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
function _unwrapblock(mb::Polyform, poly::Polyform, inner, sys)
    nvc = nv(graphrep(poly))
    P = eltype(poly.particles)
    parts = P[]
    sites = eltype(inner)[]
    prebonds = NTuple{2,UnitRange{Int}}[]
    for (b, mpart) in enumerate(mb.particles)
        off = (b - 1) * nvc
        pose = mpart.pose
        copyparts = [P(pose * p.pose, p.leading_vertex + off, p.species_index)
                     for p in poly.particles]
        for sp in copyparts
            ov, cts = _overlap_and_contacts(parts, sp, sys)
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
    f(s, _) = isunitcell(s; kwargs...) ? (hit[] = copy(s); BREAK) : ACCEPT
    polyenum(f, rules; maxsize=maxtilesize)
    return hit[]
end

_bondtype(sys, siteof, (r1, r2)) =
    findfirst(==(minmax(color(siteof[r1]), color(siteof[r2]))), sys._bondlist)

# Translations that lay an open site onto a compatible, already-aligned partner site. Both signs
# appear (the site pair swaps), so the shell search only needs non-negative coefficients.
function _candidatelatticevectors(opensites, sys)
    intmat = interactionmatrix(sys)
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
function _tilings!(record!::F, state, parts, sys, vecs, siteof, nreps, maxvecs, from) where {F}
    length(state.chosen) == maxvecs && return
    for idx in from:length(vecs)
        push!(state.chosen, idx)
        n0 = length(state.contacts)
        added = _placeshell!(state, parts, sys, vecs, siteof, nreps)
        # a shell without a single bond only builds disconnected unions; skip it
        if added !== nothing && length(state.contacts) > n0
            record!()
            _tilings!(record!, state, parts, sys, vecs, siteof, nreps, maxvecs, idx + 1)
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
function _placeshell!(state, parts, sys, vecs, siteof, nreps)
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
                first(_overlap_and_contacts(prior, sp, sys)) === true && return fail()
            end
            first(_overlap_and_contacts(shell, sp, sys)) === true && return fail()
            ov, cts = _overlap_and_contacts(parts, sp, sys)
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
