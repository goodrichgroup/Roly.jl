# Periodic witnesses: translated copies of a polyform closing open binding sites into valid
# bonds realize an infinite chain or crystal, whose per-cell composition no finite structure
# reaches. A tiling need not close every site — an infinite chain of squares closes one axis and
# leaves the other open — so partial closures are enumerated too, with completeness a property.

"""
    tilings(poly::Polyform; nreps=2)

Enumerate the periodic closures of `poly`: choices of up to `dimension` translation vectors
under which copies of `poly` form valid bonds with each other and never overlap, each vector
contributing at least one bond.

  - `nreps`: neighbor shells placed and checked per vector; a cutoff, like every witness search
  - returns a vector of `(; vectors, bondtypes, complete)`: the generating vectors, the bond
    type of every bond one unit cell contributes, and whether all open sites were closed

Candidate vectors connect interacting, aligned pairs of open sites.
"""
function tilings(poly::Polyform; nreps::Integer=2)
    sys = bindingrules(poly)
    opensites = collect_open_bindingsites(poly)
    V = fieldtype(posetype(sys), :x)
    out = NamedTuple{(:vectors, :bondtypes, :complete),Tuple{Vector{V},Vector{Int},Bool}}[]
    isempty(opensites) && return out
    vecs = _candidatelatticevectors(opensites, sys)
    isempty(vecs) && return out

    siteof = Dict(s.vertices => s for s in opensites)
    state = (consumed=Set{UnitRange{Int}}(), contacts=NTuple{2,UnitRange{Int}}[],
             placed=Vector{eltype(poly.particles)}[], chosen=Int[])
    record!() = push!(out, (vectors=vecs[state.chosen],
                            bondtypes=[_bondtype(sys, siteof, c) for c in state.contacts],
                            complete=length(state.consumed) == length(siteof)))
    _tilings!(record!, state, poly, sys, vecs, siteof, nreps, dimension(sys), 1)
    return out
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
function _tilings!(record!::F, state, poly, sys, vecs, siteof, nreps, maxvecs, from) where {F}
    length(state.chosen) == maxvecs && return
    for idx in from:length(vecs)
        push!(state.chosen, idx)
        n0 = length(state.contacts)
        added = _placeshell!(state, poly, sys, vecs, siteof, nreps)
        # a shell without a single bond only builds disconnected unions; skip it
        if added !== nothing && length(state.contacts) > n0
            record!()
            _tilings!(record!, state, poly, sys, vecs, siteof, nreps, maxvecs, idx + 1)
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
function _placeshell!(state, poly, sys, vecs, siteof, nreps)
    k = length(state.chosen)
    n0 = length(state.contacts)
    added = UnitRange{Int}[]
    shell = eltype(poly.particles)[]
    function fail()
        foreach(r -> delete!(state.consumed, r), added)
        resize!(state.contacts, n0)
        return nothing
    end

    for m in Iterators.product(ntuple(_ -> 0:nreps, k - 1)..., 1:nreps)
        t = sum(m[i] * vecs[state.chosen[i]] for i in 1:k)
        for part in poly.particles
            sp = part + t
            for prior in state.placed
                first(_overlap_and_contacts(prior, sp, sys)) === true && return fail()
            end
            first(_overlap_and_contacts(shell, sp, sys)) === true && return fail()
            ov, cts = _overlap_and_contacts(poly.particles, sp, sys)
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
