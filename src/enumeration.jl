mutable struct PolyformAux{BS<:BindingSite}
    seen::Set{NautyDiGraph}
    attachments::Vector{Tuple{BS,SpeciesSiteLoc,Int}}
end
Base.copy(polyaux::PolyformAux) = typeof(polyaux)(copy(polyaux.seen), copy(polyaux.attachments))

function ls!(k::Polyform, s::Polyform)
    copy!(k, s)
    return lower!(k)
end

function adj!(u::Polyform, v::Polyform, j::Integer, aux::PolyformAux)
    if nparticles(v) == 0
        j > nspecies(bindingrules(v)) && return nothing
        return copy!(u, Polyform(bindingrules(v), j))
    end

    j == 1 && collect_attachments!(aux.attachments, v)
    j > length(aux.attachments) && return nothing

    site, loc, t = aux.attachments[j]
    copy!(u, v)
    out = raise!(u, site, loc, t)
    (ismissing(out) || isnothing(out)) && return out

    if graphrep(out) ∈ aux.seen
        return missing
    else
        push!(aux.seen, copy(graphrep(out)))
        return out
    end
end

"""
    polyenum([f], rules::BindingRules; maxsize=Inf, maxstrs=Inf, kwargs...)

Iterate over all polyforms allowed by `rules` using _reverse search_.
Polyforms are generated up to size `maxsize`, and the enumeration terminates after `maxstrs`
structures have been generated. Use the function `f` to process generated structures and impose
constraints.

`f(s)` must take a structure `s` and return one of three signals:

- `ACCEPT` (or `true`): enumeration continues as normal.
- `REJECT` (or `false`): the offspring of the current structure will not be generated.
- `BREAK`: the enumeration terminates immediately.

To ensure well-defined behavior, the function `f` may not accept the offspring of structures
that it rejects.

`checkbounded` refuses a rule set that admits arbitrarily large structures, see
[`isbounded`](@ref), by throwing an `ArgumentError`. It defaults to `true` when neither `maxsize`
nor `maxstrs` caps the enumeration.

All other keyword arguments are passed to the underlying `reversesearch` routine.

Returns a `NamedTuple` `(; nstructures, largest_size, status)`, where `status` is the `RSStatus` of
the underlying reverse search: `Finished` if the enumeration ran to completion, and otherwise
`MaxDepthReached`, `MaxVerticesReached`, or `BreakTriggered` to indicate why it stopped early.
"""
function polyenum(f, rules::BindingRules; maxsize=Inf, maxstrs=Inf,
                  checkbounded::Bool=isinf(maxsize) && isinf(maxstrs), kwargs...)
    checkbounded && _checkbounded(rules)
    v₀ = Polyform(rules)
    BS = sitetype(rules)
    aux = PolyformAux{BS}(Set{NautyDiGraph}(), Tuple{BS,SpeciesSiteLoc,Int}[])
    rsys = RSSystem(ls!, adj!, v₀; aux)

    frs = isnothing(f) ? nothing : (v, _) -> f(v, nparticles(v))
    result = reversesearch(frs, rsys; maxdepth=maxsize, maxverts=maxstrs+1, kwargs...)
    return (; nstructures=result.nvertices - 1, largest_size=result.depth_reached, status=result.status)
end
polyenum(rules::BindingRules; kwargs...) = polyenum(nothing, rules; kwargs...)

"""
    polygen(rules::BindingRules; kwargs...)

Return a list of all polyforms allowed by `rules`, sorted by number
of particles. Keyword arguments are forwarded to [`polyenum`](@ref).
"""
function polygen(rules; kwargs...)
    v₀ = Polyform(rules)
    strs = typeof(v₀)[]
    f(s, args...) = (push!(strs, copy(s)); true)
    polyenum(f, rules; kwargs...)
    return sort!(strs; by=Roly.nparticles)
end

"""
    PolyformCount(n::Integer, largest_size, size_truncated)
    PolyformCount(trials::AbstractVector, largest_size, size_truncated)

Represent the number of polyforms allowed by a set of binding rules, or an estimate thereof.

- `n`: the number of allowed structures.
- `exact`: `true` if `n` was obtained by complete enumeration, `false` if it was estimated.
- `uncertainty`: the standard error of the estimate.
- `largest_size`: the size of the largest structure. If not exact, this is a *lower bound* of the true largest size.
- `size_truncated`: `true` if the count covers only structures up to a given cutoff size
- `trials`: the individual estimates that `n` averages; empty if exact.
"""
struct PolyformCount
    n::Float64
    exact::Bool
    uncertainty::Float64
    largest_size::Int
    size_truncated::Bool
    trials::Vector{Float64}
end

function PolyformCount(n::Integer, largest_size, size_truncated)
    PolyformCount(Float64(n), true, 0.0, largest_size, size_truncated, Float64[])
end

function PolyformCount(trials::AbstractVector{<:Real}, largest_size, size_truncated)
    uncertainty = length(trials) > 1 ? std(trials) / sqrt(length(trials)) : NaN
    return PolyformCount(mean(trials), false, uncertainty, largest_size, size_truncated, trials)
end

function Base.show(io::Core.IO, c::PolyformCount)
    bound = c.exact ? "largest=$(c.largest_size)" : "largest≥$(c.largest_size)"
    trunc = c.size_truncated ? ", truncated" : ""
    if c.exact
        return print(io, "PolyformCount[n=$(round(Int, c.n)), exact, $bound$trunc]")
    end
    n = round(c.n; sigdigits=3)
    err = round(c.uncertainty; sigdigits=2)
    return print(io, "PolyformCount[n≈$n ± $err, $bound$trunc]")
end

# Enumerate exactly at sizes 1, 2, 3, ... until number of polyforms reaches `budget`.
# Report the cumulative number of structures at each size, the largest structure seen, and why the enumeration stopped.
function _count_upto_budget(rules::BindingRules; maxsize, budget)
    counts = Int[]
    largest = 0
    while length(counts) < maxsize
        res = polyenum(rules; maxsize=length(counts)+1, maxstrs=budget)
        res.status == MaxVerticesReached && return counts, largest, res.status

        push!(counts, res.nstructures)
        largest = res.largest_size
        res.status == Finished && return counts, largest, res.status
    end
    return counts, largest, MaxDepthReached
end

# Perform one unbiased sample of the number of structures of size ≤ maxsize. Structures up to `depth` are all
# visited and contribute the known count `n0`. Deeper ones are kept with probability `pkeep`, 
# so a structure of size `n` is reached with probability `pkeep^(n-depth)` and is weighted accordingly. 
function _estimatecount(rules::BindingRules; pkeep, depth, n0, maxsize, maxsamples, rng)
    estimate = Ref(Float64(n0))
    largest = Ref(0)

    function f(_, n)
        n > largest[] && (largest[] = n)
        n <= depth && return ACCEPT
        rand(rng) < pkeep || return REJECT
        estimate[] += pkeep^(-(n - depth))
        return ACCEPT
    end

    res = polyenum(f, rules; maxsize, maxstrs=n0+maxsamples)
    return (; estimate=estimate[], largest=largest[], res.status)
end

# Optimize `pkeep` on throwaway pilot runs so the estimate more likely
# reaches the bottom of the search tree, without oversampling.
function _calibrate_pkeep(rules::BindingRules; pkeep, depth, n0, maxsize, maxsamples, rng, npilots=4, pmax=0.95)
    pkeep = clamp(pkeep, eps(float(pkeep)), pmax)
    for _ in 1:npilots
        pkeep >= pmax && return pmax
        pilot = _estimatecount(rules; pkeep, depth, n0, maxsize, maxsamples, rng)
        if pilot.status == MaxVerticesReached # if we reached maxsamples, reduce pkeep
            pkeep = pkeep / 2
        elseif pilot.largest < (isfinite(maxsize) ? maxsize : depth + 1)
            # the sample died out before the bottom of the tree, or with no cap before it got
            # anywhere past the exact part, so increase pkeep
            pkeep = min(pmax, 2pkeep)
        else
            return pkeep
        end
    end
    return min(pkeep, pmax)
end

"""
    polycount(rules::BindingRules; maxsize=Inf, exact_budget=5000, kwargs...)

Estimate the number of polyforms allowed by `rules`, returning a [`PolyformCount`](@ref).

The count is exact whenever the enumeration can fit into `exact_budget`, and is otherwise estimated by randomly
subsampling the reverse-search tree beyond the budget. Only polyforms of at most `maxsize` particles are counted.
With no `maxsize` the rule set has to be bounded, which is checked first, and the estimate then covers every
structure there is.

- `maxsize=Inf`: only count structures with at most this many particles.
- `exact_budget=5000`: how many structures to enumerate exactly before switching to estimation.
- `ntrials=5`: number of independent estimates to average.
- `eta=1.2`: mean number of surviving offspring per subsampled structure, setting the initial keep
  probability `pkeep = η/branching`. Must exceed 1, and is refined automatically.
- `pkeep=nothing`: set the keep probability directly, disabling `eta` and the automatic refinement.
- `maxsamples=10^6`: how many structures a single trial may generate beyond the exact budget.
- `rng=Random.default_rng()`: generator used for the subsampling.
- `checkbounded=isinf(maxsize)`: refuse a rule set that admits arbitrarily large structures, see
  [`isbounded`](@ref), rather than count until the budget runs out.

References:

  - Knuth, D. E. Estimating the Efficiency of Backtrack Programs. Math. Comp. 29, 122-136 (1975).
    [doi:10.1090/S0025-5718-1975-0373371-6](https://doi.org/10.1090/S0025-5718-1975-0373371-6)
  - McKay, B. D. Isomorph-Free Exhaustive Generation. J. Algorithms 26, 306-324 (1998).
    [doi:10.1006/jagm.1997.0898](https://doi.org/10.1006/jagm.1997.0898)

See also [`polyenum`](@ref), [`polygen`](@ref).
"""
function polycount(
    rules::BindingRules;
    maxsize=Inf,
    exact_budget=5_000,
    ntrials=5,
    eta=1.2,
    pkeep=nothing,
    maxsamples=10^6,
    rng=Random.default_rng(),
    checkbounded::Bool=isinf(maxsize),
)
    checkbounded && _checkbounded(rules)
    budget = max(exact_budget, nspecies(rules) + nbonds(rules) + 1) # ensure we visit at least all monomers and dimers
    counts, largest, status = _count_upto_budget(rules; maxsize, budget)

    if status != MaxVerticesReached
        return PolyformCount(counts[end], largest, status == MaxDepthReached)
    end

    # past the budget with no cap, the estimate terminates only on a bounded rule set, which the
    # check has established unless it was waived
    if isinf(maxsize) && !checkbounded
        throw(
            ArgumentError(
                "`rules` may allow structures of unbounded size, of which there are infinitely " *
                "many. Pass `maxsize` to count up to a given number of particles, or " *
                "`checkbounded=true` to establish that the rules are bounded.",
            ),
        )
    end

    depth, n0 = length(counts), counts[end]
    if isnothing(pkeep)
        counts_per_size = diff([0; counts])
        branching = counts_per_size[end] / counts_per_size[end - 1]
        pkeep = _calibrate_pkeep(rules; pkeep=max(eta, 1)/branching, depth, n0, maxsize, maxsamples, rng)
    elseif !(0 < pkeep <= 1)
        throw(ArgumentError("pkeep=$pkeep must be between 0 and 1."))
    end

    trials = [_estimatecount(rules; pkeep, depth, n0, maxsize, maxsamples, rng) for _ in 1:ntrials]
    largest = maximum(t.largest for t in trials)
    samples = [t.estimate for t in trials if t.status != MaxVerticesReached] # dont count trials that ran out of nodes

    if isempty(samples)
        throw(
            ErrorException(
                "all $ntrials trials exceeded maxsamples=$maxsamples; raise " * "`maxsamples` or lower `maxsize`"
            ),
        )
    end

    return PolyformCount(samples, largest, isfinite(maxsize))
end

# What an enumeration without a cap runs first.
function _checkbounded(rules::BindingRules)
    isbounded(rules) && return nothing
    throw(
        ArgumentError(
            "`rules` allows infinitely many structures and cannot be enumerated completely. " *
            "Pass `maxsize` or `maxstrs` to enumerate up to a bound, or set `checkbounded=false` to disable this error.",
        ),
    )
end

"""
    isbounded(rules::BindingRules; maxchainlength=_nchainstates(rules) + 2)

Return `true` if `rules` admits finitely many structures and `false` otherwise.

The answer is computed by enumerating one-dimensional chains and checking if any chain can repeat indefinitely.
`maxchainlength` sets the length of the longest chain tested. Decreasing `maxchainlength` from its default can 
introduce false positive results.
"""
function isbounded(rules::BindingRules; maxchainlength::Integer=_nchainstates(rules) + 2)
    return _repeatingchain(rules; maxchainlength) === nothing
end

# A rigid motion `g` and a polyform `C` such that `C, g(C), g²(C), …` is an infinite valid
# structure, as `(; generator, cell, period)`, or `nothing` if none exists within `maxchainlength`.
function _repeatingchain(rules::BindingRules; maxchainlength::Integer=_nchainstates(rules) + 2)
    for i in 1:nspecies(rules)
        w = _walkchains(Polyform(rules, i), maxchainlength)
        w === nothing || return w
    end
    return nothing
end

# How long a chain of `rules` can get before it must repeat a state. A state is the species at the
# chain's end, which of its sites the incoming bond uses, and with which twist: one per site a bond
# can arrive at, per distinct twist it can arrive with, over every species.
function _nchainstates(rules::BindingRules)
    total = 0
    for i in 1:nspecies(rules)
        ps = species(rules, i)
        for k in 1:nsites(ps)
            site = bindingsite(ps, k)
            phases = 0
            for loc in compatible_speciessitelocs(rules, color(site))
                mate = bindingsite(rules, loc)
                phases = max(phases, _ndistincttwists(mate, site))
            end
            total += phases
        end
    end
    return total
end

# A chain as the search holds it: the polyform, and per particle the state it was added in, the
# species together with the site and twist its incoming bond uses. The first particle's is
# `(species, 0, 0)`, having no incoming bond.
mutable struct _Chain{P<:Polyform}
    poly::P
    states::Vector{NTuple{3,Int}}
end
Base.copy(c::_Chain) = _Chain(copy(c.poly), copy(c.states))
function Base.copy!(dst::_Chain, src::_Chain)
    copy!(dst.poly, src.poly)
    copy!(dst.states, src.states)
    return dst
end

mutable struct _ChainAux{BS}
    attachments::Vector{Tuple{BS,SpeciesSiteLoc,Int}}
end
Base.copy(a::_ChainAux) = _ChainAux(copy(a.attachments))

# Every way to attach a particle at an open site of the newest particle of `poly`: the site, the
# mating site, and the twist.
function _endattachments!(out, poly::Polyform)
    rules = bindingrules(poly)
    part = poly.particles[end]
    empty!(out)
    for k in 1:nsites(part, rules)
        site = bindingsite(part, rules, k)
        _isbound_vertex(poly, part, first(site.vertices); canonidxs=false) && continue
        isinert(rules, color(site)) && continue
        for loc in compatible_speciessitelocs(rules, color(site))
            mate = bindingsite(rules, loc)
            for r in 0:(_ndistincttwists(site, mate) - 1)
                push!(out, (site, loc, r))
            end
        end
    end
    return out
end

# The `j`-th extension of the chain `v` at its newest particle, written into `u`. Attaching forms
# every geometric contact, so a chain may cross-link back onto itself; that is a valid structure,
# and the newest particle is extended regardless. No graph is ever compared, so none is canonized.
function _extendchain!(u::_Chain, v::_Chain, j::Integer, aux::_ChainAux)
    j == 1 && _endattachments!(aux.attachments, v.poly)
    j > length(aux.attachments) && return nothing
    site, loc, r = aux.attachments[j]
    copy!(u.poly, v.poly)
    ismissing(add_particle!(u.poly, site, loc, r; canonize=false)) && return missing
    copy!(u.states, v.states)
    push!(u.states, (loc.species, loc.site, r))
    return u
end

# Depth-first over directed chains of up to `maxchainlength` particles, returning the first repeat
# whose motion escapes. The states a chain has passed through are exactly the chain read from its
# start, so a repeated state is a stretch that repeats under one motion.
function _walkchains(poly::Polyform, maxchainlength)
    # one particle has no incoming bond to repeat. Said here rather than left to the search, whose
    # depth limit of zero is not honored
    maxchainlength > 1 || return nothing
    root = _Chain(poly, [(speciesindex(poly.particles[1]), 0, 0)])
    aux = _ChainAux(Tuple{sitetype(bindingrules(poly)),SpeciesSiteLoc,Int}[])
    witness = Ref{Any}(nothing)
    _backtrack(_extendchain!, root; aux, maxdepth=maxchainlength - 1) do c, _
        n = nparticles(c.poly)
        j = findlast(==(c.states[end]), view(c.states, 1:(n - 1)))
        isnothing(j) && return ACCEPT
        w = _screwescape(c.poly, j, n)
        isnothing(w) && return ACCEPT
        witness[] = w
        return BREAK
    end
    return witness[]
end

# The stretch between two occurrences of a state repeats under `g`. Whether that goes on forever
# is decided by how far the screw climbs, and if it can, by finitely many placements.
function _screwescape(poly::Polyform, i::Integer, j::Integer)
    rules = bindingrules(poly)
    cell = poly.particles[i:(j - 1)]
    g = poly.particles[j].pose * inv(poly.particles[i].pose)

    axis, h = _screwclimb(g)
    isnothing(h) && return nothing               # a pure rotation cannot escape its own annulus

    us = (dot(p.pose.x, axis) for p in cell)
    rmax = maximum(bounding_radius(species(rules, speciesindex(p))) for p in cell)
    extent = maximum(us) - minimum(us) + 2 * rmax
    gm = g
    for _ in 1:ceil(Int, extent / abs(h))
        for p in cell
            moved = gm * p
            first(_overlap_and_contacts(cell, moved, rules)) === true && return nothing
        end
        gm = gm * g
    end
    return (; generator=g, cell=cell, period=length(cell))
end

# `(axis, climb)` of a rigid motion: the direction it carries things along, and how far per
# application. `(nothing, nothing)` when it carries them nowhere: a pure rotation in 3D, any
# rotation in 2D, or the identity.
function _screwclimb(g::Pose{D,F}) where {D,F}
    tol = sqrt(eps(F))
    turning = abs(rem(rotation_angle(g.psi), 2 * F(π), RoundNearest)) > tol
    if !turning
        # with no turn there is no rotation axis, so the axis is just the direction of translation
        n = norm(g.x)
        n > tol || return nothing, nothing        # the identity: no repeat at all
        return g.x / n, n
    end
    D == 3 || return nothing, nothing
    axis = rotation_axis(g.psi)
    climb = dot(g.x, axis)
    return axis, abs(climb) > tol ? climb : nothing
end
