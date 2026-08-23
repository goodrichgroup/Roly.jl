mutable struct PolyformAux{BS<:BindingSite}
    seen::Set{NautyDiGraph}
    pairs::Vector{Tuple{BS,BindingSiteLoc,Int}}
end
Base.copy(polyaux::PolyformAux) = typeof(polyaux)(copy(polyaux.seen), copy(polyaux.pairs))

function ls!(k::Polyform, s::Polyform)
    copy!(k, s)
    return lower!(k)
end

function adj!(u::Polyform, v::Polyform, j::Integer, aux::PolyformAux)
    if nparticles(v) == 0
        j > nspecies(bindingrules(v)) && return nothing
        return copy!(u, Polyform(bindingrules(v), j))
    end

    j == 1 && collect_compatible_pairs!(aux.pairs, v)
    j > length(aux.pairs) && return nothing

    site, siteloc, r = aux.pairs[j]
    copy!(u, v)
    out = raise!(u, site, siteloc, r)
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

All other keyword arguments are passed to the underlying `reversesearch` routine.

Returns a `NamedTuple` `(; nstructures, largest_size, status)`, where `status` is the `RSStatus` of
the underlying reverse search: `Finished` if the enumeration ran to completion, and otherwise
`MaxDepthReached`, `MaxVerticesReached`, or `BreakTriggered` to indicate why it stopped early.
"""
function polyenum(f, rules::BindingRules; maxsize=Inf, maxstrs=Inf, kwargs...)
    v₀ = Polyform(rules)
    BS = BindingSite{posetype(rules),numtype(rules)}
    aux = PolyformAux{BS}(Set{NautyDiGraph}(), Tuple{BS,BindingSiteLoc,Int}[])
    rsys = RSSystem(ls!, adj!, v₀; aux)

    frs = isnothing(f) ? nothing : (v, _) -> f(v, nparticles(v))
    result = reversesearch(frs, rsys; maxdepth=maxsize, maxverts=maxstrs+1, kwargs...)
    return (; nstructures=result.nvertices - 1, largest_size=result.depth_reached, status=result.status)
end
polyenum(sys::BindingRules; kwargs...) = polyenum(nothing, sys; kwargs...)

"""
    polygen(sys::BindingRules; kwargs...)

Enumerate all polyforms allowed by `sys` and return them as a `Vector`, sorted by number
of particles. Keyword arguments are forwarded to [`polyenum`](@ref).
"""
function polygen(sys; kwargs...)
    v₀ = Polyform(sys)
    strs = typeof(v₀)[]
    f(s, args...) = (push!(strs, copy(s)); true)
    polyenum(f, sys; kwargs...)
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
function _count_upto_budget(sys::BindingRules; maxsize, budget)
    counts = Int[]
    largest = 0
    while length(counts) < maxsize
        res = polyenum(sys; maxsize=length(counts)+1, maxstrs=budget)
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
function _estimatecount(sys::BindingRules; pkeep, depth, n0, maxsize, maxsamples, rng)
    estimate = Ref(Float64(n0))
    largest = Ref(0)

    function f(_, n)
        n > largest[] && (largest[] = n)
        n <= depth && return ACCEPT
        rand(rng) < pkeep || return REJECT
        estimate[] += pkeep^(-(n - depth))
        return ACCEPT
    end

    res = polyenum(f, sys; maxsize, maxstrs=n0+maxsamples)
    return (; estimate=estimate[], largest=largest[], res.status)
end

# A subsampled search that dies out before reaching `maxsize` misses the deepest structures entirely,
# and one that explodes is unaffordable. Optimize `pkeep` on throwaway pilot runs so the estimate more likely
# reaches the bottom of the search tree, without oversampling.
function _calibrate_pkeep(sys::BindingRules; pkeep, depth, n0, maxsize, maxsamples, rng, npilots=4, pmax=0.95)
    pkeep = clamp(pkeep, eps(float(pkeep)), pmax)
    for _ in 1:npilots
        pkeep >= pmax && return pmax
        pilot = _estimatecount(sys; pkeep, depth, n0, maxsize, maxsamples, rng)
        if pilot.status == MaxVerticesReached # if we reached maxsamples, reduce pkeep
            pkeep = pkeep / 2
        elseif pilot.largest < maxsize # if we didnt reach the bottom of the tree, increase pkeep
            pkeep = min(pmax, 2pkeep)
        else
            return pkeep
        end
    end
    return min(pkeep, pmax)
end

"""
    countpolyforms(sys::BindingRules; maxsize=Inf, exact_budget=5000, kwargs...)

Estimate the number of polyforms allowed by `sys`, returning a [`PolyformCount`](@ref).

The count is exact whenever the enumeration can fit into `exact_budget`, and is otherwise estimated by randomly 
subsampling the reverse-search tree beyond the budget. Only polyforms of at most `maxsize` particles are counted, 
and finite `maxsize` is required for inexact counting.

- `maxsize=Inf`: only count structures with at most this many particles.
- `exact_budget=5000`: how many structures to enumerate exactly before switching to estimation.
- `ntrials=5`: number of independent estimates to average.
- `eta=1.2`: mean number of surviving offspring per subsampled structure, setting the initial keep
  probability `pkeep = η/branching`. Must exceed 1, and is refined automatically.
- `pkeep=nothing`: set the keep probability directly, disabling `eta` and the automatic refinement.
- `maxsamples=10^6`: how many structures a single trial may generate beyond the exact budget.
- `rng=Random.default_rng()`: generator used for the subsampling.

References:

  - Knuth, D. E. Estimating the Efficiency of Backtrack Programs. Math. Comp. 29, 122-136 (1975).
    [doi:10.1090/S0025-5718-1975-0373371-6](https://doi.org/10.1090/S0025-5718-1975-0373371-6)
  - McKay, B. D. Isomorph-Free Exhaustive Generation. J. Algorithms 26, 306-324 (1998).
    [doi:10.1006/jagm.1997.0898](https://doi.org/10.1006/jagm.1997.0898)

See also [`polyenum`](@ref), [`polygen`](@ref).
"""
function countpolyforms(
    sys::BindingRules;
    maxsize=Inf,
    exact_budget=5_000,
    ntrials=5,
    eta=1.2,
    pkeep=nothing,
    maxsamples=10^6,
    rng=Random.default_rng(),
)
    budget = max(exact_budget, nspecies(sys) + nbonds(sys) + 1) # ensure we visit at least all monomers and dimers
    counts, largest, status = _count_upto_budget(sys; maxsize, budget)

    if status != MaxVerticesReached
        return PolyformCount(counts[end], largest, status == MaxDepthReached)
    end

    if isinf(maxsize)
        throw(
            ArgumentError(
                "`sys` allows structures of unbounded size, of which there are " *
                "infinitely many. Pass an explicit `maxsize` to count the structures " *
                "up to a given number of particles.",
            ),
        )
    end

    depth, n0 = length(counts), counts[end]
    if isnothing(pkeep)
        counts_per_size = diff([0; counts])
        branching = counts_per_size[end] / counts_per_size[end - 1]
        pkeep = _calibrate_pkeep(sys; pkeep=max(eta, 1)/branching, depth, n0, maxsize, maxsamples, rng)
    elseif !(0 < pkeep <= 1)
        throw(ArgumentError("pkeep=$pkeep must be between 0 and 1."))
    end

    trials = [_estimatecount(sys; pkeep, depth, n0, maxsize, maxsamples, rng) for _ in 1:ntrials]
    largest = maximum(t.largest for t in trials)
    samples = [t.estimate for t in trials if t.status != MaxVerticesReached] # dont count trials that ran out of nodes

    if isempty(samples)
        throw(
            ErrorException(
                "all $ntrials trials exceeded maxsamples=$maxsamples; raise " * "`maxsamples` or lower `maxsize`"
            ),
        )
    end

    return PolyformCount(samples, largest, true)
end
