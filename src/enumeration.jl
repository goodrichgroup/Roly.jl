mutable struct PolyformAux{BS<:BindingSite}
    seen::Set{NautyDiGraph}
    attachments::Vector{Tuple{BS,SpeciesSiteLoc,BondTwist}}
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

    site, loc, tw = aux.attachments[j]
    copy!(u, v)
    out = raise!(u, site, loc, tw)
    (ismissing(out) || isnothing(out)) && return out

    if graphrep(out) ∈ aux.seen
        return missing
    else
        push!(aux.seen, copy(graphrep(out)))
        return out
    end
end

"""
    EnumerationStatus

Why an enumeration stopped: `Complete` when it ran out of polyforms, `MaxSizeReached` or
`MaxPolysReached` when it reached the `maxsize` or `maxpolys` cap, and `Interrupted` when the
callback returned `BREAK`.
"""
@enum EnumerationStatus Complete MaxSizeReached MaxPolysReached Interrupted

# Returns the reverse search's status in the words the enumeration was asked in.
function _status(s::RSStatus)
    s == Finished && return Complete
    s == MaxDepthReached && return MaxSizeReached
    s == MaxVerticesReached && return MaxPolysReached
    return Interrupted
end

"""
    PolyformCount

The number of polyforms a set of binding rules allows, counted or estimated, as returned by
[`polyenum`](@ref) and [`polycount`](@ref).

- `n`: the number of polyforms, an `Int` if exact and the estimated mean as a `Float64` otherwise.
- `exact`: `true` if `n` was obtained by complete enumeration, `false` if it was estimated.
- `err`: the standard error of the estimate, `0` if exact.
- `maxsize`: the size of the largest polyform counted, which is the `maxsize` cap whenever the
  enumeration reached it. If estimated, a lower bound.
- `status`: why the enumeration stopped, an [`EnumerationStatus`](@ref).
- `trials`: the individual estimates `n` averages, empty if exact.
"""
struct PolyformCount{T<:Real}
    n::T
    exact::Bool
    err::Float64
    maxsize::Int
    status::EnumerationStatus
    trials::Vector{Float64}
end

"""
    PolyformCount(n::Integer, maxsize, status)

Return the exact count of `n` polyforms, the largest of `maxsize` particles.
"""
PolyformCount(n::Integer, maxsize, status) = PolyformCount(Int(n), true, 0.0, maxsize, status, Float64[])

"""
    PolyformCount(trials::AbstractVector, maxsize, status)

Return the estimate averaging `trials`, with their standard error, the largest polyform seen having
`maxsize` particles.
"""
function PolyformCount(trials::AbstractVector{<:Real}, maxsize, status)
    err = length(trials) > 1 ? std(trials) / sqrt(length(trials)) : NaN
    return PolyformCount(mean(trials), false, err, maxsize, status, trials)
end

function Base.show(io::Core.IO, c::PolyformCount)
    count = c.exact ? "n=$(c.n), exact" : "n≈$(round(c.n; sigdigits=3)) ± $(round(c.err; sigdigits=2))"
    bound = c.exact ? "maxsize=$(c.maxsize)" : "maxsize≥$(c.maxsize)"
    return print(io, "PolyformCount[$count, $bound, $(c.status)]")
end

"""
    polyenum([f], rules::BindingRules; maxsize=Inf, maxpolys=Inf, kwargs...)

Iterate over all polyforms allowed by `rules` using _reverse search_.
Polyforms are generated up to size `maxsize`, and the enumeration terminates after `maxpolys`
polyforms have been generated. Use the function `f` to process generated polyforms and impose
constraints.

`f(p)` must take a polyforms `p` and return one of three signals:

- `ACCEPT` (or `true`): enumeration continues as normal.
- `REJECT` (or `false`): the offspring of the current polyform will not be generated.
- `BREAK`: the enumeration terminates immediately.

To ensure well-defined behavior, the function `f` may not accept the offspring of polyforms
that it rejects.

`checkbounded` refuses a rule set that admits arbitrarily large polyforms, see
[`isbounded`](@ref), by throwing an `ArgumentError`, and warns if that cannot be decided. It defaults
to `true` when neither `maxsize` nor `maxpolys` caps the enumeration.

All other keyword arguments are passed to the underlying `reversesearch` routine.

Returns a [`PolyformCount`](@ref) of what was enumerated: `n` polyforms, the largest of `maxsize`
particles, and the [`EnumerationStatus`](@ref) the run stopped with.
"""
function polyenum(f, rules::BindingRules; maxsize=Inf, maxpolys=Inf,
                  checkbounded::Bool=isinf(maxsize) && isinf(maxpolys), kwargs...)
    checkbounded && _checkbounded(rules)
    v₀ = Polyform(rules)
    BS = sitetype(rules)
    aux = PolyformAux{BS}(Set{NautyDiGraph}(), Tuple{BS,SpeciesSiteLoc,BondTwist}[])
    rsys = RSSystem(ls!, adj!, v₀; aux)

    frs = isnothing(f) ? nothing : (v, _) -> f(v, nparticles(v))
    result = reversesearch(frs, rsys; maxdepth=maxsize, maxverts=maxpolys+1, kwargs...)
    return PolyformCount(result.nvertices - 1, result.depth_reached, _status(result.status))
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

# Enumerate exactly at sizes 1, 2, 3, ... until number of polyforms reaches `budget`.
# Report the cumulative number of polyforms at each size, the largest polyform seen, and why the enumeration stopped.
function _count_upto_budget(rules::BindingRules; maxsize, budget)
    counts = Int[]
    largest = 0
    while length(counts) < maxsize
        res = polyenum(rules; maxsize=length(counts)+1, maxpolys=budget)
        res.status == MaxPolysReached && return counts, largest, res.status

        push!(counts, res.n)
        largest = res.maxsize
        res.status == Complete && return counts, largest, res.status
    end
    return counts, largest, MaxSizeReached
end

# Perform one unbiased sample of the number of polyforms of size ≤ maxsize. Polyforms up to `depth` are all
# visited and contribute the known count `n0`. Deeper ones are kept with probability `pkeep`, 
# so a polyform of size `n` is reached with probability `pkeep^(n-depth)` and is weighted accordingly. 
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

    res = polyenum(f, rules; maxsize, maxpolys=n0+maxsamples)
    return (; estimate=estimate[], largest=largest[], res.status)
end

# Optimize `pkeep` on throwaway pilot runs so the estimate more likely
# reaches the bottom of the search tree, without oversampling.
function _calibrate_pkeep(rules::BindingRules; pkeep, depth, n0, maxsize, maxsamples, rng, npilots=4, pmax=0.95)
    pkeep = clamp(pkeep, eps(float(pkeep)), pmax)
    for _ in 1:npilots
        pkeep >= pmax && return pmax
        pilot = _estimatecount(rules; pkeep, depth, n0, maxsize, maxsamples, rng)
        if pilot.status == MaxPolysReached # if we reached maxsamples, reduce pkeep
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
    polycount(rules::BindingRules; maxsize=Inf, budget=5000, kwargs...)

Estimate the number of polyforms allowed by `rules`, returning a [`PolyformCount`](@ref).

The count is exact whenever the enumeration can fit into `budget`, and is otherwise estimated by randomly
subsampling the reverse-search tree beyond the budget. Only polyforms of at most `maxsize` particles are counted.
With no `maxsize` the rule set has to be bounded, which is checked first, and the estimate then covers every
polyform there is.

- `maxsize=Inf`: only count polyforms with at most this many particles.
- `budget=5000`: how many polyforms to enumerate exactly before switching to estimation.
- `ntrials=5`: number of independent estimates to average.
- `pkeep=nothing`: the probability a polyform beyond the budget is kept. Set from the branching
  of the exact part and refined on pilot runs unless given.
- `maxsamples=10^6`: how many polyforms a single trial may generate beyond the exact budget.
- `rng=Random.default_rng()`: generator used for the subsampling.
- `checkbounded=isinf(maxsize)`: refuse a rule set that admits arbitrarily large polyforms, see
  [`isbounded`](@ref), rather than count until the budget runs out, and warn if that cannot be decided.

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
    budget=5_000,
    ntrials=5,
    pkeep=nothing,
    maxsamples=10^6,
    rng=Random.default_rng(),
    checkbounded::Bool=isinf(maxsize),
)
    checkbounded && _checkbounded(rules)
    budget = max(budget, nspecies(rules) + nbonds(rules) + 1) # ensure we visit at least all monomers and dimers
    counts, largest, status = _count_upto_budget(rules; maxsize, budget)

    if status != MaxPolysReached
        return PolyformCount(counts[end], largest, status)
    end

    # past the budget with no cap, the estimate terminates only on a bounded rule set, which the
    # check has established unless it was waived
    if isinf(maxsize) && !checkbounded
        throw(
            ArgumentError(
                "`rules` may allow infinitely many polyforms. Pass `maxsize` to count up to a given number of particles, or " *
                "`checkbounded=true` to establish that the rules are bounded.",
            ),
        )
    end

    depth, n0 = length(counts), counts[end]
    if isnothing(pkeep)
        counts_per_size = diff([0; counts])
        branching = counts_per_size[end] / counts_per_size[end - 1]
        # aim for 1.2 surviving offspring per kept polyform, which keeps the sample alive without
        # letting it grow; the pilots then adjust it
        pkeep = _calibrate_pkeep(rules; pkeep=1.2 / branching, depth, n0, maxsize, maxsamples, rng)
    elseif !(0 < pkeep <= 1)
        throw(ArgumentError("pkeep=$pkeep must be between 0 and 1."))
    end

    trials = [_estimatecount(rules; pkeep, depth, n0, maxsize, maxsamples, rng) for _ in 1:ntrials]
    largest = maximum(t.largest for t in trials)
    samples = [t.estimate for t in trials if t.status != MaxPolysReached] # dont count trials that ran out of nodes

    if isempty(samples)
        throw(
            ErrorException(
                "all $ntrials trials exceeded maxsamples=$maxsamples; raise " * "`maxsamples` or lower `maxsize`"
            ),
        )
    end

    return PolyformCount(samples, largest, isfinite(maxsize) ? MaxSizeReached : Complete)
end

# Throws if `rules` is unbounded, and warns if it cannot be shown to be bounded. Runs before any
# enumeration without a cap.
function _checkbounded(rules::BindingRules)
    bounded = isbounded(rules)
    bounded === true && return nothing
    if ismissing(bounded)
        @warn "Could not decide whether `rules` allows finitely many polyforms, so the enumeration may not " *
              "terminate. Pass `maxsize` or `maxpolys` to enumerate up to a bound or set `checkbounded=false to disable this warning.`."
        return nothing
    end
    throw(
        ArgumentError(
            "`rules` allows infinitely many polyforms and cannot be enumerated completely. " *
            "Pass `maxsize` or `maxpolys` to enumerate up to a bound, or set `checkbounded=false` to disable this error.",
        ),
    )
end
