# @enum RejectValue REJECT = 0 ACCEPT = 1 BREAK = 2
const ACCEPT = 1
const REJECT = 0
const BREAK = -1
@enum RSStatus COMPLETE = 0 MAXVERTREACHED = 1 MAXDEPTHREACHED = 2 BREAKTRIGGERED = 3


"""
    RSSystem(ls, adj, compare=Base.:(==))

An RSSystem defines a reverse-search enumeration procedure by specifying the local search function `ls(v)`, the
adjacency oracle `adj(v, j, aux)`, a comparator function (defaults to `Base.:(==)`), and a starting vertex `v₀`.
The enumeration can be carried out by calling `reversesearch(::RSSystem)`, or by iterating over an `RSIterator(::RSSystem)`.

The local search and adjacency functions are expected to adhere to the following interfaces:

- `u = ls(v)` maps an object `v` to its parent `u`, such that `adj(u, j, aux) == v` for some `j`.
    An in-place version of the form `u = ls!(u, v)` is also supported (this version must also return `u`).
- `u, Δj = adj(v, j, aux)` maps an object `v` onto its `j`th neighbor `u`, optionally making use and/or modifying the
    auxilary information stored in `aux`. In many applications, not all values for `j` will lead to a valid object, in which cases
    it is convenient to allow `adj` to internally increase `j` until the next valid neighbor is found. `adj` must therefore also
    return the index jump `Δj`; at the subsequent reverse-search step, the search for the next valid neighbor will thus begin at index `j + Δj`. 
    An in-place version of the form `u, Δj = adj!(u, v, j, aux)` is also supported (this version must also return `u, Δj`).

Note that `ls` and `adj` need to either both be in-place, or both be out-of-place.

In some enumeration problems, especially when dealing with isomorphism-free generation, it can be convenient to pass additional information to the 
adjacency oracle, for example to avoid generating isomorphic neighbors. If the adjacency oracle is defined with this in mind, `aux` can be used to pass
the initial value of this auxilary data. If `aux` is defined, then at every object `v`, a new copy of `aux` is created and passed to the oracle as 
`adj(v, 1, copy(aux))`. This copy can then be used and/or modified by `adj` while generating the neighbors of `v`.
"""
struct RSSystem{isinplace,LS,ADJ,COM,VTY,ATY}
    ls::LS              # local search, ls(v), returns v_prev
    adj::ADJ            # adjacency oracle, adj(v, j, aux) return v_next(j), Δj. May modify aux.
    compare::COM        # comparator between vertices v, v' (default Base.:(==))
    v₀::VTY             # starting vertex
    aux::ATY            # (optional) auxilary data
    function RSSystem{isinplace}(ls, adj, v₀; compare=Base.:(==), aux=nothing) where {isinplace}
        return new{isinplace, typeof(ls), typeof(adj), typeof(compare), typeof(v₀), typeof(aux)}(ls, adj, compare, v₀, aux)
    end
end
isinplace(::RSSystem{iip}) where {iip} = iip

function RSSystem(ls, adj, args...; kwargs...)
    ls_iip = SciMLBase.isinplace(ls, 2, "ls")
    adj_iip = SciMLBase.isinplace(adj, 4, "adj")

    if ls_iip != adj_iip
        error("Local search and adjacency function have incompatible call signatures. The functions need to either both be in place, or both be out of place.")
    end
    return RSSystem{ls_iip}(ls, adj, args...; kwargs...)
end

mutable struct RSState{VTY,NCT}
    v::VTY
    _temp1::Union{VTY,Nothing} # Only used for inplace assignments
    _temp2::Union{VTY,Nothing} # Only used for inplace assignments
    counter::NCT
    depth::Int
end
function RSState(v; depth=0, cached::Bool=true, aux=nothing)
    if isnothing(aux)
        counter = cached ? CachedNeighborCounter() : SimpleNeighborCounter()
    else
        counter = cached ? CachedAuxNeighborCounter(aux) : SimpleAuxNeighborCounter(aux)
    end
    return RSState(copy(v), copy(v), copy(v), counter, depth)
end
hasaux(state::RSState) = hasaux(state.counter)

function forward_traverse!(state::RSState, rsys::RSSystem{isinplace}) where {isinplace}
    state.depth == 0 && return false

    if isinplace
        prev = rsys.ls(state._temp1, state.v)
        popvertex!(state.counter, rsys, state.v, prev, state._temp2)
        copy!(state.v, prev)
    else
        prev = rsys.ls(state.v)
        popvertex!(state.counter, rsys, state.v, prev)
        state.v = prev
    end
    state.depth -= 1
    return true
end

function reverse_traverse!(state::RSState, rsys::RSSystem{isinplace}) where {isinplace}
    while true
        if isinplace
            next, Δj = rsys.adj(state._temp1, state.v, countervalue(state.counter), auxvalue(state.counter))
        else
            next, Δj = rsys.adj(state.v, countervalue(state.counter), auxvalue(state.counter))
        end
        isnothing(next) && return false
        increment!(state.counter, Δj)

        # TODO: let adj return missing if the neighbor is invalid, but more could be generated
        # then add this line
        # ismissing(next) && continue

        if isinplace
            rsys.ls(state._temp2, next)
            !rsys.compare(state._temp2, state.v) && continue
            copy!(state.v, next)
        else
            !rsys.compare(rsys.ls(next), state.v) && continue
            state.v = next
        end

        state.depth += 1
        pushvertex!(state.counter)
        return true
    end
end

abstract type AbstractNeighborCounter end
abstract type AbstractSimpleNeighborCounter <: AbstractNeighborCounter end
abstract type AbstractCachedNeighborCounter <: AbstractNeighborCounter end
auxvalue(::AbstractNeighborCounter) = nothing

mutable struct SimpleNeighborCounter <: AbstractSimpleNeighborCounter
    j::Int
end
SimpleNeighborCounter() = SimpleNeighborCounter(1)
increment!(counter::AbstractSimpleNeighborCounter, Δj) = counter.j += Δj
pushvertex!(counter::AbstractSimpleNeighborCounter) = counter.j = 1
function popvertex!(counter::AbstractSimpleNeighborCounter, rsys::RSSystem{isinplace}, v, prev, temp=nothing) where {isinplace}
    counter.j = 1

    while true
        if isinplace
            next, Δj = rsys.adj(temp, prev, countervalue(counter), auxvalue(counter))
        else
            next, Δj = rsys.adj(prev, countervalue(counter), auxvalue(counter))
        end
        counter.j += Δj
        rsys.compare(next, v) && break
    end
    return
end
countervalue(counter::AbstractSimpleNeighborCounter) = counter.j
hasaux(::SimpleNeighborCounter) = false

mutable struct SimpleAuxNeighborCounter{A} <: AbstractSimpleNeighborCounter
    j::Int
    aux::A
    const aux_init::A
end
SimpleAuxNeighborCounter(aux) = SimpleAuxNeighborCounter{typeof(aux)}(1, copy(aux), copy(aux))
pushvertex!(counter::SimpleAuxNeighborCounter) = (counter.j = 1; counter.aux = copy(counter.aux_init))
popvertex!(counter::SimpleAuxNeighborCounter, args...) = (invoke(popvertex!, Tuple{SimpleNeighborCounter, typeof.(args)...}, counter, args...); counter.aux = copy(counter.aux_init))
hasaux(::SimpleAuxNeighborCounter) = true
auxvalue(counter::SimpleAuxNeighborCounter) = counter.aux

struct CachedNeighborCounter <: AbstractCachedNeighborCounter
    js::Vector{Int}
end
CachedNeighborCounter() = CachedNeighborCounter([1])
increment!(counter::AbstractCachedNeighborCounter, Δj) = counter.js[end] += Δj
pushvertex!(counter::AbstractCachedNeighborCounter) = push!(counter.js, 1)
popvertex!(counter::AbstractCachedNeighborCounter, args...) = pop!(counter.js)
countervalue(counter::AbstractCachedNeighborCounter) = counter.js[end]
hasaux(::CachedNeighborCounter) = false

struct CachedAuxNeighborCounter{A} <: AbstractCachedNeighborCounter
    js::Vector{Int}
    aux::Vector{A}
    aux_init::A
end
CachedAuxNeighborCounter(aux) = CachedAuxNeighborCounter{typeof(aux)}([1], [copy(aux)], copy(aux))
pushvertex!(counter::CachedAuxNeighborCounter) = (push!(counter.js, 1); push!(counter.aux, copy(counter.aux_init)))
popvertex!(counter::CachedAuxNeighborCounter, args...) = (pop!(counter.js); pop!(counter.aux))
auxvalue(counter::CachedAuxNeighborCounter) = counter.aux[end]
hasaux(::CachedAuxNeighborCounter) = true


"""
    rs(f, rsys::RSSystem, state::RSState; fargs=())

Low-level reverse-search function that should rarely be called directly.
See `reversesearch` or `RSIterator` for user-friendly alternatives.
"""
function rs(f, rsys::RSSystem, state::RSState; fargs=())
    break_flag = false

    while true
        success = reverse_traverse!(state, rsys)
        if success
            signal = f(state.v, state.depth, fargs...)

            if signal == BREAK
                break_flag = true
                break
            end
            if signal == REJECT
                forward_traverse!(state, rsys)
                continue
            end
        else
            success = forward_traverse!(state, rsys)
            !success && break
        end
    end
    return break_flag
end

"""
    prs(f, rsys::RSSystem, state::RSState; depth_per_task, verts_per_task, nthreads=Threads.nthreads(), fargs=())

Low-level, parallel implementation of reverse-search. This function should rarely be called directly.
See `reversesearch` or `RSIterator` for user-friendly alternatives.
"""
function prs(f, rsys::RSSystem, state::RSState; depth_per_task, verts_per_task, nthreads=Threads.nthreads(), fargs=())
    input_queue = Channel{Union{Nothing,Tuple{typeof(state.v),Int}}}(Inf)

    nworkers = min(Threads.nthreads(), nthreads) - 1
    work_tokens = Channel{Bool}(nworkers)
    break_flag = Threads.Atomic{Bool}(false)

    put!(input_queue, (copy(state.v), state.depth))

    tasks = [@spawn _rsworker(f, rsys, input_queue, work_tokens, break_flag; depth_per_task, verts_per_task, fargs) for _ in 1:nworkers]

    while true
        sleep(0.01)

        if break_flag[] || (isempty(work_tokens) && isempty(input_queue))
            # Terminate workers
            for _ in 1:nworkers
                put!(input_queue, nothing)
            end
            break
        end
    end

    foreach(wait, tasks)
    return break_flag[] # TODO: make sure this always returns the same value as the corresponding rs() call
end

function _rsworker(f, rsys::RSSystem, input_queue, work_tokens, break_flag; depth_per_task, verts_per_task, fargs=())
    hasf = !isnothing(f)

    function fwrap(v, task_depth, start_depth, task_nv, args...)
        # If another worker already broke, also break immedetely.
        break_flag[] && return BREAK

        total_depth = task_depth + start_depth

        signal = hasf ? f(v, total_depth, args...) : ACCEPT

        if signal == BREAK
            Threads.atomic_or!(break_flag, true)
        elseif signal == ACCEPT && (task_nv[] >= verts_per_task || task_depth == depth_per_task)
            signal = REJECT

            if isinplace(rsys)
                put!(input_queue, (copy(v), total_depth))
            else
                put!(input_queue, (v, total_depth))
            end
        end
        return signal
    end

    while true
        task_nv = Ref(1)

        input = take!(input_queue)
        isnothing(input) && break
        v, start_depth = input

        put!(work_tokens, true)

        state = RSState(v; depth=0) # TODO pull out of this loop, then copy to it \\ add kwargs
        rs(fwrap, rsys, state; fargs=(start_depth, task_nv, fargs...))

        take!(work_tokens)

        if break_flag[]
            break
        end
    end
    return
end

"""
    RSIterator(rsys::RSSystem; cached=true, maxdepth=Inf)

Create an iterable from the RSSystem `rsys` that makes it convienent to iterate
over the objects generated by reverse-search, e.g. via

```
for (v, depth) in RSIterator(rsys)
    # do something with v and depth
end
```

The iterator will enumerate all objects up to a depth of `maxdepth`. For more fine-grained
control over the enumeration process, use `reversesearch`.
"""
struct RSIterator{RSYS<:RSSystem}
    rsys::RSYS
    cached::Bool
    maxdepth::Union{Int,Float64}
function RSIterator(rsys::RSSystem; cached=true, maxdepth=Inf)
        return RSIterator{typeof(rsys)}(rsys, cached, maxdepth)
    end 
end

function Base.iterate(iter::RSIterator, state::RSState)
    if state.depth == iter.maxdepth
        forward_traverse!(state, iter.rsys)
    end
    not_finished = rs((_...)->ReverseSearch.BREAK, iter.rsys, state)
    if not_finished 
        return (copy(state.v), state.depth), state
    else
        return nothing
    end
end
function Base.iterate(iter::RSIterator)
    state = RSState(iter.rsys.v₀; cached=iter.cached, aux=iter.rsys.aux)
    return (copy(state.v), state.depth), state
end


"""
    reversesearch([f], rsys::RSSystem; threaded=false, cached=true, maxdepth=Inf, maxverts=Inf, fargs=(), kwargs...)

Perform reverse-search enumeration using the adjacency oracle, local search, comparator, and starting vertex defined in `rsys`.
During the enumeration, evaluate `f(v, depth)` on each object `v` generated at a certain `depth`. Stop the enumeration if a depth of 
`maxdepth` is reached, if `maxverts` objects have been generated, or if `f(v, depth)` returns the `BREAK` signal (see below).

If `threaded=true`, the enumeration is performed in parallel and the following additional keyword arguments become available:

- `nthreads`: the number of threads to use (defaults to `Threads.nthreads()`).
- `depth_per_task`: the maximal depth a single task will explore before terminating.
- `verts_per_task`: the maximal number of objects a single task will generate before terminating.

The optimal values for `depth_per_task` and `verts_per_task` are highly problem-specific, there are no default values and some tuning is usually required 
to achieve good performance.

The `cached` keyword argument determines whether information along the current branch in the search tree should be cached, or if it needs to be 
regenerated at each forward traverse. This should usually be left as `true`, unless you are dealing with large enumerations that reach very low depths or 
run into memory issues.

The optional function `f` can be used to both process the generated objects and to steer the enumeration procedure.
`f(v, depth, args...)` must take as inputs an object `v`, the `depth` at which `v` was found, and any number of optional arguments, which will be passed 
through via the `fargs` keyword argument. `f` must return one of three signals:

- `ACCEPT` (or `true`): reverse-search continues as normal.
- `REJECT` (or `false`): the children of the current object will not be generated and the enumeration continues from the parent of the current object.
- `BREAK`: the enumeration terminates immediately.

!!! warning "Warning"

    If `threaded=true`, `f` will be called from different threads. It is your responsibility to ensure that `f` is thread-safe.

The return value contains the final status of the enumeration, the number of generated vertices, and the lowest depth reached.
"""
function reversesearch(f, rsys::RSSystem; threaded=false, cached=true, kwargs...)
    state = RSState(rsys.v₀; cached, aux=rsys.aux)
    return _reversesearch(f, rsys, state, Val(threaded); kwargs...)
end
reversesearch(rsys::RSSystem; kwargs...) = reversesearch(nothing, rsys; kwargs...)

function _reversesearch(f, rsys::RSSystem, state::RSState, ::Val{threaded}; maxdepth=Inf, maxverts=Inf, fargs=(), kwargs...) where {threaded}
    hasf = !isnothing(f)

    maxdepth_flag = threaded ? Threads.Atomic{Bool}(false) : Ref(false)
    maxvert_flag = threaded ? Threads.Atomic{Bool}(false) : Ref(false)

    nv = threaded ? Threads.Atomic{Int}(1) : Ref(1)
    lowest_depth = threaded ? Threads.Atomic{Int}(1) : Ref(1)

    function fwrap(v, depth, args...)
        if threaded
            Threads.atomic_or!(maxvert_flag, nv[] >= maxverts)
        else
            maxvert_flag[] = maxvert_flag[] || nv[] >= maxverts
        end
       
        if !maxvert_flag[]
            signal = hasf ? f(v, depth, args...) : ACCEPT
        else
            signal = BREAK
        end

        if signal == ACCEPT
            if threaded
                Threads.atomic_add!(nv, 1)
                Threads.atomic_max!(lowest_depth, depth)
            else
                nv[] += 1
                lowest_depth[] = max(lowest_depth[], depth)
            end
            if depth == maxdepth
                signal = REJECT
                maxdepth_flag[] = true
            end
        end

        return signal
    end

    rs_fn = threaded ? prs : rs
    break_flag = rs_fn(fwrap, rsys, state; fargs, kwargs...)

    if maxvert_flag[]
        result = MAXVERTREACHED 
    elseif break_flag
        result = BREAKTRIGGERED
    elseif maxdepth_flag[]
        result = MAXDEPTHREACHED
    else
        result = COMPLETE
    end

    return result, nv[], lowest_depth[]
end