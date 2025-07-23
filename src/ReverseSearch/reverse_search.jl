@enum RejectValue NOREJECT = 0 REJECT = 1 BREAK = 2
@enum RSStatus COMPLETE = 0 MAXVERTREACHED = 1 MAXDEPTHREACHED = 2 BREAKTRIGGERED = 3

struct RSSystem{isinplace,LS,ADJ,COM}
    ls::LS              # local search, ls(v), returns v_prev
    adj::ADJ            # adjacency oracle, adj(v, j, aux) return v_next(j), Δj. May modify aux.
    compare::COM        # comparator between vertices v, v' (default Base.:(==))
    RSSystem{isinplace}(ls, adj, compare) where {isinplace} = 
        new{isinplace, typeof(ls), typeof(adj), typeof(compare)}(ls, adj, compare)
end
isinplace(::RSSystem{iip}) where {iip} = iip

function RSSystem(ls, adj, compare=Base.:(==))
    ls_iip = SciMLBase.isinplace(ls, 2, "ls")
    adj_iip = SciMLBase.isinplace(adj, 4, "adj")

    if ls_iip != adj_iip
        error("Local search and adjacency function have incompatible call signatures. The functions need to either both be in place, or both be out of place.")
    end
    return RSSystem{ls_iip}(ls, adj, compare)
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

function rs(f, rsys::RSSystem, state::RSState; fargs=())
    break_flag = false

    while true
        success = reverse_traverse!(state, rsys)
        if success
            reject_val = f(state.v, state.depth, fargs...)

            if reject_val == BREAK
                break_flag = true
                break
            end
            if reject_val == REJECT
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

function _rsworker(f, rsys::RSSystem, input_queue, work_tokens, break_flag; depth_per_task, verts_per_task, fargs=())
    hasf = !isnothing(f)

    function callback(v, task_depth, start_depth, task_nv, args...)
        # If another worker already broke, also break immedetely.
        break_flag[] && return BREAK

        total_depth = task_depth + start_depth

        reject_val = hasf ? f(v, total_depth, args...) : NOREJECT

        if reject_val == BREAK
            Threads.atomic_or!(break_flag, true)
        elseif reject_val == NOREJECT && (task_nv[] >= verts_per_task || task_depth == depth_per_task)
            reject_val = REJECT

            if isinplace(rsys)
                put!(input_queue, (copy(v), total_depth))
            else
                put!(input_queue, (v, total_depth))
            end
        end
        return reject_val
    end

    while true
        task_nv = Ref(1)

        input = take!(input_queue)
        isnothing(input) && break
        v, start_depth = input

        put!(work_tokens, true)

        state = RSState(v; depth=0) # TODO pull out of this loop, then copy to it \\ add kwargs
        rs(callback, rsys, state; fargs=(start_depth, task_nv, fargs...))

        take!(work_tokens)

        if break_flag[]
            break
        end
    end
    return
end

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

struct RSIterator{RSYS<:RSSystem,VTY,A}
    rsys::RSYS
    v₀::VTY
    cached::Bool
    aux::A
    maxdepth::Union{Int,Float64}
end
function RSIterator(ls, adj, v₀; compare=Base.:(==), cached=true, aux=nothing, maxdepth=Inf)
    rsys = RSSystem(ls, adj, compare)
    return RSIterator(rsys, v₀, cached, aux, maxdepth)
end
function RSIterator(rsys::RSSystem, v₀; cached=true, aux=nothing, maxdepth=Inf)
    return RSIterator(rsys, v₀, cached, aux, maxdepth)
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
    state = RSState(iter.v₀; cached=iter.cached, aux=iter.aux)
    return (copy(state.v), state.depth), state
end

function reversesearch(f, rsys::RSSystem, v₀; threaded=false, cached=true, aux=nothing, kwargs...)
    state = RSState(v₀; cached, aux)
    return _reversesearch(f, rsys, state, Val(threaded); kwargs...)
end
reversesearch(rsys::RSSystem, v₀; kwargs...) = reversesearch(nothing, rsys, v₀; kwargs...)


function _reversesearch(f, rsys::RSSystem, state::RSState, ::Val{threaded}; maxdepth=Inf, maxverts=Inf, fargs=(), kwargs...) where {threaded}
    hasf = !isnothing(f)

    maxdepth_flag = threaded ? Threads.Atomic{Bool}(false) : Ref(false)
    maxvert_flag = threaded ? Threads.Atomic{Bool}(false) : Ref(false)

    nv = threaded ? Threads.Atomic{Int}(1) : Ref(1)
    lowest_depth = threaded ? Threads.Atomic{Int}(1) : Ref(1)

    function callback(v, depth, args...)
        if threaded
            Threads.atomic_or!(maxvert_flag, nv[] >= maxverts)
        else
            maxvert_flag[] = maxvert_flag[] || nv[] >= maxverts
        end
       
        if !maxvert_flag[]
            reject_val = hasf ? f(v, depth, args...) : NOREJECT
        else
            reject_val = BREAK
        end

        if reject_val == NOREJECT
            if threaded
                Threads.atomic_add!(nv, 1)
                Threads.atomic_max!(lowest_depth, depth)
            else
                nv[] += 1
                lowest_depth[] = max(lowest_depth[], depth)
            end
            if depth == maxdepth
                reject_val = REJECT
                maxdepth_flag[] = true
            end
        end

        return reject_val
    end

    rs_fn = threaded ? prs : rs
    break_flag = rs_fn(callback, rsys, state; fargs, kwargs...)

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