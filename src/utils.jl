function are_cutvertices(g::AbstractNautyGraph, vs::AbstractVector{<:Integer})
    n = nv(g)
    return are_cutvertices!(g, vs, zeros(Bool, n), zeros(Bool, n), zeros(Bool, n), zeros(Cint, n))
end

function are_cutvertices!(g::AbstractNautyGraph, vs::AbstractVector{<:Integer},
    is_target::AbstractVector{Bool}, forbidden::AbstractVector{Bool},
    explored::AbstractVector{Bool}, queue::AbstractVector)

    fill!(is_target, false)
    fill!(forbidden, false)
    forbidden[vs] .= true
    for v in vs
        is_target .|= NautyGraphs.adjrow(g, v)
    end
    is_target[vs] .= false

    n_targets = count(is_target)
    n_targets < 2 && return false

    v0 = findfirst(is_target)
    fill!(explored, false)
    explored[v0] = true
    n_found = 1

    queue[1] = v0
    q_start = 1
    q_end = 2

    while q_end > q_start
        v = queue[q_start]
        q_start += 1

        outneighs = NautyGraphs.adjrow(g, v)
        for neigh in eachindex(outneighs)
            outneighs[neigh] || continue
            forbidden[neigh] && continue
            if !explored[neigh]
                explored[neigh] = true
                if is_target[neigh]
                    n_found += 1
                    n_found == n_targets && return false
                end
                queue[q_end] = neigh
                q_end += 1
            end
        end
    end

    return true
end

function take_nth(itr, n, default=nothing)
    for (i, val) in enumerate(itr)
        i == n && return val
    end
    return default
end


mutable struct PathIterator{G,T}
    g::G
    v::T
    w::T
end
PathIterator(g, v, w) = PathIterator(g, promote(v, w)...)
Base.IteratorSize(::PathIterator) = Base.SizeUnknown()

function Base.iterate(pathitr::PathIterator)
    path = zeros(Int, nv(pathitr.g))
    path[pathitr.v] = 1
    path = complete_path!(path, pathitr.g, pathitr.w)    
    return copy(path), path #TODO: optimize this
end

function Base.iterate(pathitr::PathIterator, state)
    path = state
    neighs = zeros(Int, nv(pathitr.g))

    path_length = maximum(path)

    u = pathitr.w
    for _ in 1:path_length-1
        k = path[u]
        path[u] = 0

        parent = findfirst(x->x==k-1, path)
        n_neighs = NautyGraphs.outneighbors!(neighs, pathitr.g, parent)
        idx = @views searchsortedfirst(neighs[1:n_neighs], u)

        if idx < n_neighs && path[neighs[idx + 1]] == 0
            path[neighs[idx + 1]] = k
            path = complete_path!(path, pathitr.g, pathitr.w, neighs)
            return copy(path), path
        else
            k -= 1
            u = parent
        end
    end

    return nothing
end

function complete_path!(path, g, w, neighs=nothing)
    ## Completes the path by depth first traversal
    ## Neighbors are picked in order of vertex number (NOT label)
    neighs = isnothing(neighs) ? zeros(Int, length(path)) : neighs

    v = argmax(path)
    v_last = 0
    k = sum(!iszero, path) + 1

    while path[w] == 0
        n_neighs = NautyGraphs.outneighbors!(neighs, g, v)
        n0 = let i=findfirst(x->x==v_last, @view neighs[1:n_neighs]) # use searchsortedfirst
            isnothing(i) ? 1 : i+1
        end
        v_last = 0

        success = false
        for neigh in @view neighs[n0:n_neighs] # Assume neighbors are sorted
            if path[neigh] == 0
                v = neigh
                path[v] = k
                k += 1
                success = true
                break
            end
        end

        if !success
            path[v] = 0
            k -= 1
            v_last = v
            v = argmax(path)
        end
    end

    return path
end


function blockdiag!(g::AbstractNautyGraph, h::AbstractNautyGraph)
    ng, nh = nv(g), nv(h)
    add_vertices!(g, nh; vertex_labels=labels(h))
    for e in edges(h)
        add_edge!(g, e.src+ng, e.dst+ng)
    end
    return g
end