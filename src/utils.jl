"""
    is_cutset(g, vs[; target, visited, queue])

Return `true` if removing vertices `vs` from `g` disconnects the graph.

`target`, `visited`, and `queue` are buffer arrays of length `nv(g)`, which will be modified.
"""
function is_cutset(g::AbstractNautyGraph, vs::AbstractVector{<:Integer};
    target::AbstractVector{Bool}=zeros(Bool, nv(g)),
    visited::AbstractVector{Bool}=zeros(Bool, nv(g)),
    queue::AbstractVector=zeros(Cint, nv(g)))

    fill!(target, false)
    fill!(visited, false)
    visited[vs] .= true
    for v in vs
        target .|= NautyGraphs.adjrow(g, v)
    end
    target[vs] .= false

    n_targets = count(target)
    n_targets < 2 && return false

    v0 = findfirst(target)
    visited[v0] = true
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
            visited[neigh] && continue
            visited[neigh] = true
            if target[neigh]
                n_found += 1
                n_found == n_targets && return false
            end
            queue[q_end] = neigh
            q_end += 1
        end
    end

    return true
end


function blockdiag!(g::AbstractNautyGraph, h::AbstractNautyGraph)
    ng, nh = nv(g), nv(h)
    add_vertices!(g, nh; vertex_labels=labels(h))
    for e in edges(h)
        add_edge!(g, e.src+ng, e.dst+ng)
    end
    return g
end