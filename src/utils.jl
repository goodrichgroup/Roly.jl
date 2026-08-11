"""
    bfs!(neighbors, dist, queue, seeds; maxdepth=typemax(Int))

Breadth-first search from `seeds` over the adjacency defined by `neighbors(v)`, which must return
an iterable of the neighbors of `v`. The result is consumed before the next call, so `neighbors`
may fill and return a reusable buffer.

`dist` must hold -1 for every unvisited vertex on entry; any other value marks a vertex as
unvisitable. On return, reached vertices hold their distance from the seeds. Vertices at `maxdepth`
are reached but not expanded. `queue` is a scratch buffer of the same length as `dist`.
"""
function bfs!(neighbors::F, dist::AbstractVector{<:Integer}, queue::AbstractVector{<:Integer},
              seeds; maxdepth=typemax(Int)) where {F}
    qend = 1
    for s in seeds
        dist[s] = 0
        queue[qend] = s
        qend += 1
    end
    qstart = 1
    while qend > qstart
        v = queue[qstart]
        qstart += 1
        d = dist[v]
        d >= maxdepth && continue
        for w in neighbors(v)
            dist[w] == -1 || continue
            dist[w] = d + 1
            queue[qend] = w
            qend += 1
        end
    end
    return dist
end

"""
    is_cutset(g, vs[; target, dist, queue])

Return `true` if removing vertices `vs` from `g` disconnects the graph.

`target`, `dist`, and `queue` are buffer arrays of length `nv(g)`, which will be modified.
"""
function is_cutset(g::AbstractNautyGraph, vs::AbstractVector{<:Integer};
    target::AbstractVector{Bool}=zeros(Bool, nv(g)),
    dist::AbstractVector{<:Integer}=zeros(Int, nv(g)),
    queue::AbstractVector{<:Integer}=zeros(Int, nv(g)))

    fill!(target, false)
    fill!(dist, -1)
    for v in vs
        target .|= NautyGraphs.adjrow(g, v)
        dist[v] = -2   # removed vertices are unvisitable
    end
    target[vs] .= false

    count(target) < 2 && return false

    bfs!(dist, queue, (findfirst(target),)) do v
        outneighs = NautyGraphs.adjrow(g, v)
        (w for w in eachindex(outneighs) if outneighs[w])
    end

    # `vs` is a cutset iff some vertex adjacent to it was not reached from the others.
    for (w, t) in enumerate(target)
        t && dist[w] < 0 && return true
    end
    return false
end

function blockdiag!(g::AbstractNautyGraph, h::AbstractNautyGraph)
    ng, nh = nv(g), nv(h)
    add_vertices!(g, nh; vertex_labels=labels(h))
    for e in edges(h)
        add_edge!(g, e.src+ng, e.dst+ng)
    end
    return g
end
