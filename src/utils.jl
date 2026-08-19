"""
    sat_overlap(axes, corners1, pose1, corners2, pose2, skin)

Return `true` unless some axis in `axes` separates the two convex bodies, i.e. `true` if they
overlap, with a `skin` of clearance counting as separated.

Two convex bodies are disjoint exactly when some axis exists on which their projections do not
overlap, and it is enough to test a finite candidate set. In 2D the edge normals of both polygons suffice.
In 3D it takes both solids' face normals and the cross products of their edge directions,
which catch the edge-on-edge configurations no face normal separates.

Axes need not be normalized or even nonzero: each is scaled here, and a degenerate one (from
parallel edges, say) carries no information and is skipped.
"""
function sat_overlap(axes, corners1, pose1, corners2, pose2, skin::Real)
    for axis in axes
        n2 = dot(axis, axis)
        n2 < eps(typeof(n2)) && continue
        a = axis / sqrt(n2)
        lo1, hi1 = extrema(dot(a, pose1 * c) for c in corners1)
        lo2, hi2 = extrema(dot(a, pose2 * c) for c in corners2)
        (hi2 < lo1 + skin || hi1 < lo2 + skin) && return false
    end
    return true
end

"""
    edgenormals(corners, pose)

The normals of a 2D polygon's edges, in world coordinates: candidate separating
axes for [`sat_overlap`](@ref). Only the direction matters, so the sign is not fixed.
"""
edgenormals(corners, pose) = (
    let e = pose.psi * (corners[mod1(i + 1, length(corners))] - corners[i])
        SVector(-e[2], e[1])
    end
    for i in eachindex(corners)
)

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