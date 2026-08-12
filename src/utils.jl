"""
    sat_overlap(axes, corners1, pose1, corners2, pose2, skin)

Return `true` unless some axis in `axes` separates the two convex bodies, i.e. `true` if they
overlap, with a `skin` of clearance counting as separated.

Two convex bodies are disjoint exactly when some axis exists on which their projections do not
overlap, and it is enough to test a finite candidate set — which is all that differs between
dimensions, so it is the caller's to supply. In 2D the edge normals of both polygons suffice;
in 3D it takes both solids' face normals *and* the cross products of their edge directions,
which catch the edge-on-edge configurations no face normal separates.

Axes need not be normalised or even nonzero: each is scaled here, and a degenerate one (from
parallel edges, say) carries no information and is skipped. Normalising is not cosmetic —
`skin` is a length, so comparing it against projections along an unnormalised axis would scale
the tolerance by that axis's magnitude.
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
    translate_overlap(normals, offsets, t, skin)

Return `true` unless a centrally symmetric convex body and an identically oriented copy of
itself, offset by `t`, are separated by a `skin` of clearance. `normals` and `offsets` are the
body's outward face normals and the distances from its centre to those faces; `t` is the offset
expressed in the body's own frame, so all three are in body coordinates and no rotation happens
here.

Exact, and O(#faces) where [`sat_overlap`](@ref) is O(#faces × #corners) with a quadratic
candidate set behind it. Two bodies `K` and `K + t` meet exactly when `t` lies in the difference
body `K ⊕ (−K)`, whose faces are those of `K` together with their opposites — and that is `2K`
when `K` is centrally symmetric, so the whole question is whether `t` clears each face plane at
twice its offset. Both conditions matter: the shapes must be the same, and the body must be
centrally symmetric, or the difference body is not `2K` and the offsets are wrong.
"""
function translate_overlap(normals, offsets, t, skin::Real)
    for (nrm, d) in zip(normals, offsets)
        abs(dot(nrm, t)) >= 2d - skin && return false
    end
    return true
end

"""
    edgenormals(corners, pose)

The outward-ish normals of a 2D polygon's edges, in world coordinates: candidate separating
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