"""
    _shrunkface(p::Polyhedron, i, pose, shrink)

The corners of face `i`, transformed by `pose` and pulled towards the face centroid by the
fraction `shrink`. Shrinking leaves a visible seam between neighboring faces, which is what
keeps the per-site coloring readable once several particles are drawn side by side.
"""
function _shrunkface(p::Polyhedron, i::Integer, pose::Pose, shrink::Real)
    c = facecentroid(p, i)
    return [Point3f(pose * (c + (1 - shrink) * (corners(p)[v] - c))) for v in facevertices(p, i)]
end

"""
    _facemesh(p::Polyhedron, pose, shrink, colors)

Triangulate every face of `p` as a fan around its centroid, giving all points of a face that
face's color. Returns the point list, an `ntriangles × 3` index matrix, and the point colors.
"""
function _facemesh(p::Polyhedron, pose::Pose, shrink::Real, colors)
    pts = Point3f[]
    ptcolors = similar(colors, 0)
    tris = Matrix{Int}(undef, sum(i -> length(facevertices(p, i)), 1:nfaces(p)), 3)

    t = 0
    for i in 1:nfaces(p)
        rim = _shrunkface(p, i, pose, shrink)
        base = length(pts)
        push!(pts, Point3f(pose * facecentroid(p, i)))
        append!(pts, rim)
        append!(ptcolors, fill(colors[i], length(rim) + 1))
        for k in eachindex(rim)
            t += 1
            tris[t, :] .= (base + 1, base + 1 + k, base + 1 + mod1(k + 1, length(rim)))
        end
    end
    return pts, tris, ptcolors
end

function plot_particlespecies!(
    ax,
    spcs::PolyhedronParticleSpecies{F},
    pose::Pose=Pose{3,F}();
    site_color=nothing,
    species_index=nothing,
    sys=nothing,
    shrink=0.04,
    strokewidth=1.5,
    strokecolor=:black,
    kwargs...,
) where {F}
    _, _, colors = _resolve_colors(spcs, species_index, sys, site_color)
    p = shape(spcs)

    pts, tris, ptcolors = _facemesh(p, pose, shrink, colors)
    mesh!(ax, pts, tris; color=ptcolors, kwargs...)

    if strokewidth > 0
        segs = Point3f[]
        for i in 1:nfaces(p)
            rim = _shrunkface(p, i, pose, shrink)
            for k in eachindex(rim)
                push!(segs, rim[k], rim[mod1(k + 1, length(rim))])
            end
        end
        linesegments!(ax, segs; color=strokecolor, linewidth=strokewidth)
    end
    return nothing
end
