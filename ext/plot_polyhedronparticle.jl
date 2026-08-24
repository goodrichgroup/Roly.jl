"""
    _faceloop(p::Polyhedron, i, pose)

The corners of face `i` of `p`, transformed by `pose`, in winding order.
"""
function _faceloop(p::Polyhedron, i::Integer, pose::Pose)
    return [Point3f(pose * corners(p)[v]) for v in facevertices(p, i)]
end

"""
    _facemesh!(pts, tris, ptcolors, p::Polyhedron, pose, colors, border, bordercolors)

Append the triangulation of `p` at `pose` to the given buffers: every face is drawn in its own
color, with a rim band of `bordercolors` covering the outer `border` fraction of the face.

Each face is a fan around its centroid out to an inner ring, plus a band of quads from that
ring to the true corners. The outer boundary is the face's real corners, so the solid stays
closed and adjacent faces meet exactly; the band only recolors the outer sliver.

The outline is part of the mesh rather than a separate `linesegments!` on purpose. A backend
that composites plot by plot — CairoMakie — draws a separate line plot over the whole mesh
regardless of depth, so outlines on the far side of a solid would show through. Inside the
mesh the band is sorted with the triangles it belongs to.
"""
function _facemesh!(pts, tris, ptcolors, p::Polyhedron, pose::Pose, colors, border::Real, bordercolors)
    for i in 1:nfaces(p)
        rim = _faceloop(p, i, pose)
        n = length(rim)
        c = Point3f(pose * facecentroid(p, i))
        inner = [c + (1 - border) * (r - c) for r in rim]

        # Centroid, then the inner ring in the face color; then a second copy of the inner
        # ring and the true corners in the border color. The ring is duplicated so the band
        # has a crisp edge instead of a gradient into the face color.
        base = length(pts)
        push!(pts, c)
        append!(pts, inner)
        append!(ptcolors, fill(colors[i], n + 1))
        append!(pts, inner)
        append!(pts, rim)
        append!(ptcolors, fill(bordercolors[i], 2n))

        fan(k) = base + 1 + k
        dup(k) = base + 1 + n + k
        out(k) = base + 1 + 2n + k
        for k in 1:n
            k2 = mod1(k + 1, n)
            push!(tris, (base + 1, fan(k), fan(k2)))
            push!(tris, (dup(k), out(k), out(k2)))
            push!(tris, (dup(k), out(k2), dup(k2)))
        end
    end
    return nothing
end

function particlemesh(
    spcs::PolyhedronParticleSpecies{F},
    pose::Pose=Pose{3,F}();
    site_color=nothing,
    speciesindex=nothing,
    rules=nothing,
    alpha=0.8,
    border=0.035,
) where {F}
    _, _, colors, _ = _resolve_colors(
        spcs, speciesindex, rules, site_color; bond_tint=FACE_TINT, inert_color=nothing
    )
    facecolors = [RGBAf(c, alpha) for c in colors]
    bordercolors = [RGBAf(_shade(c, 0.35), alpha) for c in colors]

    pts, tris, ptcolors = Point3f[], NTuple{3,Int}[], RGBAf[]
    _facemesh!(pts, tris, ptcolors, polyhedron(spcs), pose, facecolors, border, bordercolors)
    return pts, tris, ptcolors
end

function plot_particlespecies!(
    ax, spcs::PolyhedronParticleSpecies{F}, pose::Pose=Pose{3,F}(); shading=NoShading, kwargs...
) where {F}
    pts, tris, ptcolors = particlemesh(spcs, pose; kwargs...)
    # Flat shading, because the point of the pale tints is that a face's color identifies it.
    # Diffuse lighting would darken faces by their orientation and swamp that.
    mesh!(ax, pts, _facematrix(tris); color=ptcolors, shading)
    return nothing
end
