function _disk_marker(r::Real)
    return BezierPath([
        MoveTo(Point2f(r, 0)),
        EllipticalArc(Point2f(0, 0), r, r, 0, 0, 2π),
        ClosePath(),
    ])
end

function plot_particlespecies!(
    ax,
    spcs::PatchyParticleSpecies{2,F},
    pose::Pose=Pose{2,F}();
    site_color=nothing,
    speciesindex=nothing,
    sys=nothing,
    site_radius=0.25spcs.r,
    strokewidth=2,
    kwargs...,
) where {F}
    n = nsites(spcs)
    r = spcs.r

    si, _, colors, _ = _resolve_colors(spcs, speciesindex, sys, site_color)
    cx, cy = pose.x[1], pose.x[2]

    # The body is a much fainter wash of the species color, so the patches read on top of it.
    scatter!(ax, [cx], [cy];
             marker=_disk_marker(r), markersize=1, markerspace=:data,
             color=_tint(species_basecolor(si), PATCHY_BODY_TINT),
             strokecolor=:black, strokewidth, kwargs...)

    xs = [(pose * bindingsites(spcs, i).pose.x)[1] for i in 1:n]
    ys = [(pose * bindingsites(spcs, i).pose.x)[2] for i in 1:n]
    scatter!(ax, xs, ys;
             marker=_disk_marker(site_radius), markersize=1, markerspace=:data,
             color=colors, strokecolor=:black, strokewidth=1, alpha=0.8)

    return nothing
end

"""
    _spheremesh!(pts, tris, cols, center, r, color; nlat, nlon)

Append a UV-tessellated sphere to the given mesh buffers, with `color` shaded per vertex by
[`_lambert`](@ref). The pole rings are degenerate, which costs a few zero-area triangles and
saves special-casing them.
"""
function _spheremesh!(pts, tris, cols, center, r::Real, color; nlat::Int=18, nlon::Int=28)
    base = length(pts)
    for i in 0:nlat, j in 0:(nlon - 1)
        θ, φ = π * i / nlat, 2π * j / nlon
        nrm = Vec3f(sin(θ) * cos(φ), sin(θ) * sin(φ), cos(θ))
        push!(pts, Point3f(center + r * nrm))
        push!(cols, _lambert(color, nrm))
    end

    idx(i, j) = base + i * nlon + mod(j, nlon) + 1
    for i in 0:(nlat - 1), j in 0:(nlon - 1)
        push!(tris, (idx(i, j), idx(i + 1, j), idx(i + 1, j + 1)))
        push!(tris, (idx(i, j), idx(i + 1, j + 1), idx(i, j + 1)))
    end
    return nothing
end

function particlemesh(
    spcs::PatchyParticleSpecies{3,F},
    pose::Pose=Pose{3,F}();
    site_color=nothing,
    speciesindex=nothing,
    sys=nothing,
    site_radius=0.25spcs.r,
    patch_alpha=0.85,
) where {F}
    si, _, colors, _ = _resolve_colors(spcs, speciesindex, sys, site_color)
    pts, tris, cols = Point3f[], NTuple{3,Int}[], RGBAf[]

    # The body stays opaque — it is a solid particle, and a see-through sphere shows its own
    # far side as a bright crescent. It is instead a much fainter wash of the species color,
    # so the patches read on top of it.
    body = RGBAf(_tint(species_basecolor(si), PATCHY_BODY_TINT), 1)
    _spheremesh!(pts, tris, cols, Point3f(pose.x), Float32(spcs.r), body)
    for i in 1:nsites(spcs)
        x = Point3f(pose * bindingsites(spcs, i).pose.x)
        _spheremesh!(
            pts, tris, cols, x, Float32(site_radius), RGBAf(colors[i], patch_alpha);
            nlat=12, nlon=18,
        )
    end
    return pts, tris, cols
end

function plot_particlespecies!(
    ax, spcs::PatchyParticleSpecies{3,F}, pose::Pose=Pose{3,F}(); shading=NoShading, kwargs...
) where {F}
    pts, tris, cols = particlemesh(spcs, pose; kwargs...)
    mesh!(ax, pts, _facematrix(tris); color=cols, shading)
    return nothing
end
