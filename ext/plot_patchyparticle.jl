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
    species_index=nothing,
    sys=nothing,
    site_radius=0.25spcs.r,
    strokewidth=2,
    kwargs...,
) where {F}
    n = nsites(spcs)
    r = spcs.r

    _, pal, colors = _resolve_colors(spcs, species_index, sys, site_color)
    cx, cy = pose.x[1], pose.x[2]

    scatter!(ax, [cx], [cy];
             marker=_disk_marker(r), markersize=1, markerspace=:data,
             color=pal[1], strokecolor=:black, strokewidth, kwargs...)

    xs = [(pose * bindingsites(spcs, i).pose.x)[1] for i in 1:n]
    ys = [(pose * bindingsites(spcs, i).pose.x)[2] for i in 1:n]
    scatter!(ax, xs, ys;
             marker=_disk_marker(site_radius), markersize=1, markerspace=:data,
             color=colors, strokecolor=:black, strokewidth=1, alpha=0.8)

    return nothing
end

function plot_particlespecies!(
    ax,
    spcs::PatchyParticleSpecies{3,F},
    pose::Pose=Pose{3,F}();
    site_color=nothing,
    species_index=nothing,
    sys=nothing,
    site_radius=0.25spcs.r,
    kwargs...,
) where {F}
    n = nsites(spcs)
    _, pal, colors = _resolve_colors(spcs, species_index, sys, site_color)

    mesh!(ax, Sphere(Point3f(pose.x), Float32(spcs.r)); color=pal[1], kwargs...)
    for i in 1:n
        x = Point3f(pose * bindingsites(spcs, i).pose.x)
        mesh!(ax, Sphere(x, Float32(site_radius)); color=colors[i])
    end
    return nothing
end
