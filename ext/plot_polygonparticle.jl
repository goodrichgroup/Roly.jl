function plot_particlespecies!(ax, spcs::PolygonParticleSpecies{F}, pose::Pose=Pose{2,F}();
                               site_color=nothing, species_index=nothing, sys=nothing,
                               cornerradius=nothing, strokewidth=3, kwargs...) where {F}
    n = nsites(spcs)
    a = 2spcs.rmin * tan(π / n)
    isnothing(cornerradius) && (cornerradius = 0.2a)

    _, _, colors = _resolve_colors(spcs, species_index, sys, site_color)
    return _draw_ngon!(ax, pose.x[1], pose.x[2], rotation_angle(pose.psi);
                       n, a, cornerradius, color=colors, strokewidth, kwargs...)
end

function _ngon_segment(n::Integer, a::Real, cornerradius::Real)
    ψ = π / n
    R = a/2 * cot(ψ)

    r = cornerradius
    x = r * tan(ψ)
    y = a - 2x

    start = Point(r * (tan(ψ) - sin(ψ)), -r * (1 - cos(ψ)))
    com = Point(a / 2, -R)
    return BezierPath([
        MoveTo(start - com),
        EllipticalArc(Point(x, -r) - com, r, r, 0, π/2 + ψ, π/2),
        LineTo(Point(x+y, 0) - com),
        EllipticalArc(Point(x+y, -r) - com, r, r, 0, π/2, π/2 - ψ),
        LineTo(Point(x+y, -r) - com),
        LineTo(Point(x, -r) - com),
        ClosePath()])
end

function _draw_ngon!(ax, x, y, ψ; n, a, cornerradius, kwargs...)
    rotation = [ψ + π - 2π/n * i for i in 0:n-1]
    marker = _ngon_segment(n, a, cornerradius)
    return scatter!(ax, fill(x, n), fill(y, n);
                    marker, rotation, markersize=1, markerspace=:data, kwargs...)
end
