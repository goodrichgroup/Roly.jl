function plot_particlespecies!(ax, spcs::PolygonParticleSpecies{F}, pose::Pose=Pose{2,F}();
                               sitecolor=nothing, speciesindex=nothing, rules=nothing,
                               cornerradius=nothing, thickness=nothing, strokewidth=3, alpha=1,
                               strokecolor=(:black, alpha), kwargs...) where {F}
    n = nsites(spcs)
    a = 2spcs.rmin * tan(π / n)
    isnothing(cornerradius) && (cornerradius = 0.3spcs.rmin)
    isnothing(thickness) && (thickness = 0.12_ngon_extent(n, spcs.rmin))
    thickness = min(thickness, spcs.rmin - cornerradius)     # the inner rounding has to fit

    _, _, colors, _ = _resolve_colors(spcs, speciesindex, rules, sitecolor)
    return _draw_ngon!(ax, pose.x[1], pose.x[2], rotation_angle(pose.psi);
                       n, a, cornerradius, thickness, color=colors, strokewidth, strokecolor, alpha, kwargs...)
end

function _ngon_extent(n::Integer, rmin::Real)
    R = rmin / cos(π / n)
    θs = [π - 2π * i / n + s * π / n for i in 0:(n - 1) for s in (-1, 1)]
    xs, ys = R .* cos.(θs), R .* sin.(θs)
    return max(maximum(xs) - minimum(xs), maximum(ys) - minimum(ys))
end

function _ngon_segment(n::Integer, a::Real, cornerradius::Real, thickness::Real)
    ψ = π / n
    R = a / 2 * cot(ψ)
    r, t = cornerradius, thickness
    com = Point(a / 2, -R)
    u(θ) = Point(cos(θ), sin(θ))
    o1, o2 = Point(r * tan(ψ), -r), Point(a - r * tan(ψ), -r)
    i1, i2 = Point((t + r) * tan(ψ), -t - r), Point(a - (t + r) * tan(ψ), -t - r)
    return BezierPath([
        MoveTo(o1 + r * u(π / 2 + ψ) - com),
        EllipticalArc(o1 - com, r, r, 0, π / 2 + ψ, π / 2),
        LineTo(o2 + r * u(π / 2) - com),
        EllipticalArc(o2 - com, r, r, 0, π / 2, π / 2 - ψ),
        LineTo(i2 + r * u(π / 2 - ψ) - com),
        EllipticalArc(i2 - com, r, r, 0, π / 2 - ψ, π / 2),
        LineTo(i1 + r * u(π / 2) - com),
        EllipticalArc(i1 - com, r, r, 0, π / 2, π / 2 + ψ),
        ClosePath(),
    ])
end

function _draw_ngon!(ax, x, y, ψ; n, a, cornerradius, thickness, kwargs...)
    rotation = [ψ + π - 2π / n * i for i in 0:(n - 1)]
    marker = _ngon_segment(n, a, cornerradius, thickness)
    return scatter!(ax, fill(x, n), fill(y, n);
                    marker, rotation, markersize=1, markerspace=:data, kwargs...)
end
