function plot_particlespecies!(ax, spcs::PolygonParticleSpecies{F}, pose::Pose=Pose{2,F}();
                               sitecolor=nothing, speciesindex=nothing, rules=nothing,
                               cornerradius=nothing, thickness=nothing, strokewidth=3, alpha=1,
                               strokecolor=(:black, alpha), kwargs...) where {F}
    isnothing(cornerradius) && (cornerradius = 0.3spcs.rmin)
    isnothing(thickness) && (thickness = 0.12_extent(spcs.corners))
    thickness = min(thickness, spcs.rmin - cornerradius)     # the inner rounding has to fit

    _, _, colors, _ = _resolve_colors(spcs, speciesindex, rules, sitecolor)
    return _draw_polygon!(ax, pose.x[1], pose.x[2], rotation_angle(pose.psi); corners=spcs.corners,
                          cornerradius, thickness, color=colors, strokewidth, strokecolor, alpha, kwargs...)
end

# The larger side of the bounding box of `corners`.
_extent(corners) = maximum(d -> maximum(c[d] for c in corners) - minimum(c[d] for c in corners), 1:2)

# Half the angle the boundary turns through at each corner, which is π/n at every corner of a
# regular n-gon.
function _halfturns(corners)
    n = length(corners)
    return map(1:n) do i
        e1, e2 = corners[i] - corners[mod1(i - 1, n)], corners[mod1(i + 1, n)] - corners[i]
        return atan(abs(cross(e1, e2)), dot(e1, e2)) / 2
    end
end

# The band along the edge from corner `A` to corner `B`, with the inside of the polygon to the right
# of that direction, in the frame of the polygon. The band is drawn in the frame of the edge, with
# `A` at the origin and the inside below, and runs to the bisectors of its two corners, where the
# rounding of each corner turns by `ψA` and `ψB`.
function _edgesegment(A, B, ψA, ψB, cornerradius::Real, thickness::Real)
    r, t = cornerradius, thickness
    L = norm(B - A)
    θ = atan(B[2] - A[2], B[1] - A[1])
    s, c = sincos(θ)
    T(p) = Point(A[1] + c * p[1] - s * p[2], A[2] + s * p[1] + c * p[2])
    u(φ) = Point(cos(φ), sin(φ))
    o1, o2 = Point(r * tan(ψA), -r), Point(L - r * tan(ψB), -r)
    i1, i2 = Point((t + r) * tan(ψA), -t - r), Point(L - (t + r) * tan(ψB), -t - r)
    return BezierPath([
        MoveTo(T(o1 + r * u(π / 2 + ψA))),
        EllipticalArc(T(o1), r, r, 0, θ + π / 2 + ψA, θ + π / 2),
        LineTo(T(o2 + r * u(π / 2))),
        EllipticalArc(T(o2), r, r, 0, θ + π / 2, θ + π / 2 - ψB),
        LineTo(T(i2 + r * u(π / 2 - ψB))),
        EllipticalArc(T(i2), r, r, 0, θ + π / 2 - ψB, θ + π / 2),
        LineTo(T(i1 + r * u(π / 2))),
        EllipticalArc(T(i1), r, r, 0, θ + π / 2, θ + π / 2 + ψA),
        ClosePath(),
    ])
end

# One band per edge, in site order, each drawn so that the inside lies to its right.
function _draw_polygon!(ax, x, y, ψ; corners, cornerradius, thickness, kwargs...)
    n = length(corners)
    halfturns = _halfturns(corners)
    counterclockwise = sum(cross(corners[i], corners[mod1(i + 1, n)]) for i in 1:n) > 0
    markers = map(1:n) do i
        j = mod1(i + 1, n)
        return counterclockwise ?
               _edgesegment(corners[j], corners[i], halfturns[j], halfturns[i], cornerradius, thickness) :
               _edgesegment(corners[i], corners[j], halfturns[i], halfturns[j], cornerradius, thickness)
    end
    return scatter!(ax, fill(x, n), fill(y, n);
                    marker=markers, rotation=ψ, markersize=1, markerspace=:data, kwargs...)
end
