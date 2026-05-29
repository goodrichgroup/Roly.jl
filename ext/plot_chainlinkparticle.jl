function plot_particlespecies!(ax, spcs::ChainLinkParticleSpecies{F}, pose::Pose=Pose{2,F}();
                               site_color=nothing, species_index=nothing, sys=nothing,
                               strokewidth=3, kwargs...) where {F}
    n = nsites(spcs)
    a = spcs.a

    if !isnothing(sys) && isnothing(species_index)
        species_index = findfirst(==(spcs), species(sys))
    else
        species_index = 1
    end
    pal = species_palette(species_index, n)

    if isnothing(site_color)
        if !isnothing(sys)
            inert_sites = [i for i in 1:n if isinert(sys, (species_index, i))]
            site_color = (_, site_idx) -> site_idx ∈ inert_sites ? INERT_COLOR : pal[site_idx]
        else
            site_color = (_, site_idx) -> pal[site_idx]
        end
    end

    colors = [site_color(species_index, i) for i in 1:n]
    return _draw_chainlink!(ax, pose.x[1], pose.x[2], rotation_angle(pose.psi);
                            a, color=colors, strokewidth, kwargs...)
end

function _semicircle(a::Real)
    r = a / 2
    return BezierPath([
        MoveTo(Point(0, 0)),
        LineTo(Point(0, r)),
        EllipticalArc(Point(0, 0), r, r, 0, π/2, 3π/2),
        ClosePath()])
end

function _draw_chainlink!(ax, x, y, ψ; a, kwargs...)
    rotation = ψ .+ [0, π]
    marker = _semicircle(a)
    return scatter!(ax, fill(x, 2), fill(y, 2);
                    marker, rotation, markersize=1, markerspace=:data, kwargs...)
end
