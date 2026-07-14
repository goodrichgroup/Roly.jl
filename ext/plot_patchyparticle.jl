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

    cx, cy = pose.x[1], pose.x[2]

    scatter!(ax, [cx], [cy];
             marker=_disk_marker(r), markersize=1, markerspace=:data,
             color=pal[1], strokecolor=:black, strokewidth, kwargs...)

    xs     = [(pose * bindingsites(spcs, i).pose.x)[1] for i in 1:n]
    ys     = [(pose * bindingsites(spcs, i).pose.x)[2] for i in 1:n]
    colors = [site_color(species_index, i) for i in 1:n]
    scatter!(ax, xs, ys;
             marker=_disk_marker(site_radius), markersize=1, markerspace=:data,
             color=colors, strokecolor=:black, strokewidth=1, alpha=0.8)

    return nothing
end
