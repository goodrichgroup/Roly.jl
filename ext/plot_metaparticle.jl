function plot_particlespecies!(ax, spcs::MetaParticleSpecies, pose::Pose=Pose{2}();
                               site_color=nothing, species_index=1, sys=nothing, kwargs...)

    n = nsites(spcs)

    original_sys = assemblysystem(spcs.polyform)
    active_site_indices = spcs.active_site_indices
    active_sitelocs = collect(Iterators.flatten(
        Roly.color2siteloc(original_sys, color(bindingsites(spcs.polyform, i))) for i in active_site_indices))

    if !isnothing(sys)
        species_index = findfirst(==(spcs), species(sys))
    end

    if isnothing(site_color)
        if !isnothing(sys)
            pal = species_palette(species_index, n)
            inert_sites = [i for i in 1:n if isinert(sys, (species_index, i))]
        else
            pal = species_palette(1, n)
            inert_sites = Int[]
        end

        site_color = function (spcs_idx, part_site_idx)
            siteloc = (spcs_idx, part_site_idx)
            siteloc ∉ active_sitelocs && return INERT_COLOR
            meta_site_idx = findfirst(==(siteloc), active_sitelocs)
            meta_site_idx ∈ inert_sites && return INERT_COLOR
            return pal[meta_site_idx]
        end
    end

    return plot_polyform!(ax, spcs.polyform, pose; site_color, kwargs...)
end
