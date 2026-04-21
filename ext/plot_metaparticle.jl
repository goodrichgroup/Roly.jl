function plot_particlespecies!(ax, spcs::MetaParticleSpecies, pose::Pose=Pose{2}();
                               site_color=nothing, species_index=nothing, sys=nothing, kwargs...)

    n = nsites(spcs)

    original_sys = assemblysystem(spcs.polyform)
    active_sites_polyformframe = [bindingsites(spcs.polyform, i) for i in spcs.active_site_indices]
    metasite_index(site) = findfirst(s -> s.vertices == site.vertices, active_sites_polyformframe)

    if !isnothing(sys) && isnothing(species_index)
        species_index = findfirst(==(spcs), species(sys))
    else
        species_index = 1
    end
    pal = species_palette(species_index, n)

    if isnothing(site_color)
        if !isnothing(sys)
            inert_sites = [i for i in 1:n if isinert(sys, (species_index, i))]
        else
            inert_sites = Int[]
        end
        site_color = function (spcs_idx, site_idx)
            site_idx ∈ inert_sites && return INERT_COLOR
            return pal[site_idx]
        end
    end

    for part in spcs.polyform.particles
        ps = species(original_sys, part.species_index)
        part_pose = isnothing(pose) ? part.pose : pose * part.pose
        
        scolor = function (_, site_idx)
            mindex = metasite_index(bindingsites(part, original_sys, site_idx))
            isnothing(mindex) && return INERT_COLOR
            return site_color(species_index, mindex)
        end

        plot_particlespecies!(ax, ps, part_pose; sys, site_color=scolor, species_index, kwargs...)
    end

    return
end
