module MakieExt

using Roly
using Makie
import Roly: species, assemblysystem, canonical_vertices, particle, polyformplot, polyformplot!, render

@recipe PolyformPlot (poly, ) begin
    assemblysystem = nothing
    pose = nothing
    Makie.mixin_generic_plot_attributes()...
end

function Makie.plot!(p::PolyformPlot{<:Tuple{<:Polyform}})
    plot_polyform!(p, p.poly[], p.pose[])
    return p
end
function Makie.plot!(p::PolyformPlot{<:Tuple{<:ParticleSpecies}})
    plot_particlespecies!(p, p.poly[], p.pose[])
    return p
end

function render(p; hidedecorations=true, kwargs...)
    f = Figure()
    ax = Axis(f[1, 1], aspect=DataAspect())
    hidedecorations && hidedecorations!(ax)
    polyformplot!(ax, p)
    return f
end


"""
    plot_particlespecies!(ax, spcs::ParticleSpecies, pose::Pose; site_color, kwargs...)

Draw the particle species `spcs` at `pose` onto `ax`.

`site_color` is an optional callback `site_index -> color`
that controls per-site coloring. When omitted, a default palette is used.
"""
function plot_particlespecies!(ax, spcs::ParticleSpecies, args...; kwargs...)
    error("$(typeof(spcs)) does not implement `plot_particlespecies`!")
end

"""
    plot_polyform!(ax, poly::Polyform, pose=nothing; kwargs...)

Draw all particles of `poly` onto `ax`, coloring inert sites with `INERT_COLOR`.
"""
function plot_polyform!(ax, poly::Polyform, pose=nothing; kwargs...)
    sys = assemblysystem(poly)

    for part in poly.particles
        ps = species(sys, part.species_index)
        part_pose = isnothing(pose) ? part.pose : pose * part.pose
        plot_particlespecies!(ax, ps, part_pose; assemblysystem=sys, kwargs...)
    end
    return
end

include("palette.jl")
include("plot_polygonparticle.jl")

end # module
