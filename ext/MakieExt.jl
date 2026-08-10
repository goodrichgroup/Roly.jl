module MakieExt

using Roly
using Makie
import Roly: species, bindingrules, canonical_vertices, polyformplot, polyformplot!, render

@recipe PolyformPlot (poly, ) begin
    bindingrules = nothing
    pose = nothing
    Makie.mixin_generic_plot_attributes()...
end

function Makie.plot!(p::PolyformPlot{<:Tuple{<:Polyform}})
    plot_polyform!(p, p.poly[], p.pose[])
    return p
end
function Makie.plot!(p::PolyformPlot{<:Tuple{<:ParticleSpecies}})
    # With no pose given, let each species method supply its own identity-pose default.
    pose = p.pose[]
    isnothing(pose) ? plot_particlespecies!(p, p.poly[]) : plot_particlespecies!(p, p.poly[], pose)
    return p
end

"""
    render(p; hidedecorations=true, kwargs...)

Draw a `Polyform` or `ParticleSpecies` into a fresh figure, picking a 2D or 3D axis to match
its dimension.

3D output needs a backend with a depth buffer — GLMakie or WGLMakie. CairoMakie sorts
primitives instead of depth-testing them, so faces of a polyform will occlude each other
incorrectly.
"""
function render(p; hidedecorations=true, kwargs...)
    f = Figure()
    ax = dimension(p) == 3 ? Axis3(f[1, 1], aspect=:data) : Axis(f[1, 1], aspect=DataAspect())
    hidedecorations && hidedecorations!(ax)
    polyformplot!(ax, p; kwargs...)
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
    sys = bindingrules(poly)

    for part in poly.particles
        ps = species(sys, part.species_index)
        part_pose = isnothing(pose) ? part.pose : pose * part.pose
        plot_particlespecies!(ax, ps, part_pose; sys, kwargs...)
    end
    return
end

include("palette.jl")
include("plot_polygonparticle.jl")
include("plot_polyhedronparticle.jl")
include("plot_patchyparticle.jl")

end # module
