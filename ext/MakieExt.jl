module MakieExt

using Roly
using Makie
using LinearAlgebra: dot, normalize
import Roly: species, bindingrules, polyformplot, polyformplot!, render
import Roly: corners, faces, facevertices, nfaces, facecentroid, facenormal, polyhedron, bounding_radius

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
    # With no pose given, let each species method supply its own identity-pose default. The
    # rules are optional too, and decide which sites are drawn as bonding.
    pose = p.pose[]
    args = isnothing(pose) ? (p.poly[],) : (p.poly[], pose)
    plot_particlespecies!(p, args...; rules=p.bindingrules[])
    return p
end

# Documented on the stub in `Roly`, so that Documenter finds it without loading this extension.
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

`site_color` is an optional callback `(speciesindex, site_index) -> color` that controls
per-site coloring. When omitted, the species' palette is used, with sites that no bond can
use greyed out.
"""
function plot_particlespecies!(ax, spcs::ParticleSpecies, args...; kwargs...)
    error("$(typeof(spcs)) does not implement `plot_particlespecies`!")
end

"""
    particlemesh(spcs::ParticleSpecies, pose; kwargs...)

Return `(points, triangles, colors)` for `spcs` drawn at `pose`, or `nothing` for a species
with no triangle-mesh form.

This exists so that a whole polyform can go into a *single* `mesh!` call. A backend that
composites plot by plot rather than depth-testing — CairoMakie — draws whole plots in
insertion order, so one plot per particle means particles occlude each other by construction
order rather than by depth. Within one mesh CairoMakie does sort the triangles, so merging
fixes it. GLMakie has a depth buffer and is unaffected either way.
"""
particlemesh(::ParticleSpecies, args...; kwargs...) = nothing

# Triangle triples as the `ntriangles × 3` index matrix `mesh!` wants.
_facematrix(tris::Vector{NTuple{3,Int}}) = [tris[i][j] for i in eachindex(tris), j in 1:3]

"""
    plot_polyform!(ax, poly::Polyform, pose=nothing; kwargs...)

Draw all particles of `poly` onto `ax`, coloring inert sites with `INERT_COLOR`.

Species that provide a [`particlemesh`](@ref) are merged into one mesh so they depth-sort
against each other; the rest are drawn one plot per particle.
"""
function plot_polyform!(ax, poly::Polyform, pose=nothing; kwargs...)
    rules = bindingrules(poly)
    pts, tris, cols = Point3f[], NTuple{3,Int}[], RGBAf[]

    for part in poly.particles
        ps = species(rules, part.speciesindex)
        part_pose = isnothing(pose) ? part.pose : pose * part.pose
        geom = particlemesh(ps, part_pose; rules, kwargs...)
        if isnothing(geom)
            plot_particlespecies!(ax, ps, part_pose; rules, kwargs...)
            continue
        end
        p, t, c = geom
        offset = length(pts)
        append!(pts, p)
        append!(cols, c)
        for (a, b, d) in t
            push!(tris, (a + offset, b + offset, d + offset))
        end
    end

    isempty(tris) || mesh!(ax, pts, _facematrix(tris); color=cols, shading=NoShading)
    return
end

include("palette.jl")
include("plot_polygonparticle.jl")
include("plot_polyhedronparticle.jl")
include("plot_patchyparticle.jl")

end # module
