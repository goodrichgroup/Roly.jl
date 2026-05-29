struct ChainLinkParticleSpecies{F,B<:BindingSite} <: ParticleSpecies{2,F,Pose{2,F,Angle2d{F}}}
    g::NautyDiGraph
    sites::Vector{B}
    n::Int
    a::F
    skin::F
end
function ChainLinkParticleSpecies(n::Real, a::F=1.0; colors=1:2, labels=colors) where {F<:Real}
    if n <= 1 || isinf(n)
        n = 0
    end
    tol = sqrt(eps(F)) * a

    ϕ = n == 0 ? 0 : π / n

    back = BindingSite(Pose(SVector{2,F}(pol2cart(a/2, π - ϕ)), Angle2d{F}(π - ϕ)), colors[1], 1:1, tol, tol / a)
    front = BindingSite(Pose(SVector{2,F}(pol2cart(a/2, 0 + ϕ)), Angle2d{F}(0 + ϕ)), colors[2], 2:3, tol, tol / a)
    sites = [back, front]

    g = NautyDiGraph(cycle_digraph(3); vertex_labels=[labels; labels[2]])
    return ChainLinkParticleSpecies{F,eltype(sites)}(g, sites, n, a, tol)
end

function Base.show(io::Core.IO, ps::ChainLinkParticleSpecies)
    return print(io, "$(dimension(ps))d ChainLinkParticleSpecies with $(nsites(ps)) sites")
end

function Base.copy(clp::ChainLinkParticleSpecies)
    return ChainLinkParticleSpecies(copy(clp.g), copy(clp.sites), clp.n, clp.a, clp.skin)
end

dimension(::ChainLinkParticleSpecies) = 2
graphrep(p::ChainLinkParticleSpecies) = p.g
nsites(::ChainLinkParticleSpecies) = 2
bindingsites(p::ChainLinkParticleSpecies, i::Integer) = p.sites[i]
function setcolors!(p::ChainLinkParticleSpecies, colors::AbstractVector{<:Integer})
    length(colors) != nsites(p) && throw(ArgumentError("incorrect number of colors"))
    for k in eachindex(p.sites)
        s = p.sites[k]
        p.sites[k] = BindingSite(s.pose, colors[k], s.vertices, s.touching_tolerance, s.alignment_tolerance)
    end
    return nothing
end

isconvex(::ChainLinkParticleSpecies) = true

function could_contact(
    p1::SpeciesAndPose{<:ChainLinkParticleSpecies}, p2::SpeciesAndPose{<:ChainLinkParticleSpecies}; kwargs...
)
    return true
end
function overlap(
    p1::SpeciesAndPose{<:ChainLinkParticleSpecies}, p2::SpeciesAndPose{<:ChainLinkParticleSpecies}; kwargs...
)
    spcs1, pose1 = p1
    spcs2, pose2 = p2
    skin_sum = spcs1.skin + spcs2.skin
    return norm(pose1.x - pose2.x) < (spcs1.a + spcs2.a) / 2 - skin_sum
end
