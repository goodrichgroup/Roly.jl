mutable struct SphereParticleSpecies{D,F,P<:Pose{D,F},B<:BindingSite} <: ParticleSpecies{D,F,P}
    g::NautyDiGraph
    sites::Vector{B}
    r::F
end
function SphereParticleSpecies(r::F, site_poses::AbstractVector{P};
                                colors=1:length(site_poses),
                                labels=colors) where {D,F<:Real,P<:Pose{D,F}}
    n = length(site_poses)
    touching_tol = sqrt(eps(F)) * r
    alignment_tol = n > 0 ? touching_tol / maximum(norm(p.x) for p in site_poses) : touching_tol
    B = BindingSite{P,F}
    sites = B[BindingSite(site_poses[i], colors[i], i:i, touching_tol, alignment_tol) for i in 1:n]
    error("the graph encoding depends on the distribution of binding sites")
    g_base = n > 0 ? complete_digraph(n) : SimpleDiGraph(0)
    g = NautyDiGraph(g_base; vertex_labels=collect(labels))
    return SphereParticleSpecies{D,F,P,B}(g, sites, r)
end

Base.show(io::Core.IO, ps::SphereParticleSpecies{D}) where {D} =
    print(io, "$(D)d SphereParticleSpecies with $(nsites(ps)) sites")

Base.copy(p::SphereParticleSpecies) = typeof(p)(copy(p.g), copy(p.sites), p.r)

dimension(::SphereParticleSpecies{D}) where {D} = D
graphrep(p::SphereParticleSpecies) = p.g
nsites(p::SphereParticleSpecies) = length(p.sites)
bindingsites(p::SphereParticleSpecies, i::Integer) = p.sites[i]
symmetrynumber(::SphereParticleSpecies) = 1
isconvex(::SphereParticleSpecies) = true

can_skip_overlap_check(::SphereParticleSpecies, ::SphereParticleSpecies) = true

function setcolors!(p::SphereParticleSpecies, colors::AbstractVector{<:Integer})
    length(colors) != nsites(p) && throw(ArgumentError("incorrect number of colors"))
    for k in eachindex(p.sites)
        s = p.sites[k]
        p.sites[k] = BindingSite(s.pose, colors[k], s.vertices, s.touching_tolerance, s.alignment_tolerance)
    end
    return
end

function could_contact(p1::SpeciesAndPose{<:SphereParticleSpecies},
                       p2::SpeciesAndPose{<:SphereParticleSpecies}; kwargs...)
    return overlap(p1, p2; kwargs...)
end

function overlap(p1::SpeciesAndPose{<:SphereParticleSpecies},
                 p2::SpeciesAndPose{<:SphereParticleSpecies}; kwargs...)
    spcs1, pose1 = p1
    spcs2, pose2 = p2
    return norm(pose1.x - pose2.x) < spcs1.r + spcs2.r
end
