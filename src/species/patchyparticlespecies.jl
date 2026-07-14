"""
    PatchyParticleSpecies{D,F,B}

A `D`-dimensional sphere of radius `r` with binding sites on its surface (a disk in 2D).
"""
struct PatchyParticleSpecies{D,F,B<:BindingSite} <: ParticleSpecies{D,F,B}
    g::NautyDiGraph
    sites::Vector{B}
    r::F
    skin::F
end

"""
    PatchyParticleSpecies(g, r, poses; colors=eachindex(poses))

General constructor. `g` must have ≥ 3 vertices and contain a directed cycle; vertex labels
must be pre-set on `g`. Site `i` occupies graph vertex `i` (`vertices = i:i`). For truncated
sites (multiple graph vertices per site) build the graph and `BindingSite`s manually.
"""
function PatchyParticleSpecies(
    g::NautyDiGraph,
    r::Real,
    patch_positions,
    patch_twists=zeros(float(typeof(r)), length(patch_positions));
    colors=1:length(patch_positions)
)
    D = length(first(patch_positions))
    F = float(eltype(first(patch_positions)))
    r = F(r)

    tol = sqrt(eps(F)) * r
    sites = [
        BindingSite(normal_pose(patch_positions[i], patch_twists[i]), colors[i], i:i, tol, tol / r) for
        i in eachindex(patch_positions, patch_twists, colors)
    ]
    return PatchyParticleSpecies{D,F,eltype(sites)}(g, sites, r, tol)
end

function PatchyDisk(angles, r=1; colors=1:length(angles), labels=colors)
    F = float(eltype(angles))
    r = F(r)
    n = length(angles)
    tol = sqrt(eps(F)) * r
    positions = [SVector(r * cos(F(phi)), r * sin(F(phi))) for phi in angles]

    if n == 2
        g = NautyDiGraph(cycle_digraph(4); vertex_labels=[labels[1], labels[1], labels[2], labels[2]])
        sites = [BindingSite(normal_pose(positions[1], F(0)), colors[1], 1:2, tol, tol / r),
                 BindingSite(normal_pose(positions[2], F(0)), colors[2], 3:4, tol, tol / r)]
    else
        g = NautyDiGraph(cycle_digraph(n); vertex_labels=labels)
        sites = [BindingSite(normal_pose(positions[i], F(0)), colors[i], i:i, tol, tol / r) for i in 1:n]
    end
    return PatchyParticleSpecies{2,F,eltype(sites)}(g, sites, r, tol)
end

function Base.show(io::Core.IO, ps::PatchyParticleSpecies{D}) where {D}
    return print(io, "$(D)d PatchyParticleSpecies with $(nsites(ps)) sites")
end

Base.copy(ps::PatchyParticleSpecies) =
    typeof(ps)(copy(ps.g), copy(ps.sites), ps.r, ps.skin)

dimension(::PatchyParticleSpecies{D}) where {D} = D
graphrep(ps::PatchyParticleSpecies) = ps.g
nsites(ps::PatchyParticleSpecies) = length(ps.sites)
bindingsites(ps::PatchyParticleSpecies, i::Integer) = ps.sites[i]
isconvex(::PatchyParticleSpecies) = true

function setcolors!(ps::PatchyParticleSpecies, colors::AbstractVector{<:Integer})
    length(colors) != nsites(ps) && throw(ArgumentError("incorrect number of colors"))
    for k in eachindex(ps.sites)
        s = ps.sites[k]
        ps.sites[k] = BindingSite(s.pose, colors[k], s.vertices, s.touching_tolerance, s.alignment_tolerance)
    end
    return nothing
end

function could_contact(
    p1::SpeciesAndPose{<:PatchyParticleSpecies}, p2::SpeciesAndPose{<:PatchyParticleSpecies}; kwargs...
)
    return true # this check would be identical with overlap, so no need to do it twice
end

function overlap(
    p1::SpeciesAndPose{<:PatchyParticleSpecies}, p2::SpeciesAndPose{<:PatchyParticleSpecies}; kwargs...
)
    spcs1, pose1 = p1
    spcs2, pose2 = p2
    return norm(pose1.x - pose2.x) < (spcs1.r + spcs2.r) - (spcs1.skin + spcs2.skin)
end

_bounding_radius(ps::PatchyParticleSpecies) = ps.r