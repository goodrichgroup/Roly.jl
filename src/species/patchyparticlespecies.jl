"""
    PatchyParticleSpecies{D,F,B}

A `D`-dimensional sphere of radius `r` with binding sites on its surface (a disk in 2D).
"""
struct PatchyParticleSpecies{D,F,B<:BindingSite} <: ParticleSpecies{D,B}
    g::NautyDiGraph
    sites::Vector{B}
    r::F
    skin::F
end

"""
    PatchyParticleSpecies(g, r, poses; colors=eachindex(poses))

General constructor. `g` must have one vertex per patch; vertex labels must be pre-set on
`g`, and its automorphism group must match the symmetry of `patch_positions` — the graph is
yours to build, but [`_check_encoding`](@ref) still verifies it describes the arrangement. Site `i` occupies graph vertex `i` (`vertices = i:i`). For truncated
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
        BindingSite(normal_pose(patch_positions[i], patch_twists[i]), colors[i], i:i, tol, tol / r, 1) for
        i in eachindex(patch_positions, patch_twists, colors)
    ]
    return _check_encoding(PatchyParticleSpecies{D,F,eltype(sites)}(g, sites, r, tol))
end

"""
    PatchyDisk(angles, r=1; colors=1:length(angles), labels=colors)

A 2D disk of radius `r` with one binding site placed at each angle in `angles` (radians,
measured counterclockwise from the +x axis). `colors` and `labels` optionally override the
default per-site color and vertex label assignments.
"""
function PatchyDisk(angles, r=1; colors=1:length(angles), labels=colors)
    F = float(eltype(angles))
    r = F(r)
    n = length(angles)
    tol = sqrt(eps(F)) * r
    positions = [SVector(r * cos(F(phi)), r * sin(F(phi))) for phi in angles]

    g, ranges = cycleencoding(n; labels)
    sites = [BindingSite(normal_pose(positions[i], F(0)), colors[i], ranges[i], tol, tol / r, length(ranges[i])) for i in 1:n]
    return _check_encoding(PatchyParticleSpecies{2,F,eltype(sites)}(g, sites, r, tol))
end

"""
    PatchySphere(p::Polyhedron, r=1; colors=1:nfaces(p), labels=colors, encoding=:auto)
    PatchySphere(sym::Symbol, n=0, r=1; kwargs...)

A 3D sphere of radius `r` carrying one patch per face of `p`, so that the patches inherit the
solid's rotation group. The second form names a rotation group instead of a solid, resolved
by [`Polyhedron`](@ref): `PatchySphere(:T)`, `PatchySphere(:O)`, `PatchySphere(:D, 5)`.

Patches sit where the face centroid directions pierce the sphere, and share the graph
encoding, the labelling rules and the binding site frame convention of
[`PolyhedronParticleSpecies`](@ref) — so a polyhedron species and a patchy sphere built from
the same solid have interchangeable encodings and can share one set of `BindingRules`.
"""
function PatchySphere(
    p::Polyhedron{F}, r::Real=1; colors=1:nfaces(p), labels=colors, encoding::Symbol=:auto
) where {F}
    n = nfaces(p)
    length(colors) == n ||
        throw(ArgumentError("expected $n colors, one per face, got $(length(colors))"))
    encoding in (:auto, :dart, :cycle) ||
        throw(ArgumentError("encoding must be :auto, :dart or :cycle, got :$encoding"))

    usecycle = encoding === :cycle || (encoding === :auto && allunique(labels))
    g, ranges = usecycle ? cycleencoding(n; labels) : dartencoding(p; labels)

    r = F(r)
    tol = sqrt(eps(F)) * r
    P = Pose{3,F,RotMatrix3{F}}
    sites = Vector{BindingSite{P,F}}(undef, n)
    for i in 1:n
        c = facecentroid(p, i)
        # On a sphere the patch normal has to be radial, so the frame is built from the
        # centroid direction and the tangential part of the direction to the first edge.
        ex = normalize(c)
        v = edgemidpoint(p, i, 1) - c
        ez = normalize(v - dot(v, ex) * ex)
        psi = RotMatrix3{F}(hcat(ex, cross(ez, ex), ez))
        sites[i] = BindingSite(P(r * ex, psi), colors[i], ranges[i], tol, tol / r, length(ranges[i]))
    end
    return _check_encoding(PatchyParticleSpecies{3,F,eltype(sites)}(g, sites, r, tol))
end

PatchySphere(sym::Symbol, n::Integer=0, r::Real=1; a=1.0, kwargs...) =
    PatchySphere(Polyhedron(sym, n; a), r; kwargs...)

function Base.show(io::Core.IO, ps::PatchyParticleSpecies{D}) where {D}
    return print(io, "$(D)d PatchyParticleSpecies with $(nsites(ps)) sites")
end

Base.copy(ps::PatchyParticleSpecies) =
    typeof(ps)(copy(ps.g), copy(ps.sites), ps.r, ps.skin)


graphrep(ps::PatchyParticleSpecies) = ps.g
nsites(ps::PatchyParticleSpecies) = length(ps.sites)
bindingsites(ps::PatchyParticleSpecies, i::Integer) = ps.sites[i]
isconvex(::PatchyParticleSpecies) = true

function setcolors!(ps::PatchyParticleSpecies, colors::AbstractVector{<:Integer})
    length(colors) != nsites(ps) && throw(ArgumentError("incorrect number of colors"))
    for k in eachindex(ps.sites)
        s = ps.sites[k]
        ps.sites[k] = BindingSite(s.pose, colors[k], s.vertices, s.touching_tolerance, s.alignment_tolerance, s.gauge)
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

bounding_radius(ps::PatchyParticleSpecies) = ps.r