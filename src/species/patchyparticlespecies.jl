"""
    PatchyParticleSpecies{D,F,B}

A `D`-dimensional sphere of radius `r` with binding sites on its surface.
"""
struct PatchyParticleSpecies{D,F,B<:BindingSite} <: ParticleSpecies{D,B}
    g::NautyDiGraph
    sites::Vector{B}
    r::F
    skin::F
end

"""
    PatchyParticleSpecies(g::NautyDiGraph, r, poses; colors=eachindex(poses))

Construct a `D`-dimensional sphere of radius `r` with binding sites on its surface.

General constructor, requiring manual encoding of the particle's symmetry into a `NautyDiGraph`.
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
    n = length(patch_positions)

    tol = sqrt(eps(F)) * r
    poses = [normal_pose(patch_positions[i], patch_twists[i]) for i in eachindex(patch_positions)]
    # One vertex per patch leaves no room for a turn about the patch normal, so the gauge is 1
    # whatever the arrangement; a patch that needs its own symmetry needs `dartencoding`.
    gauges = ones(Int, n)
    labs = siteorbits(poses, gauges, collect(colors))
    stabs = sitestabilisers(poses, gauges, labs)
    setlabels!(g, collect(Cint, labs))

    sites = [BindingSite(poses[i], colors[i], i:i, tol, tol / r, 1, stabs[i]) for i in 1:n]
    return _check_encoding(PatchyParticleSpecies{D,F,eltype(sites)}(g, sites, r, tol))
end

"""
    PatchyDisk(angles, r=1; colors=1:length(angles))

A 2D disk of radius `r` with one binding site placed at each angle in `angles` (radians,
measured counterclockwise from the +x axis).

`colors` assigns interaction colors to the patches, and the symmetry follows from it: two
patches are interchangeable exactly when a rotation carries one onto the other and they are the
same color (see [`siteorbits`](@ref)).
"""
function PatchyDisk(angles, r=1; colors=1:length(angles))
    F = float(eltype(angles))
    r = F(r)
    n = length(angles)
    tol = sqrt(eps(F)) * r
    positions = [SVector(r * cos(F(phi)), r * sin(F(phi))) for phi in angles]
    poses = [normal_pose(positions[i], F(0)) for i in 1:n]

    # A 2D site has no turn about its in-plane normal, so its gauge is 1 throughout.
    gauges = ones(Int, n)
    labels = siteorbits(poses, gauges, collect(colors))
    stabs = sitestabilisers(poses, gauges, labels)

    g, ranges = cycleencoding(n; labels)
    sites = [BindingSite(poses[i], colors[i], ranges[i], tol, tol / r, 1, stabs[i]) for i in 1:n]
    return _check_encoding(PatchyParticleSpecies{2,F,eltype(sites)}(g, sites, r, tol))
end

"""
    PatchySphere(p::Polyhedron, r=1; colors=1:nfaces(p), locking=true, twists=0)
    PatchySphere(group::RotationGroup, r=1; a=1.0, kwargs...)

A 3D sphere of radius `r` carrying one patch per face of polyhedron `p`, so that the patches inherit the
polyhedron's rotation group.

See [`PolyhedronParticleSpecies`](@ref) for documentation of the keyword arguments.
"""
PatchySphere(p::Polyhedron, r::Real=1; colors=1:nfaces(p), locking=true, twists=0) =
    _patchysphere(p, r, colors, locking, twists, nothing)

PatchySphere(group::RotationGroup, r::Real=1; a=1.0, kwargs...) =
    PatchySphere(Polyhedron(group; a), r; kwargs...)

function _patchysphere(p::Polyhedron{F}, r::Real, colors, locking, twists, usecycle::Union{Nothing,Bool}) where {F}
    n = nfaces(p)
    length(colors) == n ||
        throw(ArgumentError("expected $n colors, one per face, got $(length(colors))"))

    r = F(r)
    tol = sqrt(eps(F)) * r
    # On a sphere the patch normal has to be radial, so the frame is built from the centroid
    # direction and the tangential part of the direction to the face's first edge. Everything
    # else about the construction is shared with the polyhedron species; see `_facesites`.
    P = Pose{3,F,RotMatrix3{F}}
    patchposes(fs) = map(eachindex(fs)) do i
        c = facecentroid(p, i)
        ex = normalize(c)
        v = (corners(p)[fs[i][1]] + corners(p)[fs[i][2]]) / 2 - c
        ez = normalize(v - dot(v, ex) * ex)
        P(r * ex, RotMatrix3{F}(hcat(ex, cross(ez, ex), ez)))
    end

    g, sites = _facesites(p, patchposes, colors, _perface(locking, n, "locking flags"),
                          _perface(twists, n, "twists"), usecycle, tol, tol / r)
    return _check_encoding(PatchyParticleSpecies{3,F,eltype(sites)}(g, sites, r, tol))
end

function Base.show(io::Core.IO, ps::PatchyParticleSpecies{D}) where {D}
    return print(io, "$(D)d PatchyParticleSpecies with $(nsites(ps)) sites")
end

Base.copy(ps::PatchyParticleSpecies) =
    typeof(ps)(copy(ps.g), copy(ps.sites), ps.r, ps.skin)


graphrep(ps::PatchyParticleSpecies) = ps.g
nsites(ps::PatchyParticleSpecies) = length(ps.sites)
bindingsites(ps::PatchyParticleSpecies, i::Integer) = ps.sites[i]
isconvex(::PatchyParticleSpecies) = true

function could_contact(
    ::SpeciesAndPose{<:PatchyParticleSpecies}, ::SpeciesAndPose{<:PatchyParticleSpecies}; kwargs...
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