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
    colors=1:length(patch_positions),
)
    D = length(first(patch_positions))
    F = float(eltype(first(patch_positions)))
    r = F(r)
    n = length(patch_positions)

    tol = sqrt(eps(F)) * r
    poses = [normal_pose(patch_positions[i], patch_twists[i]) for i in eachindex(patch_positions)]
    # One vertex per patch leaves no room for a turn about the patch normal, so the sitesym is 1
    # whatever the arrangement; a patch that needs its own symmetry needs `dartencoding`.
    sitesyms = ones(Int, n)
    labs = siteorbits(poses, sitesyms, collect(colors))
    stabs = sitestabilizers(poses, sitesyms, labs)
    setlabels!(g, collect(Cint, labs))

    sites = [BindingSite(poses[i], colors[i], i:i, tol, tol / r, 1, stabs[i]) for i in 1:n]
    return check_encoding(PatchyParticleSpecies{D,F,eltype(sites)}(g, sites, r, tol))
end

"""
    PatchyDisk(angles, r=1; colors=1:length(angles))

A 2D disk of radius `r` with one binding site placed at each angle in `angles` (radians,
measured counterclockwise from the +x axis).
"""
function PatchyDisk(angles, r=1; colors=1:length(angles))
    F = float(eltype(angles))
    r = F(r)
    n = length(angles)
    tol = sqrt(eps(F)) * r
    positions = [SVector(r * cos(F(phi)), r * sin(F(phi))) for phi in angles]
    poses = [normal_pose(positions[i], F(0)) for i in 1:n]

    # A 2D site has no turn about its in-plane normal, so its sitesym is 1 throughout.
    sitesyms = ones(Int, n)
    labels = siteorbits(poses, sitesyms, collect(colors))
    stabs = sitestabilizers(poses, sitesyms, labels)

    g, ranges = cycleencoding(n; labels)
    sites = [BindingSite(poses[i], colors[i], ranges[i], tol, tol / r, 1, stabs[i]) for i in 1:n]
    return check_encoding(PatchyParticleSpecies{2,F,eltype(sites)}(g, sites, r, tol))
end

"""
    PatchySphere(p::Polyhedron, r=1; colors=1:nfaces(p), locking=true, twists=0)

A 3D sphere of radius `r` carrying one patch per face of polyhedron `p`, so that the patches inherit the
polyhedron's rotation group.

See [`PolyhedronParticleSpecies`](@ref) for documentation of the keyword arguments.
"""
function PatchySphere(p::Polyhedron, r::Real=1; colors=1:nfaces(p), locking=true, twists=0)
    _patchysphere(p, r, colors, locking, twists, nothing)
end

"""
    _reseat_radially(pose, r)

Move a face's binding site frame onto the sphere of radius `r`: the site sits at `r` along the
centroid direction with a radial normal, keeping its twist reference by re-orthogonalising the
face's local z. The two frames coincide whenever the centroid direction is the face normal, as
on any Platonic solid.
"""
function _reseat_radially(pose::Pose{3,F}, r::Real) where {F}
    ex = normalize(pose.x)
    ez = pose.psi[:, 3]
    ez = normalize(ez - dot(ez, ex) * ex)
    return Pose{3,F,RotMatrix3{F}}(F(r) * ex, RotMatrix3{F}(hcat(ex, cross(ez, ex), ez)))
end

function _patchysphere(p::Polyhedron{F}, r::Real, colors, locking, twists, usecycle::Union{Nothing,Bool}) where {F}
    n = nfaces(p)
    length(colors) == n || throw(ArgumentError("expected $n colors, one per face, got $(length(colors))"))

    r = F(r)
    tol = sqrt(eps(F)) * r
    patchposes(fs) = [_reseat_radially(q, r) for q in _faceposes(p, fs)]

    g, sites = _facesites(
        p,
        patchposes,
        colors,
        _expandperface(locking, n, "locking flags"),
        _expandperface(twists, n, "twists"),
        usecycle,
        tol,
        tol / r,
    )
    return check_encoding(PatchyParticleSpecies{3,F,eltype(sites)}(g, sites, r, tol))
end

function Base.show(io::Core.IO, ps::PatchyParticleSpecies{D}) where {D}
    return print(io, "$(D)d PatchyParticleSpecies with $(nsites(ps)) sites")
end

Base.copy(ps::PatchyParticleSpecies) = typeof(ps)(copy(ps.g), copy(ps.sites), ps.r, ps.skin)

graphrep(ps::PatchyParticleSpecies) = ps.g
nsites(ps::PatchyParticleSpecies) = length(ps.sites)
bindingsites(ps::PatchyParticleSpecies, i::Integer) = ps.sites[i]
isconvex(::PatchyParticleSpecies) = true

function could_contact(::SpeciesAndPose{<:PatchyParticleSpecies}, ::SpeciesAndPose{<:PatchyParticleSpecies}; kwargs...)
    return true # this check would be identical with overlap, so no need to do it twice
end

function overlap(p1::SpeciesAndPose{<:PatchyParticleSpecies}, p2::SpeciesAndPose{<:PatchyParticleSpecies}; kwargs...)
    spcs1, pose1 = p1
    spcs2, pose2 = p2
    return norm(pose1.x - pose2.x) < (spcs1.r + spcs2.r) - (spcs1.skin + spcs2.skin)
end

bounding_radius(ps::PatchyParticleSpecies) = ps.r