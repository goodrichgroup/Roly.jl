"""
    PolyhedronParticleSpecies{F}

A 3D convex polyhedron with one binding site at the centre of each face.
"""
struct PolyhedronParticleSpecies{F,B<:BindingSite} <: ParticleSpecies{3,B}
    g::NautyDiGraph
    sites::Vector{B}
    shape::Polyhedron{F}
    normals::Vector{SVector{3,F}}
    edgedirections::Vector{SVector{3,F}}
    rmin::F
    rmax::F
    skin::F
end

"""
    PolyhedronParticleSpecies(p::Polyhedron; colors=1:nfaces(p), locking=true, encoding=:auto)

Build a particle species from the polyhedron `p`, with one binding site at each face centroid.

`colors` assigns interaction colors to the binding sites, and is all that needs saying: the
particle's symmetry follows from it, since two faces are interchangeable exactly when a
rotation of the solid carries one onto the other and they are the same color (see
[`siteorbits`](@ref)).

So `colors=1:nfaces(p)` gives every face its own identity and a symmetry number of 1;
`colors=fill(1, nfaces(p))` makes them all alike and recovers the solid's full rotation group;
and colouring the caps of a cube apart from its sides leaves the subgroup that preserves that
split.

`locking` says whether a site holds its partner in the orientation its frame names, and takes
either one flag for the whole species or one per face. The default, `true`, is the ordinary
reading of an oriented binding site, and leaves a bond with a single registration unless the
particle's own symmetry makes the frame ambiguous. Setting a face rotation-free instead admits
every orientation the face geometrically permits: a triangular prism with square sides is only
2-fold about them, so `locking=true` bonds its prisms coplanar, while freeing a side face also
allows the neighbour stood on its side. See [`nregistrations`](@ref).
"""
function PolyhedronParticleSpecies(
    p::Polyhedron{F}; colors=1:nfaces(p), locking=true, encoding::Symbol=:auto
) where {F}
    n = nfaces(p)
    length(colors) == n ||
        throw(ArgumentError("expected $n colors, one per face, got $(length(colors))"))
    encoding in (:auto, :dart, :cycle) ||
        throw(ArgumentError("encoding must be :auto, :dart or :cycle, got :$encoding"))
    locking = _locking(locking, n)

    rmin = inradius(p)
    rmax = bounding_radius(p)
    tol = sqrt(eps(F)) * rmax

    # The twist references are settled in two passes. The solid has already picked each face's
    # first corner from its own shape; that is enough to compute the labelling, since it fixes
    # the frames up to a gauge turn and nothing here looks finer. Knowing the labelling then
    # gives the symmetry group, and re-winding along it pins the references up to `stab`, which
    # is what bonds need. See `_propagate_faces`.
    fs = faces(p)
    poses = _faceposes(p, fs)
    gauges = facegauge(p)
    labels = siteorbits(poses, gauges, collect(colors))
    fs = _propagate_faces(corners(p), fs, labels)
    poses = _faceposes(p, fs)
    stabs = sitestabilisers(poses, gauges, labels)

    usecycle = encoding === :cycle ||
        (encoding === :auto && _cycle_suffices(_twistfreedoms(gauges, stabs, locking), labels))
    g, ranges = usecycle ? cycleencoding(n; labels) : dartencoding(fs; labels)

    sites = [BindingSite(poses[i], colors[i], ranges[i], tol, tol / rmin, gauges[i], stabs[i],
                         locking[i]) for i in 1:n]

    return _check_encoding(PolyhedronParticleSpecies{F,eltype(sites)}(
        g, sites, p, facenormals(p), _edgedirections(p), rmin, rmax, tol
    ))
end

"""
    _faceposes(p::Polyhedron, fs)

The binding site frames of `p`'s faces, given face corner lists `fs`: local x along the
outward normal, local z pointing at the midpoint of each face's first edge.

`fs` is passed separately rather than read from `p` because [`_propagate_faces`](@ref) re-winds
the lists, and the frames and the encoding have to be built from the same ones.
"""
function _faceposes(p::Polyhedron{F}, fs::Vector{Vector{Int}}) where {F}
    P = Pose{3,F,RotMatrix3{F}}
    return map(eachindex(fs)) do i
        x = facecentroid(p, i)
        ex = facenormal(p, i)
        ez = normalize((corners(p)[fs[i][1]] + corners(p)[fs[i][2]]) / 2 - x)
        P(x, RotMatrix3{F}(hcat(ex, cross(ez, ex), ez)))
    end
end

# Separating axis candidates contributed by the edges. Only the direction matters, and only
# up to sign, so parallel edges are collapsed: a cube contributes 3 rather than 12.
function _edgedirections(p::Polyhedron{F}) where {F}
    dirs = SVector{3,F}[]
    for f in faces(p), k in eachindex(f)
        d = normalize(corners(p)[f[mod1(k + 1, length(f))]] - corners(p)[f[k]])
        any(e -> isapprox(abs(dot(e, d)), 1; atol=sqrt(eps(F))), dirs) && continue
        push!(dirs, d)
    end
    return dirs
end

function Base.show(io::Core.IO, ps::PolyhedronParticleSpecies)
    return print(io, "$(dimension(ps))d PolyhedronParticleSpecies with $(nsites(ps)) sites")
end

function Base.copy(ps::PolyhedronParticleSpecies)
    return typeof(ps)(
        copy(ps.g), copy(ps.sites), ps.shape, ps.normals, ps.edgedirections, ps.rmin, ps.rmax, ps.skin
    )
end

graphrep(p::PolyhedronParticleSpecies) = p.g
nsites(p::PolyhedronParticleSpecies) = length(p.sites)
bindingsites(p::PolyhedronParticleSpecies, i::Integer) = p.sites[i]
isconvex(::PolyhedronParticleSpecies) = true
bounding_radius(ps::PolyhedronParticleSpecies) = ps.rmax

"""
    shape(ps::PolyhedronParticleSpecies)

Return the [`Polyhedron`](@ref) the species was built from.
"""
shape(ps::PolyhedronParticleSpecies) = ps.shape
corners(ps::PolyhedronParticleSpecies) = corners(ps.shape)

function setcolors!(p::PolyhedronParticleSpecies, colors::AbstractVector{<:Integer})
    _recolor!(p, p.sites, colors)
    return nothing
end

function _SAT_overlap(
    p1::SpeciesAndPose{<:PolyhedronParticleSpecies}, p2::SpeciesAndPose{<:PolyhedronParticleSpecies}
)
    spcs1, pose1 = p1
    spcs2, pose2 = p2

    # The 3D separating axis theorem needs more candidate axes than the 2D one: the face
    # normals of both solids, plus the cross products of their edge directions, which catch
    # the edge-on-edge configurations no face normal separates.
    skin_sum = spcs1.skin + spcs2.skin
    corners1, corners2 = corners(spcs1), corners(spcs2)

    function separated(axis)
        n2 = dot(axis, axis)
        n2 < eps(typeof(n2)) && return false
        axis = axis / sqrt(n2)
        lo1, hi1 = extrema(dot(axis, pose1 * c) for c in corners1)
        lo2, hi2 = extrema(dot(axis, pose2 * c) for c in corners2)
        return hi2 < lo1 + skin_sum || hi1 < lo2 + skin_sum
    end

    for nrm in spcs1.normals
        separated(pose1.psi * nrm) && return false
    end
    for nrm in spcs2.normals
        separated(pose2.psi * nrm) && return false
    end
    for e1 in spcs1.edgedirections, e2 in spcs2.edgedirections
        separated(cross(pose1.psi * e1, pose2.psi * e2)) && return false
    end
    return true
end

function overlap(
    p1::SpeciesAndPose{<:PolyhedronParticleSpecies},
    p2::SpeciesAndPose{<:PolyhedronParticleSpecies};
    kwargs...,
)
    spcs1, pose1 = p1
    spcs2, pose2 = p2
    d = norm(pose1.x - pose2.x)
    # check outer and inner spheres first 
    d >= spcs1.rmax + spcs2.rmax && return false
    d < (spcs1.rmin + spcs2.rmin) - (spcs1.skin + spcs2.skin) && return true
    return _SAT_overlap(p1, p2)
end

"""
    UnitTetrahedron

A regular tetrahedron with unit-length edges and one binding site per face.
"""
const UnitTetrahedron = PolyhedronParticleSpecies(Tetrahedron())

"""
    UnitCube

A cube with unit-length edges and one binding site per face.
"""
const UnitCube = PolyhedronParticleSpecies(Cube())

"""
    UnitOctahedron

A regular octahedron with unit-length edges and one binding site per face.
"""
const UnitOctahedron = PolyhedronParticleSpecies(Octahedron())

"""
    UnitDodecahedron

A regular dodecahedron with unit-length edges and one binding site per face.
"""
const UnitDodecahedron = PolyhedronParticleSpecies(Dodecahedron())

"""
    UnitIcosahedron

A regular icosahedron with unit-length edges and one binding site per face.
"""
const UnitIcosahedron = PolyhedronParticleSpecies(Icosahedron())

"""
    UnitPyramid(n, a=1.0; h=a, kwargs...)

A pyramid over a regular `n`-gon with edge length `a`, with one binding site per face.
Rotation group `C_n`.
"""
UnitPyramid(n::Integer, a::Real=1.0; h::Real=a, kwargs...) =
    PolyhedronParticleSpecies(Pyramid(n, a; h); kwargs...)

"""
    UnitPrism(n, a=1.0; h=a, kwargs...)

A prism over a regular `n`-gon with edge length `a`, with one binding site per face.
Rotation group `D_n`, except for `UnitPrism(4)` whose default height makes it a cube.
"""
UnitPrism(n::Integer, a::Real=1.0; h::Real=a, kwargs...) =
    PolyhedronParticleSpecies(Prism(n, a; h); kwargs...)

"""
    UnitAntiprism(n, a=1.0; kwargs...)

A uniform antiprism over a regular `n`-gon with edge length `a`, with one binding site per
face. Rotation group `D_n`, except for `UnitAntiprism(3)` which is a regular octahedron.
"""
UnitAntiprism(n::Integer, a::Real=1.0; kwargs...) = PolyhedronParticleSpecies(Antiprism(n, a); kwargs...)
