"""
    PolyhedronParticleSpecies{F}

A 3D convex polyhedron with one binding site per face.
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
    PolyhedronParticleSpecies(p::Polyhedron; colors=1:nfaces(p), locking=true, twists=0, encoding=:auto)

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

`twists` turns a face's binding site about its own normal, by an angle in radians, and likewise
takes one value or one per face. This is the bond *registry*: which relative orientation a bond
means, as opposed to how many it admits. Note that turning both faces of a bond by the same
amount does not cancel, it turns the partner by twice that — the offset appears on both sides
of the face-to-face flip, and `Δ·Rx(-θ) = Rx(θ)·Δ`.

Any angle is allowed. Whole dart steps, `2π/degree`, are taken by rotating the face's corner
list so that the frame and the vertex numbering move together; whatever is left over turns the
frame alone. Twisting one face of a symmetry orbit differently from its fellows splits the
orbit, lowering the symmetry number — deliberately breaking a symmetry is what this is for, and
the graph records the break. A twist shared across an orbit leaves the symmetry intact, since
turns about a site's own normal commute with its stabiliser.
"""
function PolyhedronParticleSpecies(
    p::Polyhedron{F}; colors=1:nfaces(p), locking=true, twists=0, encoding::Symbol=:auto
) where {F}
    n = nfaces(p)
    length(colors) == n ||
        throw(ArgumentError("expected $n colors, one per face, got $(length(colors))"))
    encoding in (:auto, :dart, :cycle) ||
        throw(ArgumentError("encoding must be :auto, :dart or :cycle, got :$encoding"))

    rmin = inradius(p)
    rmax = bounding_radius(p)
    tol = sqrt(eps(F)) * rmax

    g, sites = _facesites(p, fs -> _faceposes(p, fs), colors,
                          _perface(locking, n, "locking flags"),
                          _perface(twists, n, "twists"), encoding, tol, tol / rmin)

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

function overlap(
    p1::SpeciesAndPose{<:PolyhedronParticleSpecies},
    p2::SpeciesAndPose{<:PolyhedronParticleSpecies};
    kwargs...,
)
    spcs1, pose1 = p1
    spcs2, pose2 = p2
    # Outer and inner spheres first; both are measured from the origin, which the centred
    # corners of a `Polyhedron` make the solid's centroid.
    d = norm(pose1.x - pose2.x)
    d >= spcs1.rmax + spcs2.rmax && return false
    d < (spcs1.rmin + spcs2.rmin) - (spcs1.skin + spcs2.skin) && return true

    # Separating axes. 3D needs more candidates than 2D: the face normals of both solids, plus
    # the cross products of their edge directions, which catch the edge-on-edge configurations
    # no face normal separates.
    axes = Iterators.flatten((
        (pose1.psi * nrm for nrm in spcs1.normals),
        (pose2.psi * nrm for nrm in spcs2.normals),
        (cross(pose1.psi * e1, pose2.psi * e2)
         for e1 in spcs1.edgedirections, e2 in spcs2.edgedirections),
    ))
    return sat_overlap(axes, corners(spcs1), pose1, corners(spcs2), pose2,
                       spcs1.skin + spcs2.skin)
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
