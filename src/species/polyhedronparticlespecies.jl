"""
    PolyhedronParticleSpecies{F}

A 3D convex polyhedron with one binding site per face.
"""
struct PolyhedronParticleSpecies{F,B<:BindingSite} <: ParticleSpecies{3,B}
    g::NautyDiGraph
    sites::Vector{B}
    shape::Polyhedron{F}
    normals::Vector{SVector{3,F}}
    offsets::Vector{F}
    edgedirections::Vector{SVector{3,F}}
    centrosymmetric::Bool
    rmin::F
    rmax::F
    skin::F
end

"""
    PolyhedronParticleSpecies(p::Polyhedron; colors=1:nfaces(p), locking=true, twists=0)

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

The graph encoding is chosen rather than asked for: the sparse one whenever it provably carries
everything the dart encoding would, see `_cycle_suffices`.
"""
PolyhedronParticleSpecies(p::Polyhedron; colors=1:nfaces(p), locking=true, twists=0) =
    _polyhedronspecies(p, colors, locking, twists, nothing)

# `usecycle` forces an encoding; see `_facesites`. Internal, and the reason the public
# constructor above is a one-liner.
function _polyhedronspecies(p::Polyhedron{F}, colors, locking, twists,
                            usecycle::Union{Nothing,Bool}) where {F}
    n = nfaces(p)
    length(colors) == n ||
        throw(ArgumentError("expected $n colors, one per face, got $(length(colors))"))

    rmin = inradius(p)
    rmax = bounding_radius(p)
    tol = sqrt(eps(F)) * rmax

    g, sites = _facesites(p, fs -> _faceposes(p, fs), colors,
                          _perface(locking, n, "locking flags"),
                          _perface(twists, n, "twists"), usecycle, tol, tol / rmin)

    normals = facenormals(p)
    offsets = [dot(normals[i], facecentroid(p, i)) for i in 1:n]
    return _check_encoding(PolyhedronParticleSpecies{F,eltype(sites)}(
        g, sites, p, normals, offsets, _edgedirections(p), _iscentrosymmetric(p), rmin, rmax, tol
    ))
end

# Is the solid carried onto itself by inversion through its centre? Corners are centred, so this
# is just "every corner has an opposite". Cheap to ask once, and it is what makes the difference
# body of two translates equal to the solid at twice the size; see `translate_overlap`.
function _iscentrosymmetric(p::Polyhedron{F}) where {F}
    cs = corners(p)
    atol = sqrt(eps(F)) * maximum(norm, cs)
    return all(c -> any(c2 -> isapprox(-c, c2; atol), cs), cs)
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
    # The shape and everything derived from it are shared, not copied: they are immutable, and
    # `overlap` uses `shape ===` to recognize two particles of one species as translates.
    return typeof(ps)(
        copy(ps.g), copy(ps.sites), ps.shape, ps.normals, ps.offsets, ps.edgedirections,
        ps.centrosymmetric, ps.rmin, ps.rmax, ps.skin
    )
end

graphrep(p::PolyhedronParticleSpecies) = p.g
nsites(p::PolyhedronParticleSpecies) = length(p.sites)
bindingsites(p::PolyhedronParticleSpecies, i::Integer) = p.sites[i]
isconvex(::PolyhedronParticleSpecies) = true
bounding_radius(ps::PolyhedronParticleSpecies) = ps.rmax

"""
    _tiles(ps::PolyhedronParticleSpecies)

True for a cube whose binding sites are aligned with its own edges. See [`_tiles`](@ref) for
what this claims and what it is worth.

Cubes only, out of the several solids that tile space, because a cube is the one where nothing
else has to be checked. All six faces are congruent, so any bond the rules permit joins two
faces that are actually flush — a box with an `a×a` face and an `a×b` face would let a bond
place two boxes overlapping at a mismatched face, and no longer tiling. All six stabilisers are
4, so every bond has a single registration, and `locking` cannot open a second one that tips a
neighbour off the lattice the way a prism's square side face would.

That leaves the sites. A bond's relative rotation is fixed by the two frames, and it maps the
cubic lattice to itself exactly when each frame is built on the cube's own edges — so each
site's twist reference must point along an edge direction. Whole dart steps keep it there, and
`twists` of a fraction of a step do not, which is precisely when the partner arrives turned off
the lattice and free to overlap something.
"""
function _tiles(ps::PolyhedronParticleSpecies{F}) where {F}
    p = ps.shape
    atol = sqrt(eps(F)) * ps.rmax
    ncorners(p) == 8 && nfaces(p) == 6 && ps.centrosymmetric || return false
    # Equidistant corners and equidistant faces: among the boxes, that is the cube.
    all(c -> isapprox(norm(c), ps.rmax; atol), corners(p)) || return false
    all(d -> isapprox(d, ps.rmin; atol), ps.offsets) || return false
    # And every site's twist reference lies along an edge, so no bond leaves the lattice.
    return all(1:nsites(ps)) do i
        ez = bindingsites(ps, i).pose.psi[:, 3]
        any(e -> isapprox(abs(dot(e, ez)), 1; atol=sqrt(eps(F))), ps.edgedirections)
    end
end

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
    skin = spcs1.skin + spcs2.skin
    t = pose2.x - pose1.x
    # Outer and inner spheres first; both are measured from the origin, which the centred
    # corners of a `Polyhedron` make the solid's centroid.
    d = norm(t)
    d >= spcs1.rmax + spcs2.rmax && return false
    d < (spcs1.rmin + spcs2.rmin) - skin && return true

    # Two identically oriented particles of one centrally symmetric species are translates of
    # each other, and the whole question is then one point against the faces of the doubled
    # solid — O(#faces), against a candidate set that is quadratic in the edges. This is the
    # case a lattice is made of, and cubes, boxes, even prisms and most Platonic solids qualify.
    spcs1.shape === spcs2.shape && spcs1.centrosymmetric && _sameorientation(pose1, pose2) &&
        return translate_overlap(spcs1.normals, spcs1.offsets, pose1.psi \ t, skin)

    # Separating axes. 3D needs more candidates than 2D: the face normals of both solids, plus
    # the cross products of their edge directions, which catch the edge-on-edge configurations
    # no face normal separates.
    axes = Iterators.flatten((
        (pose1.psi * nrm for nrm in spcs1.normals),
        (pose2.psi * nrm for nrm in spcs2.normals),
        (cross(pose1.psi * e1, pose2.psi * e2)
         for e1 in spcs1.edgedirections, e2 in spcs2.edgedirections),
    ))
    return sat_overlap(axes, corners(spcs1), pose1, corners(spcs2), pose2, skin)
end

# Do the two particles carry the same orientation, so that one is a translate of the other?
# Tight on purpose: the fast path it guards is exact only for an exact translate, and a pair
# that misses it merely pays the full axis set.
#
# Central symmetry is not a technicality here, and dropping it was tried and measured wrong. The
# faces of a Minkowski difference body `K ⊕ (−K)` come from three sources, and only two of them
# are faces of `K`: the third is *edge against edge*, which is why separating axes need cross
# products at all. A tetrahedron's difference body is a cuboctahedron, whose six square faces
# belong to no face of the tetrahedron — so for a solid that is not centrally symmetric, its own
# face normals do not decide even the aligned case. When `K = −K` the difference body is `2K`
# and the third source contributes nothing, which is exactly the condition tested here.
@inline _sameorientation(pose1::Pose{D,F}, pose2::Pose{D,F}) where {D,F} =
    isapprox(pose1.psi, pose2.psi; atol=sqrt(eps(F)), rtol=zero(F))

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
