"""
    RotationGroup

A proper (or *chiral*) point group: the rotations carrying a rigid body onto itself, without
the reflections a full point group would also contain. These are the only groups a particle can
have, which is why Roly names them and not the point groups.

The five families are [`Cyclic`](@ref), [`Dihedral`](@ref), [`Tetrahedral`](@ref),
[`Octahedral`](@ref) and [`Icosahedral`](@ref). A `RotationGroup` is a way of *asking* for a
shape — `Polyhedron(Octahedral())` is a cube — rather than something computed from one; use
[`rotationgroup`](@ref) to get the actual rotations of a solid.

Spelled out rather than abbreviated to the customary `C`, `D`, `T`, `O`, `I`, all of which are
taken: `T` is the universal type-parameter convention, `I` is `LinearAlgebra.I`, and `D` is
Roly's dimension parameter throughout (`ParticleSpecies{D,B}`, `Pose{D,F}`). Exporting the
short names would shadow them on `using Roly`.
"""
abstract type RotationGroup end

"""
    Cyclic(n)

The rotation group `C_n`: `n` turns about a single axis, and nothing else. Order `n`, realized
by [`Pyramid`](@ref).
"""
struct Cyclic{N} <: RotationGroup end
Cyclic(n::Integer) = Cyclic{Int(n)}()

"""
    Dihedral(n)

The rotation group `D_n`: the `n` turns about a principal axis together with the `n` half turns
about axes perpendicular to it. Order `2n`, realized by [`Prism`](@ref) and [`Antiprism`](@ref).
"""
struct Dihedral{N} <: RotationGroup end
Dihedral(n::Integer) = Dihedral{Int(n)}()

"""
    Tetrahedral()

The rotation group `T` of a regular tetrahedron, of order 12.
"""
struct Tetrahedral <: RotationGroup end

"""
    Octahedral()

The rotation group `O` shared by the cube and the regular octahedron, of order 24.
"""
struct Octahedral <: RotationGroup end

"""
    Icosahedral()

The rotation group `I` shared by the regular dodecahedron and icosahedron, of order 60.
"""
struct Icosahedral <: RotationGroup end

Base.show(io::Core.IO, ::Cyclic{N}) where {N} = print(io, "Cyclic($N)")
Base.show(io::Core.IO, ::Dihedral{N}) where {N} = print(io, "Dihedral($N)")

"""
    grouporder(group::RotationGroup)

Return the number of rotations in `group`. A solid realizing it has this symmetry number when
all of its faces are alike, which is what `length(rotationgroup(Polyhedron(group)))` measures
the hard way.
"""
grouporder(::Cyclic{N}) where {N} = N
grouporder(::Dihedral{N}) where {N} = 2N
grouporder(::Tetrahedral) = 12
grouporder(::Octahedral) = 24
grouporder(::Icosahedral) = 60

"""
    Polyhedron{F}

A convex polyhedron, stored as a list of `corners` and a list of `faces`.

Each face is a list of indices into `corners`, wound counter-clockwise as seen from
*outside* the solid. The winding is what fixes the orientation of the graph encoding, and is
checked on construction, as is convexity; each face is then rotated to start at a canonical
corner, see [`_canonical_faces`](@ref). Species refine that further with
[`_propagate_faces`](@ref), which needs the labelling and so cannot happen here.

**Corners are translated to put their centroid at the origin**, and everything downstream may
assume it. Several quantities are measured from the origin rather than from the solid —
[`bounding_radius`](@ref) is `maximum(norm, corners)` and [`inradius`](@ref) the smallest face
centroid norm, both of which feed the overlap fast path — and rotations are about the origin
too. Normalising once at the door is what makes those readings mean what they say; the
alternative is every one of them recomputing a centroid, and being silently wrong if one
forgets.

Convexity is not decoration either. [`PolyhedronParticleSpecies`](@ref) reports `isconvex` for
every solid and tests overlap by separating axes, which is a convex-only argument, and
[`site_symmetry`](@ref) relies on a solid being the intersection of its faces' half-spaces.
"""
struct Polyhedron{F<:AbstractFloat}
    corners::Vector{SVector{3,F}}
    faces::Vector{Vector{Int}}

    function Polyhedron(corners::Vector{SVector{3,F}}, faces::Vector{Vector{Int}}) where {F}
        _check_winding(length(corners), faces)
        corners = _recenter(corners)
        _check_convex(corners, faces)
        return new{F}(corners, _canonical_faces(corners, faces))
    end
end

"""
    Polyhedron(corners, faces)
    Polyhedron(corners)

Construct a polyhedron from its `corners` and `faces`. Faces must be wound counter-clockwise
as seen from outside the solid; this is verified, since a reversed face would silently
produce a graph encoding with the wrong symmetry group.

If `faces` is omitted, they are derived from the corners by finding all supporting planes of
the convex hull, which requires the solid to be convex. That is the form to reach for, and
together with the library constructors ([`Cube`](@ref), [`Prism`](@ref), …) covers most uses.

Giving `faces` explicitly is the expert path, and what it buys is control the derived form
does not offer: a face's *first corner* is the twist reference of the binding site placed on
it, so the ordering within each face is a modelling choice. It is only a starting point —
[`_canonical_faces`](@ref) rotates each face to an intrinsic corner and species then re-wind
along their symmetry group — so it is not a way to set a bond registry, which is what the
`twists` keyword of [`PolyhedronParticleSpecies`](@ref) is for.
"""
function Polyhedron(corners::AbstractVector{<:AbstractVector}, faces)
    F = float(eltype(first(corners)))
    return Polyhedron(
        SVector{3,F}[SVector{3,F}(c) for c in corners], Vector{Int}[collect(Int, f) for f in faces]
    )
end
function Polyhedron(corners::AbstractVector{<:AbstractVector})
    F = float(eltype(first(corners)))
    cs = SVector{3,F}[SVector{3,F}(c) for c in corners]
    return Polyhedron(cs, _derive_faces(cs))
end

"""
    _check_winding(ncorners, faces)

Verify that `faces` describes a closed, consistently oriented surface: every directed edge
occurs exactly once, and its reverse occurs in exactly one other face. This is precisely the
condition that makes the edge pairing of [`dartencoding`](@ref) well defined.

Also that every corner is used by some face. A corner strictly inside the hull lies on no
supporting plane, so [`_derive_faces`](@ref) never mentions it and it would otherwise be
carried along silently — contributing nothing to the shape while still counting towards
[`bounding_radius`](@ref) and shifting the centroid the corners are recentred on.
"""
function _check_winding(ncorners::Integer, faces::AbstractVector{<:AbstractVector{Int}})
    length(faces) < 3 && throw(ArgumentError("a polyhedron needs at least 3 faces"))
    used = falses(ncorners)
    for f in faces, v in f
        checkbounds(Bool, used, v) && (used[v] = true)
    end
    all(used) || throw(ArgumentError(
        "corner$(count(!, used) > 1 ? "s" : "") $(join(findall(!, used), ", ")) " *
        "$(count(!, used) > 1 ? "are" : "is") used by no face; a corner strictly inside the " *
        "hull is not part of the solid, and one on the boundary means the face list is incomplete"
    ))
    seen = Dict{Tuple{Int,Int},Int}()
    for (i, f) in enumerate(faces)
        length(f) < 3 && throw(ArgumentError("face $i has fewer than 3 corners"))
        allunique(f) || throw(ArgumentError("face $i repeats a corner"))
        all(in(1:ncorners), f) || throw(ArgumentError("face $i indexes a nonexistent corner"))
        for k in eachindex(f)
            e = (f[k], f[mod1(k + 1, length(f))])
            haskey(seen, e) && throw(ArgumentError(
                "directed edge $e occurs in faces $(seen[e]) and $i; faces must all be wound counter-clockwise seen from outside"
            ))
            seen[e] = i
        end
    end
    for (e, i) in seen
        haskey(seen, reverse(e)) || throw(ArgumentError(
            "edge $e of face $i is not shared with a second face; the surface is not closed"
        ))
    end
    return nothing
end

"""
    _check_convex(corners, faces)

Verify that every corner lies on the inner side of every face's plane, so that the solid is
the intersection of its faces' half-spaces. See [`Polyhedron`](@ref) for what depends on it.
"""
function _check_convex(corners::Vector{SVector{3,F}}, faces::AbstractVector{<:AbstractVector{Int}}) where {F}
    atol = sqrt(eps(F)) * maximum(norm, corners)
    for (i, f) in enumerate(faces)
        c = sum(corners[v] for v in f) / length(f)
        nrm = zero(SVector{3,F})
        for k in eachindex(f)
            nrm += cross(corners[f[k]] - c, corners[f[mod1(k + 1, length(f))]] - c)
        end
        nrm = normalize(nrm)
        d = dot(nrm, c)
        for (j, x) in enumerate(corners)
            dot(nrm, x) <= d + atol || throw(ArgumentError(
                "corner $j lies outside the plane of face $i, so the solid is not convex; " *
                "Roly's polyhedra must be convex, since overlap is tested by separating axes"
            ))
        end
    end
    return nothing
end

"""
    _canonical_faces(corners, faces)

Return `faces` with each corner list cyclically rotated so that the face begins at an
intrinsically chosen corner, fixing the twist reference of that face's binding site.

A site's local z points at the midpoint of its face's *first* edge, so which corner comes
first is the site's twist reference. That reference is not physical, but it is not free
either: [`site_symmetry`](@ref) compares two faces' frames up to a turn by `2π/gauge`, so
faces that a rotation of the solid carries onto one another must pick *corresponding* first
edges. A square face has gauge 4 and every corner is corresponding, which is why the choice
never mattered for cubes; a rectangle has gauge 2 and degree 4, its corners falling into two
classes a quarter turn apart, and picking a long edge on one face and a short edge on another
hides the symmetry that relates them.

The choice is made from the face alone: rotate to the lexicographically least cyclic word of
`(edge length, interior angle)`. That word determines a planar polygon up to congruence and
its cyclic symmetries are exactly the face's own rotations, so the minimising positions form
exactly one gauge orbit. Congruent faces have identical words, so their minimising positions
correspond up to a gauge turn — which is the tolerance `site_symmetry` allows. A regular face
ties everywhere and the rotation is a no-op.

This is only half of what the twist references have to satisfy, and the weaker half; see
[`_propagate_faces`](@ref), which runs after it and overrides it wherever the two disagree.
"""
function _canonical_faces(corners::Vector{SVector{3,F}}, faces::Vector{Vector{Int}}) where {F}
    atol = sqrt(eps(F)) * maximum(norm, corners)
    return map(faces) do f
        k = length(f)
        edge(m) = norm(corners[f[mod1(m + 1, k)]] - corners[f[mod1(m, k)]])
        function turn(m)
            a = corners[f[mod1(m - 1, k)]] - corners[f[mod1(m, k)]]
            b = corners[f[mod1(m + 1, k)]] - corners[f[mod1(m, k)]]
            return acos(clamp(dot(a, b) / (norm(a) * norm(b)), -one(F), one(F)))
        end
        # Does the word starting at s come before the one starting at t?
        function isbefore(s, t)
            for j in 0:(k - 1), (x, y) in ((edge(s + j), edge(t + j)), (turn(s + j), turn(t + j)))
                isapprox(x, y; atol) || return x < y
            end
            return false
        end
        best = 1
        for s in 2:k
            isbefore(s, best) && (best = s)
        end
        circshift(f, 1 - best)
    end
end

"""
    _propagate_faces(corners, faces, labels)

Return `faces` with each corner list cyclically rotated so that faces related by a
label-preserving rotation of the solid carry twist references related by that same rotation.

This is the condition [`_canonical_faces`](@ref) cannot reach and the one bonds actually need.
That rule reads a face's own shape, so it can only line references up *to the face's gauge* —
every dart of a square is like every other, so it has nothing to say. But two frames differing
by a gauge turn describe the same bond only when that turn is also a symmetry of the whole
particle, i.e. lies in `stab`. Where `gauge > stab` the leftover turns are real, and choosing
them inconsistently across a face orbit makes otherwise equivalent bonds differ.

A triangular prism with square sides is the case in point: `gauge` 4, `stab` 2. Left to the
face rule, one side face's reference comes out a quarter turn from its neighbours' — legal by
gauge, not by stab — and bonding face 2 to face 3 then turns the neighbour out of the plane
where face 2 to face 4 does not, so the prisms stop tiling as polyiamonds.

The group has to be the *label-preserving* rotations, not all of them. A rotation outside it
relates two faces by an element of the solid's stabiliser rather than the species', which is
too coarse by exactly the amount that goes wrong above: a cube with four sticky sides and two
caps has all six faces in one orbit under its 24 rotations, but only the 8 preserving
sides-from-caps are symmetries of the particle, and propagating under the 24 leaves adjacent
side faces a quarter turn apart.

Correctness, and that the slack is exactly `stab` rather than merely `gauge`, is proved in
`docs/src/orientation.md`. The ordering it forces is worth noting: labels are needed to know the
group, and propagation only ever moves a dart within its gauge orbit, so `facegauge` and
`siteorbits` — which compare frames up to gauge — may be computed first and are unaffected.
"""
function _propagate_faces(cs::Vector{SVector{3,F}}, faces::Vector{Vector{Int}},
                          labels) where {F}
    atol = sqrt(eps(F)) * maximum(norm, cs)
    centroid(f) = sum(cs[v] for v in f) / length(f)
    midpoint(f, k) = (cs[f[k]] + cs[f[mod1(k + 1, length(f))]]) / 2

    # The rotations that carry every face onto one of the same label. Composition and inverses
    # preserve that, so these form a group, which the proof needs.
    faceof(x) = findfirst(i -> isapprox(centroid(faces[i]), x; atol), eachindex(faces))
    group = filter(_rotationgroup(cs, faces)) do Q
        all(eachindex(faces)) do i
            j = faceof(Q * centroid(faces[i]))
            !isnothing(j) && labels[j] == labels[i]
        end
    end
    length(group) == 1 && return faces

    faces = [copy(f) for f in faces]
    for j in 2:length(faces)
        cj = centroid(faces[j])
        for i in 1:(j - 1)
            labels[i] == labels[j] || continue
            # A symmetry carrying face i's centroid onto face j's carries the face itself, and
            # so carries its first edge midpoint onto one of face j's.
            k = findfirst(Q -> isapprox(Q * centroid(faces[i]), cj; atol), group)
            isnothing(k) && continue
            target = group[k] * midpoint(faces[i], 1)
            d = findfirst(m -> isapprox(midpoint(faces[j], m), target; atol), eachindex(faces[j]))
            isnothing(d) && continue
            faces[j] = circshift(faces[j], 1 - d)
            break
        end
    end
    return faces
end

"""
    _derive_faces(corners)

Find the faces of the convex hull of `corners` by testing every corner triple for a
supporting plane, then winding each face counter-clockwise about its outward normal.

Measured from the corner centroid, so that a plane's offset `d` is its signed distance from the
inside of the solid: `d > 0` says the normal points outward, and `dot(nrm, r) <= d` for every
corner says the plane supports the hull. (The corners are not yet centred here — the
`Polyhedron` constructor does that afterwards, on the way in.)

Every corner on a face's plane joins that face, including one sitting mid-edge, which gives a
face with three collinear corners. That is a face of degree `k+1` where the shape has `k`
sides, so it raises `nedges` and the dart count with it — arguably what a combinatorial
encoding should do with a subdivided face, but worth knowing before subdividing one.
"""
function _derive_faces(corners::Vector{SVector{3,F}}) where {F}
    n = length(corners)
    center = sum(corners) / n
    rel = [c - center for c in corners]
    atol = sqrt(eps(F)) * maximum(norm, rel)

    faces = Vector{Int}[]
    planes = Tuple{SVector{3,F},F}[]
    for i in 1:(n - 2), j in (i + 1):(n - 1), k in (j + 1):n
        nrm = cross(rel[j] - rel[i], rel[k] - rel[i])
        # Collinear triple: no plane through it, and any plane it lies on is found by another.
        norm(nrm) < atol && continue
        nrm = normalize(nrm)
        d = dot(nrm, rel[i])
        # Point the normal away from the centroid, which is inside the hull.
        d < 0 && ((nrm, d) = (-nrm, -d))
        # Supporting plane of the hull?
        all(r -> dot(nrm, r) <= d + atol, rel) || continue
        any(pl -> isapprox(pl[1], nrm; atol) && isapprox(pl[2], d; atol), planes) && continue

        push!(planes, (nrm, d))
        push!(faces, _wind_counterclockwise!(
            corners, findall(r -> abs(dot(nrm, r) - d) <= atol, rel), nrm
        ))
    end
    return faces
end

"""
    _wind_counterclockwise!(corners, idxs, nrm)

Sort `idxs` — the corners lying on one supporting plane — into counter-clockwise order seen
from outside, i.e. looking down `nrm`. In place, since the caller hands over a fresh `findall`.

The angle about the face centroid is what orders them, and `atan` is the honest way to get it:
a half-plane comparator avoids the transcendental but is more code and easier to get wrong,
for no measurable gain at these sizes.
"""
function _wind_counterclockwise!(corners::Vector{SVector{3,F}}, idxs::Vector{Int},
                                 nrm::SVector{3,F}) where {F}
    c = sum(corners[i] for i in idxs) / length(idxs)
    u = normalize(corners[first(idxs)] - c)
    v = cross(nrm, u)
    return sort!(idxs; by=i -> atan(dot(corners[i] - c, v), dot(corners[i] - c, u)))
end

"""
    corners(p::Polyhedron)

Return the corner positions of `p`. Named `corners` rather than `vertices` to keep the
Graphs.jl meaning of `vertices` free for the graph encodings.
"""
corners(p::Polyhedron) = p.corners

"""
    faces(p::Polyhedron)

Return the faces of `p` as lists of corner indices, each wound counter-clockwise seen from
outside.
"""
faces(p::Polyhedron) = p.faces

"""
    facevertices(p::Polyhedron, i)

Return the corner indices of the `i`th face of `p`.
"""
facevertices(p::Polyhedron, i::Integer) = p.faces[i]

ncorners(p::Polyhedron) = length(p.corners)

"""
    nfaces(p::Polyhedron)

Return the number of faces of `p`, which is the number of binding sites of a species built
from it.
"""
nfaces(p::Polyhedron) = length(p.faces)

"""
    nedges(p::Polyhedron)

Return the number of edges of `p`. The dart encoding has `2 * nedges(p)` vertices.
"""
nedges(p::Polyhedron) = sum(length, p.faces) ÷ 2

facedegree(p::Polyhedron, i::Integer) = length(p.faces[i])

Base.eltype(::Type{<:Polyhedron{F}}) where {F} = F
Base.eltype(p::Polyhedron) = eltype(typeof(p))

"""
    facecentroid(p::Polyhedron, i)
    facecentroids(p::Polyhedron)

Return the centroid of the `i`th face of `p`, or of all faces. This is where a binding site
of a [`PolyhedronParticleSpecies`](@ref) sits.
"""
facecentroid(p::Polyhedron, i::Integer) = sum(p.corners[v] for v in p.faces[i]) / facedegree(p, i)
facecentroids(p::Polyhedron) = [facecentroid(p, i) for i in 1:nfaces(p)]

"""
    facenormal(p::Polyhedron, i)
    facenormals(p::Polyhedron)

Return the outward unit normal of the `i`th face of `p`, or of all faces. The direction
follows from the counter-clockwise winding of the face.
"""
function facenormal(p::Polyhedron{F}, i::Integer) where {F}
    f = p.faces[i]
    c = facecentroid(p, i)
    nrm = zero(SVector{3,F})
    for k in eachindex(f)
        nrm += cross(p.corners[f[k]] - c, p.corners[f[mod1(k + 1, length(f))]] - c)
    end
    return normalize(nrm)
end
facenormals(p::Polyhedron) = [facenormal(p, i) for i in 1:nfaces(p)]

"""
    edgemidpoint(p::Polyhedron, i, k)

Return the midpoint of edge `k` of face `i` of `p`. Both arguments are indices into `p`, the
first over its faces and the second over that face's own edges, of which edge `k` is the one
running from the face's `k`th corner to its `k+1`th — not an edge from corner `i` to corner
`k`. The `k = 1` midpoint is the twist reference of the face's binding site.
"""
function edgemidpoint(p::Polyhedron, i::Integer, k::Integer)
    f = p.faces[i]
    return (p.corners[f[k]] + p.corners[f[mod1(k + 1, length(f))]]) / 2
end

"""
    bounding_radius(p::Polyhedron)

Return the distance from the origin to the farthest corner of `p`.
"""
bounding_radius(p::Polyhedron) = maximum(norm, p.corners)

"""
    inradius(p::Polyhedron)

Return the distance from the origin to the closest face centroid of `p`.
"""
inradius(p::Polyhedron) = minimum(norm(facecentroid(p, i)) for i in 1:nfaces(p))

"""
    minedgelength(p::Polyhedron)

Return the length of the shortest edge of `p`. Named for what it returns: the edges of a
`Polyhedron` need not all be the same length, and only the library solids built to an edge
length `a` have a single one.
"""
function minedgelength(p::Polyhedron)
    return minimum(
        norm(p.corners[f[k]] - p.corners[f[mod1(k + 1, length(f))]]) for f in p.faces for k in eachindex(f)
    )
end

"""
    rotationgroup(p::Polyhedron)

Return the proper rotations that map `p` onto itself, as `RotMatrix3`es about the corner
centroid. `length(rotationgroup(p))` is the solid's true symmetry number.

Every symmetry maps the first dart to some dart, and that correspondence determines the
rotation, so it suffices to test the `2 * nedges(p)` candidates this generates.
"""
rotationgroup(p::Polyhedron) = _rotationgroup(corners(p), faces(p))

# Same, on the raw corner and face lists, so the Polyhedron constructor can use it before there
# is a Polyhedron. Rotations are about the origin, which the centred-corners invariant makes the
# corner centroid; see `Polyhedron`.
function _rotationgroup(cs::Vector{SVector{3,F}}, faces::Vector{Vector{Int}}) where {F}
    atol = sqrt(eps(F)) * maximum(norm, cs)

    # An orthonormal frame attached to dart k of face i, i.e. to that face's kth edge: along
    # the edge, along the outward normal, and their cross product.
    function dartframe(i, k)
        f = faces[i]
        e1 = normalize(cs[f[mod1(k + 1, length(f))]] - cs[f[k]])
        e2 = _facenormal(cs, f)
        return hcat(e1, e2, cross(e1, e2))
    end

    Mref = dartframe(1, 1)
    group = RotMatrix3{F}[]
    for i in eachindex(faces), k in eachindex(faces[i])
        # every potential symmetry rotation maps the ref frame to some other dartframe, and is therefore
        # given by the matrix
        R = dartframe(i, k) * Mref'
        # R is a symmetry if every corners is mapped to some other corner
        all(c -> any(c2 -> isapprox(R * c, c2; atol), cs), cs) || continue
        push!(group, RotMatrix3{F}(R))
    end
    return group
end

function _facenormal(cs::Vector{SVector{3,F}}, f::Vector{Int}) where {F}
    c = sum(cs[v] for v in f) / length(f)
    nrm = zero(SVector{3,F})
    for k in eachindex(f)
        nrm += cross(cs[f[k]] - c, cs[f[mod1(k + 1, length(f))]] - c)
    end
    return normalize(nrm)
end

"""
    geometriclabels(p::Polyhedron)

Return one label per face, grouping faces that are equivalent under [`rotationgroup`](@ref).

Passing these to [`dartencoding`](@ref) yields the solid's true rotation group, whereas
`labels=fill(1, nfaces(p))` yields the *combinatorial* symmetry of the face lattice, which
can be larger. For example, `Pyramid(3)` is combinatorially a tetrahedron and would report 12 instead of
its actual 3.
"""
function geometriclabels(p::Polyhedron)
    cents = facecentroids(p)
    c0 = sum(corners(p)) / ncorners(p)
    cents = [c - c0 for c in cents]
    atol = sqrt(eps(eltype(p))) * maximum(norm, cents)

    labels = zeros(Int, nfaces(p))
    group = rotationgroup(p)
    next = 0
    for i in eachindex(labels)
        labels[i] == 0 || continue
        next += 1
        for R in group
            j = findfirst(c -> isapprox(R * cents[i], c; atol), cents)
            isnothing(j) || (labels[j] = next)
        end
    end
    return labels
end

"""
    facegauge(p::Polyhedron, i)
    facegauge(p::Polyhedron)

Return the order of face `i`'s own rotational symmetry about its outward normal, or one entry
per face.

This is the `gauge` of a binding site placed on that face: a face invariant under turns of
`2π/gauge` has that many equally good twist references, so nothing about the particle in
isolation may depend on which was picked. It divides the face degree and equals it only for a
regular face — a rectangular face has degree 4 but gauge 2.

It is a property of the face alone, not of the solid it belongs to. A triangular prism is only
2-fold about a square side face, but the face is still a square, so its gauge is 4.
"""
function facegauge(p::Polyhedron{F}, i::Integer) where {F}
    f = facevertices(p, i)
    k = length(f)
    c = facecentroid(p, i)
    nrm = facenormal(p, i)
    rel = [corners(p)[v] - c for v in f]
    atol = sqrt(eps(F)) * maximum(norm, rel)
    # A rotational symmetry of a k-gon shifts its corner ring cyclically, and the shift by s is
    # realized by the turn 2πs/k about the normal.
    return count(0:(k - 1)) do s
        R = AngleAxis(2F(π) * s / k, nrm[1], nrm[2], nrm[3])
        all(j -> isapprox(R * rel[j], rel[mod1(j + s, k)]; atol), 1:k)
    end
end
facegauge(p::Polyhedron) = [facegauge(p, i) for i in 1:nfaces(p)]

function Base.show(io::Core.IO, p::Polyhedron)
    return print(io, "Polyhedron[v=$(ncorners(p)), e=$(nedges(p)), f=$(nfaces(p))]")
end
function Base.show(io::Core.IO, ::MIME"text/plain", p::Polyhedron{F}) where {F}
    println(io, "Polyhedron{$F}:")
    println(io, " - corners: \t$(ncorners(p))")
    println(io, " - edges: \t$(nedges(p))")
    print(io, " - faces: \t$(nfaces(p)) $(_facesummary(p))")
end
function _facesummary(p::Polyhedron)
    degs = sort!(unique(facedegree(p, i) for i in 1:nfaces(p)))
    return "(" * join(("$(count(i -> facedegree(p, i) == d, 1:nfaces(p)))×$d-gon" for d in degs), ", ") * ")"
end

"""
    dartencoding(p::Polyhedron; labels=1:nfaces(p))
    dartencoding(faces::Vector{Vector{Int}}; labels=1:length(faces))
    dartencoding(group::RotationGroup; a=1.0, labels=...)

Return `(g, ranges)`, the dart encoding of `p` and the graph vertices belonging to each face.

The second form takes the face lists directly, for species that have re-wound them with
[`_propagate_faces`](@ref): a site's vertex range has to start at the dart its frame points
at, so the encoding and the poses must be built from the same lists. The third names a
[`RotationGroup`](@ref) instead of a solid, and encodes whichever solid [`Polyhedron`](@ref)
realizes it with.

A *dart* is a half-edge: one of the two directed traversals of a polyhedron edge, the one
belonging to a given face. Each edge therefore has two darts, one per adjoining face, a face of
degree `k` owns `k` of them, and the graph has `2 * nedges(p)` vertices in total. The encoding
consists of

1. a directed `k`-cycle through each face's own darts, following the face's counter-clockwise
   winding, and
2. a bidirectional edge joining the two darts that sit on each shared polyhedron edge.

The automorphism group of the resulting `g` is the polyhedron's rotational symmetry group.
`labels` assigns one symmetry label per face, inherited by that face's darts: all labels
distinct gives a symmetry number of 1, all labels equal gives the full rotation group, and
merging some faces gives the subgroup preserving that labelling.

For the Platonic solids this is a Cayley digraph of that group: the rotations act freely and
transitively on the darts, so the darts *are* the group elements and `|Aut| = 2 * nedges(p)`
— 12 = 2·6, 24 = 2·12, 60 = 2·30. It stops being one as soon as the solid is less symmetric
than its face lattice: a prism has `|G| = 2n` but `2E = 6n`.
"""
dartencoding(p::Polyhedron; labels=1:nfaces(p)) = dartencoding(faces(p); labels)
dartencoding(group::RotationGroup; a=1.0, kwargs...) = dartencoding(Polyhedron(group; a); kwargs...)

function dartencoding(fs::Vector{Vector{Int}}; labels=1:length(fs))
    length(labels) == length(fs) ||
        throw(ArgumentError("expected $(length(fs)) labels, one per face, got $(length(labels))"))

    # graph vertices for each face of the polyhedron
    ranges = Vector{UnitRange{Int}}(undef, length(fs))
    ndarts = 0
    for (i, f) in enumerate(fs)
        ranges[i] = (ndarts + 1):(ndarts + length(f))
        ndarts += length(f)
    end

    vertex_labels = Cint[labels[i] for (i, f) in enumerate(fs) for _ in f]
    g = NautyDiGraph(ndarts; vertex_labels)

    # Dart k of face i sits on the directed edge f[k] -> f[k+1]. Its partner is the dart of
    # the adjacent face carrying the reversed edge
    dart_of_edge = Dict{Tuple{Int,Int},Int}()
    for (i, f) in enumerate(fs), k in eachindex(f)
        dart_of_edge[(f[k], f[mod1(k + 1, length(f))])] = first(ranges[i]) + k - 1
    end

    for (i, f) in enumerate(fs)
        v0 = first(ranges[i])
        for k in eachindex(f)
            v = v0 + k - 1
            add_edge!(g, v, v0 + mod(k, length(f)))
            partner = dart_of_edge[(f[mod1(k + 1, length(f))], f[k])]
            if v < partner
                add_edge!(g, v, partner)
                add_edge!(g, partner, v)
            end
        end
    end
    return g, ranges
end

"""
    cycleencoding(nsites; labels=1:nsites)

Return `(g, ranges)`, the cycle encoding of a particle with `nsites` binding sites: a single
directed cycle carrying one vertex per site.

This is the 2D polygon encoding, and it is also valid in 3D under the conditions
[`_cycle_suffices`](@ref) states. It is *not* valid in general.
"""
function cycleencoding(nsites::Integer; labels=1:nsites)
    length(labels) == nsites ||
        throw(ArgumentError("expected $nsites labels, one per site, got $(length(labels))"))
    nsites < 1 && throw(ArgumentError("a particle needs at least one binding site"))

    g = NautyDiGraph(cycle_digraph(nsites); vertex_labels=collect(Cint, labels))
    return g, [i:i for i in 1:nsites]
end

"""
    _twistfreedoms(gauges, stabs, locking)
    _perface(x, n, what)

Per-site [`twistfreedom`](@ref) from the two symmetry counts and the locking flags, and the
normalisation of a keyword given either as a single value for the whole species or as one value
per site.
"""
_twistfreedoms(gauges, stabs, locking) = [l ? s : g for (g, s, l) in zip(gauges, stabs, locking)]

function _perface(x, n::Integer, what::AbstractString)
    x isa Union{Bool,Real} && return fill(x, n)
    length(x) == n ||
        throw(ArgumentError("expected $n $what, one per face, got $(length(x))"))
    return collect(x)
end

"""
    _facesites(p, poseof, colors, locking, twists, usecycle, touching_tol, alignment_tol)

Build the graph and binding sites of a species carrying one site per face of `p`, shared by
[`PolyhedronParticleSpecies`](@ref) and [`PatchySphere`](@ref). They differ only in where a
site sits and how its frame is built, which is what `poseof` supplies: it maps a list of face
corner lists to the poses of the corresponding sites.

`usecycle` picks the encoding: `nothing` takes the cheap one whenever [`_cycle_suffices`](@ref)
says it is equivalent, which is what the species constructors pass, and `true`/`false` force
one. Forcing is not a modelling choice and so is not a keyword of any public constructor — it
exists for the tests, which check that the two agree wherever both are valid and that the cheap
one is rejected where it is not.

The order is forced. A face's first corner is its site's twist reference, and settling it takes
three steps, each needing the one before:

1. The solid arrives with each face already started at an intrinsically chosen corner
   ([`_canonical_faces`](@ref)), which pins the references up to each face's `gauge`.
2. That is enough to derive the labelling, since [`siteorbits`](@ref) compares frames only up
   to `gauge`. Knowing the labelling gives the symmetry group.
3. [`_propagate_faces`](@ref) then re-winds along that group, pinning the references up to
   `stab`, which is the strength bonds need.

`twists` is applied last, as a per-face angle about the site's own normal. It is folded into
the key `siteorbits` groups by, alongside the color, so that twisting one face of an orbit
differently from its fellows splits that orbit rather than silently breaking the equivariance
propagation just established. Splitting is always safe: the graph then distinguishes more, not
less; and a twist shared across an orbit preserves equivariance, since turns about a site's own
normal commute with its stabiliser.
"""
function _facesites(p::Polyhedron{F}, poseof, colors, locking, twists,
                    usecycle::Union{Nothing,Bool}, touching_tol::Real,
                    alignment_tol::Real) where {F}
    n = nfaces(p)
    gauges = facegauge(p)
    labels = siteorbits(poseof(faces(p)), gauges, collect(zip(colors, twists)))

    fs = _propagate_faces(corners(p), faces(p), labels)
    # A twist is an angle about the site's own normal. Whole dart steps are taken by rotating
    # the face's corner list, which moves the frame and the vertex numbering together and keeps
    # `contact_pairing`'s anchor on a pair of coincident darts; the remainder is applied to the
    # frame alone. That remainder is fine: the anchor has to be consistent and to tell
    # registrations apart, not to mark a physical coincidence, and both survive.
    steps = [round(Int, t * length(f) / (2F(π))) for (f, t) in zip(fs, twists)]
    fs = [circshift(f, -mod(m, length(f))) for (f, m) in zip(fs, steps)]
    poses = [pose * RotX(F(t) - 2F(π) * m / length(f))
             for (pose, t, m, f) in zip(poseof(fs), twists, steps, fs)]
    stabs = sitestabilisers(poses, gauges, labels)

    cyclic = something(usecycle, _cycle_suffices(_twistfreedoms(gauges, stabs, locking), labels))
    g, ranges = cyclic ? cycleencoding(n; labels) : dartencoding(fs; labels)

    sites = [BindingSite(poses[i], colors[i], ranges[i], touching_tol, alignment_tol,
                         gauges[i], stabs[i], locking[i]) for i in 1:n]
    return g, sites
end

"""
    _cycle_suffices(twistfreedoms, labels)

Whether one graph vertex per site carries everything the encoding has to carry, so that
[`cycleencoding`](@ref) can stand in for [`dartencoding`](@ref).

Two things have to fit. The *symmetry number* must come out right, which needs all `labels`
distinct: no rotation can then preserve the labelling, so the answer is 1 whatever structure
the graph has internally, and a bare cycle gives 1. And the *bonds* must be distinguishable,
which needs every site's twist freedom to be 1. A one-vertex site pins no turn about its own
normal, so a bond to it is a single graph edge however the partner is turned; a site with `q`
registrations would have all `q` canonicalise alike and be merged, losing `q - 1` real
structures rather than merely encoding them coarsely. See [`contact_pairing`](@ref).

Distinct labels already force every stabiliser to 1, so for locking sites the second condition
comes free and this reduces to the first. It bites only on rotation-free sites, which are
asking for exactly the registrations a single vertex cannot record.
"""
_cycle_suffices(qs, labels) = allunique(labels) && all(isone, qs)


"""
    Polyhedron(group::RotationGroup; a=1.0)

Return a polyhedron realizing the proper rotation `group`, with edge length `a`:

| `group`         | solid            | order |
|:----------------|:-----------------|:------|
| `Cyclic(n)`     | `Pyramid(n)`     | `n`   |
| `Dihedral(n)`   | `Prism(n)`       | `2n`  |
| `Tetrahedral()` | `Tetrahedron()`  | 12    |
| `Octahedral()`  | `Cube()`         | 24    |
| `Icosahedral()` | `Dodecahedron()` | 60    |

One realization per group, since a group does not determine a shape. The alternatives
`Octahedron()`, `Icosahedron()` and `Antiprism(n)` have the same groups and are constructed by
name.
"""
Polyhedron(::Cyclic{N}; a=1.0) where {N} = Pyramid(N, a)
Polyhedron(::Dihedral{N}; a=1.0) where {N} = Prism(N, a)
Polyhedron(::Tetrahedral; a=1.0) = Tetrahedron(a)
Polyhedron(::Octahedral; a=1.0) = Cube(a)
Polyhedron(::Icosahedral; a=1.0) = Dodecahedron(a)

################################################################################
# Solid library
################################################################################

# The corner centroid is fixed by every symmetry, so placing it at the origin makes the whole
# rotation group act about the particle's own origin.
_recenter(cs) = (c0 = sum(cs) / length(cs); [c - c0 for c in cs])

function _rescale(cs, a)
    s = a / minimum(norm(cs[i] - cs[j]) for i in 1:(length(cs) - 1) for j in (i + 1):length(cs))
    return _recenter([s * c for c in cs])
end

"""
    Tetrahedron(a=1.0)

A regular tetrahedron with edge length `a`. Proper rotation group `T`, of order 12.
"""
function Tetrahedron(a::Real=1.0)
    F = float(typeof(a))
    cs = SVector{3,F}[(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]
    return Polyhedron(_rescale(cs, a))
end

"""
    Cube(a=1.0)

A cube with edge length `a`. Proper rotation group `O`, of order 24.
"""
function Cube(a::Real=1.0)
    F = float(typeof(a))
    cs = SVector{3,F}[(x, y, z) for x in (-1, 1) for y in (-1, 1) for z in (-1, 1)]
    return Polyhedron(_rescale(cs, a))
end

"""
    Octahedron(a=1.0)

A regular octahedron with edge length `a`. Proper rotation group `O`, of order 24.
"""
function Octahedron(a::Real=1.0)
    F = float(typeof(a))
    cs = SVector{3,F}[(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]
    return Polyhedron(_rescale(cs, a))
end

"""
    Icosahedron(a=1.0)

A regular icosahedron with edge length `a`. Proper rotation group `I`, of order 60.
"""
function Icosahedron(a::Real=1.0)
    F = float(typeof(a))
    φ = (1 + sqrt(F(5))) / 2
    cs = SVector{3,F}[]
    for s1 in (-1, 1), s2 in (-1, 1)
        push!(cs, SVector{3,F}(0, s1, s2 * φ), SVector{3,F}(s1, s2 * φ, 0), SVector{3,F}(s2 * φ, 0, s1))
    end
    return Polyhedron(_rescale(cs, a))
end

"""
    Dodecahedron(a=1.0)

A regular dodecahedron with edge length `a`. Proper rotation group `I`, of order 60.
"""
function Dodecahedron(a::Real=1.0)
    F = float(typeof(a))
    φ = (1 + sqrt(F(5))) / 2
    cs = SVector{3,F}[SVector{3,F}(x, y, z) for x in (-1, 1) for y in (-1, 1) for z in (-1, 1)]
    for s1 in (-1, 1), s2 in (-1, 1)
        push!(
            cs,
            SVector{3,F}(0, s1 / φ, s2 * φ),
            SVector{3,F}(s1 / φ, s2 * φ, 0),
            SVector{3,F}(s2 * φ, 0, s1 / φ),
        )
    end
    return Polyhedron(_rescale(cs, a))
end

"""
    Pyramid(n, a=1.0; h=a)

A pyramid over a regular `n`-gon with edge length `a` and apex height `h`. Proper rotation
group `C_n`, of order `n`.
"""
function Pyramid(n::Integer, a::Real=1.0; h::Real=a)
    n >= 3 || throw(ArgumentError("a pyramid needs n >= 3, got n=$n"))
    F = float(promote_type(typeof(a), typeof(h)))
    r = F(a) / (2 * sin(F(π) / n))
    cs = SVector{3,F}[SVector{3,F}(r * cos(2F(π) * k / n), r * sin(2F(π) * k / n), 0) for k in 0:(n - 1)]
    push!(cs, SVector{3,F}(0, 0, F(h)))
    return Polyhedron(_recenter(cs))
end

"""
    Prism(n, a=1.0; h=a)

A prism over a regular `n`-gon with edge length `a` and height `h`. Proper rotation group
`D_n`, of order `2n`. `Prism(4)` is a cube.
"""
function Prism(n::Integer, a::Real=1.0; h::Real=a)
    n >= 3 || throw(ArgumentError("a prism needs n >= 3, got n=$n"))
    F = float(promote_type(typeof(a), typeof(h)))
    r = F(a) / (2 * sin(F(π) / n))
    cs = SVector{3,F}[]
    for s in (-1, 1), k in 0:(n - 1)
        push!(cs, SVector{3,F}(r * cos(2F(π) * k / n), r * sin(2F(π) * k / n), s * F(h) / 2))
    end
    return Polyhedron(cs)
end

"""
    Antiprism(n, a=1.0)

A uniform antiprism over a regular `n`-gon with edge length `a`. Proper rotation group `D_n`,
of order `2n`. `Antiprism(3)` is an octahedron.
"""
function Antiprism(n::Integer, a::Real=1.0)
    n >= 3 || throw(ArgumentError("an antiprism needs n >= 3, got n=$n"))
    F = float(typeof(a))
    r = F(a) / (2 * sin(F(π) / n))
    h = F(a) * sqrt(1 - 1 / (4 * cos(F(π) / 2n)^2))
    cs = SVector{3,F}[]
    for k in 0:(n - 1)
        push!(cs, SVector{3,F}(r * cos(2F(π) * k / n), r * sin(2F(π) * k / n), -h / 2))
    end
    for k in 0:(n - 1)
        θ = 2F(π) * (k + F(1) / 2) / n
        push!(cs, SVector{3,F}(r * cos(θ), r * sin(θ), h / 2))
    end
    return Polyhedron(cs)
end

################################################################################
# Encoding validation
################################################################################

"""
    _siteturns(psi, gauge)

The `gauge` orientations a site is indistinguishable between: turns by `2π/gauge` about its
own outward normal, which `normal_pose` puts on the site's local x axis.

In 2D there is no rotation about an in-plane axis, so a site has exactly one orientation and
this is a singleton by construction rather than by arithmetic.
"""
_siteturns(psi::Rotation{3,F}, gauge::Integer) where {F} =
    (psi * RotX(F(2π) * m / gauge) for m in 0:(gauge - 1))
_siteturns(psi::Rotation{2}, ::Integer) = (psi,)

"""
    site_symmetry(ps::ParticleSpecies)

Return the number of rotations about the particle origin that map every binding site of `ps`
onto a binding site with the same symmetry label, matching orientation as well as position.

Frames need only agree up to the receiving site's own turns, see [`_siteturns`](@ref).
Demanding them on the nose would report 6 for a cube rather than 24, since a quarter turn
about a face normal maps that face to itself but sends its first dart to the second.

The turn count comes from each site's `gauge`, never from its graph vertex count. Those
coincide for the dart encoding, where a face gets one vertex per dart *because* it is that
symmetric, but taking the vertex count would let an encoding certify itself: one vertex per
face declares a pentagonal base 1-fold, exact frame matching then returns 1, the graph also
says 1, and a combination that should be rejected passes.

A rotation is pinned by where it sends one site's frame, so every symmetry arises from exactly
one `(site, turn)` pair that could receive site 1. That also means the answer is always finite:
a frame has no continuous stabiliser, so even a single site, or two antipodal ones, give a
definite count.
"""
site_symmetry(ps::ParticleSpecies) = length(_site_symmetries(_sitedata(ps)...))

# Poses, gauges and the key each site is matched by. `site_symmetry` keys on the graph's
# labels, since it asks what the graph claims; `siteorbits` keys on colors, since it asks what
# the arrangement is.
function _sitedata(ps::ParticleSpecies)
    n = nsites(ps)
    sites = [bindingsites(ps, i) for i in 1:n]
    labs = labels(graphrep(ps))
    return ([s.pose for s in sites], [s.gauge for s in sites],
            [labs[first(s.vertices)] for s in sites])
end

"""
    _site_symmetries(poses, gauges, keys)

Return the site permutations induced by the rotations about the particle origin that carry
every site onto one with the same `keys` entry, matching position and orientation.

!!! warning "Particles only"
    This rests on a particle's sites having *distinct positions*, which makes a rotation
    determined by where it sends site 1, and makes the site map injective for free. Neither
    holds for an assembled [`Polyform`](@ref): a bond puts two sites in the same place, so one
    rotation is recorded once per site coincident with site 1 — reporting 16 for a two-cube
    assembly whose shape has 8 — and a non-injective map can satisfy the pointwise conditions.
    Counting an assembly's symmetries needs deduplication by rotation and an explicit
    bijectivity test, and then the permutation is no longer well defined, which is why the
    check lives in `test/symmetry.jl` rather than here.
"""
function _site_symmetries(poses, gauges, keys)
    n = length(poses)
    tol = sqrt(eps(eltype(typeof(first(poses)))))
    atol = tol * maximum(norm(p.x) for p in poses)

    function permutation(Q)
        perm = zeros(Int, n)
        for i in 1:n
            j = findfirst(1:n) do j
                keys[j] == keys[i] &&
                    isapprox(Q * poses[i].x, poses[j].x; atol) &&
                    any(psi -> isapprox(Q * poses[i].psi, psi; atol=tol),
                        _siteturns(poses[j].psi, gauges[j]))
            end
            isnothing(j) && return nothing
            perm[i] = j
        end
        return perm
    end

    perms = Vector{Int}[]
    for a in 1:n
        keys[a] == keys[1] || continue
        for psi in _siteturns(poses[a].psi, gauges[a])
            perm = permutation(psi * inv(poses[1].psi))
            isnothing(perm) || push!(perms, perm)
        end
    end
    return perms
end

"""
    siteorbits(poses, gauges, colors)

Group the sites into the orbits of the rotations that preserve the *colored* arrangement, and
return one orbit index per site.

This is the labelling a species should carry: two sites are interchangeable exactly when some
rotation of the particle carries one onto the other and they are the same color. Deriving it
means labels never have to be given — placing colors on sites is the physical statement, and
the symmetry follows from the geometry. It also makes the labelling automatically no finer than
the coloring, so a graph can never merge two sites that bond differently.
"""
function siteorbits(poses, gauges, colors)
    n = length(poses)
    orbit = collect(1:n)
    for perm in _site_symmetries(poses, gauges, colors), i in 1:n
        lo, hi = minmax(orbit[i], orbit[perm[i]])
        hi == lo && continue
        replace!(orbit, hi => lo)
    end
    # Compact to 1:k so the result can be used as graph labels directly.
    ids = sort!(unique(orbit))
    return [searchsortedfirst(ids, o) for o in orbit]
end

"""
    sitestabilisers(ps::ParticleSpecies)
    sitestabilisers(poses, gauges, keys)

Return, per site, how many of the particle's own symmetries leave that site where it is.

Where [`site_symmetry`](@ref) counts all of them, this counts the ones fixing each site — the
turns about that site's normal that carry the whole particle onto itself. It is at most the
site's `gauge`, and usually less: a triangular prism is 2-fold about a square side face, so
that face has gauge 4 but a stabiliser of 2, while a cube's face has both equal to 4.

The difference is what decides how many *distinct* ways a partner can attach there, see
[`nregistrations`](@ref). Species compute this once at construction and store it on each
`BindingSite`, since it is fixed as soon as the labelling is; the `ps` method reads that back.

Keyed on the graph's labels, which is exact because a labelling is required to be at least as
fine as the coloring (see [`_check_labelling`](@ref)). A turn that leaves every label alone
therefore leaves every color alone too, so it changes neither which structure this is nor which
bonds the attached particle offers next.
"""
sitestabilisers(ps::ParticleSpecies) = [bindingsites(ps, i).stab for i in 1:nsites(ps)]

function sitestabilisers(poses, gauges, keys)
    perms = _site_symmetries(poses, gauges, keys)
    return [count(perm -> perm[i] == i, perms) for i in eachindex(poses)]
end

"""
    _recolor!(sites, g, colors)

Give `sites` the interaction colors `colors`, and bring the labelling and the stabilisers back
into step with them.

A coloring is the whole statement: which sites are interchangeable follows from it and the
geometry, and so does how many ways a partner can attach there. Leaving either behind would let
a recolored species keep a symmetry it no longer has — the bug `setcolors!` used to have, since
it rewrote colors alone.

Labels are written through each site's `vertices`, which is sound because a species' graph is
never canonised in place; see [`symmetrynumber`](@ref).

The graph's *structure* is not rebuilt, only its labels, so a recoloring needing a different
encoding cannot be applied in place, and the species method says so rather than leaving a
species that misreports its own symmetry. A cube built with distinct colors is a bare 6-cycle,
and coloring its faces alike afterwards would ask that cycle to report 24.
"""
function _recolor!(ps::ParticleSpecies, sites::AbstractVector{<:BindingSite}, colors)
    # The check needs the new labelling in place to run, so apply first and undo on failure:
    # a species left half-recolored would report a symmetry its graph does not have, which is
    # the very thing being guarded against.
    oldsites, oldlabels = copy(sites), copy(labels(graphrep(ps)))
    _recolor!(sites, graphrep(ps), colors)
    try
        _check_encoding(ps)
    catch err
        err isa ArgumentError || rethrow()
        copy!(sites, oldsites)
        setlabels!(graphrep(ps), oldlabels)
        throw(ArgumentError(
            "this recoloring changes the particle's symmetry by more than its graph encoding " *
            "can express, so it cannot be applied in place; build the species again with the " *
            "new colors instead. ($(err.msg))"
        ))
    end
    return nothing
end

function _recolor!(sites::AbstractVector{<:BindingSite}, g::NautyDiGraph, colors)
    length(colors) == length(sites) || throw(ArgumentError("incorrect number of colors"))
    poses = [s.pose for s in sites]
    gauges = [s.gauge for s in sites]
    orbits = siteorbits(poses, gauges, collect(colors))
    stabs = sitestabilisers(poses, gauges, orbits)

    labs = labels(g)
    for i in eachindex(sites)
        sites[i] = setstab(setcolor(sites[i], colors[i]), stabs[i])
        for v in sites[i].vertices
            labs[v] = orbits[i]
        end
    end
    setlabels!(g, labs)
    return nothing
end

"""
    _check_labelling(ps::ParticleSpecies)

Throw unless `ps`'s labelling is at least as fine as its coloring: sites sharing a label must
share a color.

Labels say which sites the particle cannot tell apart; colors say which bonds a site offers.
A labelling coarser than the coloring asserts both at once — that two sites are the same, and
that they behave differently — and two separate mechanisms then read the wrong one. The graph
records only the label, so structures whose futures differ get merged; and
[`_propagate_faces`](@ref) aligns twist references using the label-preserving rotations, which
for a coarser labelling are not symmetries of the *colored* particle at all, so the bond
registry it settles on is one the coloring never asked for.

Labels derived from colors by [`siteorbits`](@ref) satisfy this by construction, and usually
strictly: two sites are put in one orbit only if they share a color *and* a rotation carries
one onto the other, so geometry splits a color class wherever nothing relates its members. A
labelling made finer still is the one good reason to set labels by hand — it breaks a symmetry
on purpose, which is always safe, since the graph then distinguishes more rather than less.
Only coarser is rejected, and the fix is to say the thing in colors instead.

Being insensitive is not the same as being sound. A species with a coarse labelling and a bond
table of self-bonds `(i, i)` can come out right anyway, since such a bond never depends on the
registry between two *different* faces. Bonds between different faces do, and those break.
"""
function _check_labelling(ps::ParticleSpecies)
    labs = labels(graphrep(ps))
    seen = Dict{eltype(labs),Int}()
    for i in 1:nsites(ps)
        b = bindingsites(ps, i)
        l = labs[first(b.vertices)]
        c = get!(seen, l, color(b))
        c == color(b) || throw(ArgumentError(
            "sites sharing symmetry label $l have colors $c and $(color(b)); a labelling must " *
            "be at least as fine as the coloring, since the graph and the twist references " *
            "follow the label while the bonds follow the color. Drop the labels and give the " *
            "sites the coloring you mean — the labelling is derived from it."
        ))
    end
    return ps
end

"""
    _check_encoding(ps::ParticleSpecies)

Throw if `ps`'s graph claims a different symmetry than its binding sites actually have, or if
its labelling is coarser than its coloring (see [`_check_labelling`](@ref)).

`symmetrynumber(ps) == site_symmetry(ps)` is what makes a graph a correct encoding of a
particle, and it does *not* follow from having used [`cycleencoding`](@ref) or
[`dartencoding`](@ref) — either can be applied to an arrangement it does not describe. Too
large a symmetry number merges structures that are really distinct; too small a one splits
structures that are really the same. Both corrupt enumeration silently, so the species
constructors that build their own graph check here instead.

The symmetry check is one-sided in a way worth relying on: a twist reference chosen
inconsistently across faces can only make `site_symmetry` too *small*, never too large. A
rotation matching every site's position and normal maps every face's supporting plane to
another, hence maps the solid — an intersection of half-spaces, by convexity — onto itself,
so it is a genuine symmetry. Spurious symmetries are therefore impossible, and this check
catches the only direction [`_canonical_faces`](@ref) can fail in.
"""
function _check_encoding(ps::ParticleSpecies)
    _check_labelling(ps)
    geometric = site_symmetry(ps)
    graph = symmetrynumber(ps)
    graph == geometric || throw(ArgumentError(
        "graph encoding claims a symmetry number of $graph, but the binding sites of this " *
        "$(dimension(ps))d species have a rotational symmetry of $geometric. " *
        (graph > geometric ?
         "The labels declare sites equivalent that no rotation maps onto each other." :
         "The graph distinguishes sites that a rotation does map onto each other; the " *
         "encoding does not describe this site arrangement.")
    ))
    return ps
end
