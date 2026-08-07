"""
    Polyhedron{F}

A convex polyhedron, stored as a list of `corners` and a list of `faces`.

Each face is a list of indices into `corners`, wound counter-clockwise as seen from
*outside* the solid. The winding is what fixes the orientation of the graph encoding, and
is checked on construction.
"""
struct Polyhedron{F<:AbstractFloat}
    corners::Vector{SVector{3,F}}
    faces::Vector{Vector{Int}}

    function Polyhedron(corners::Vector{SVector{3,F}}, faces::Vector{Vector{Int}}) where {F}
        _check_winding(length(corners), faces)
        # Twist references, weakest rule first so the stronger one can override it.
        return new{F}(corners, _align_mates(corners, _align_axis(corners, faces)))
    end
end

"""
    Polyhedron(corners, faces)
    Polyhedron(corners)

Construct a polyhedron from its `corners` and `faces`. Faces must be wound counter-clockwise
as seen from outside the solid; this is verified, since a reversed face would silently
produce a graph encoding with the wrong symmetry group.

If `faces` is omitted, they are derived from the corners by finding all supporting planes of
the convex hull, which requires the solid to be convex.
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
"""
function _check_winding(ncorners::Integer, faces::AbstractVector{<:AbstractVector{Int}})
    length(faces) < 3 && throw(ArgumentError("a polyhedron needs at least 3 faces"))
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
    _principal_axis(corners)

Return the solid's unique principal axis, or `nothing` if it has no distinguished one.

Taken from the inertia tensor of the corner set: a solid with a single rotation axis of order
3 or more (a prism, antiprism or pyramid) has one distinct eigenvalue and a degenerate pair,
whereas the cube and the other high-symmetry solids have three equal eigenvalues and no
preferred direction.
"""
function _principal_axis(corners::Vector{SVector{3,F}}) where {F}
    c0 = sum(corners) / length(corners)
    tensor = zero(SMatrix{3,3,F})
    for c in corners
        r = c - c0
        tensor += dot(r, r) * one(SMatrix{3,3,F}) - r * transpose(r)
    end
    vals, vecs = eigen(Symmetric(Matrix(tensor)))
    atol = sqrt(eps(F)) * maximum(abs, vals)
    for k in 1:3
        i, j = filter(!=(k), 1:3)
        isapprox(vals[i], vals[j]; atol) && !isapprox(vals[k], vals[i]; atol) &&
            return SVector{3,F}(view(vecs, :, k))
    end
    return nothing
end

"""
    _align_axis(corners, faces)

Point every face's first edge as nearly along the solid's principal axis as it can, so that
all the site frames share one local z direction.

This is the twist convention for faces with no translation mate — the side faces of an
odd prism, say, no two of which are parallel. A bond between two such faces turns the
attached particle by π about the site's local z (that is what [`standard_rotation`](@ref)
does), so making those axes *parallel* means the turns all commute and a loop of even length
composes back to the identity. Triangular prisms then close the six-triangle ring and
enumerate as polyiamonds; with the local z axes left unrelated, that ring is lost.

A no-op when the solid has no distinguished axis, and overridden by [`_align_mates`](@ref)
wherever a translation mate exists, which is the stronger convention.
"""
function _align_axis(corners::Vector{SVector{3,F}}, faces::Vector{Vector{Int}}) where {F}
    axis = _principal_axis(corners)
    isnothing(axis) && return faces

    faces = [copy(f) for f in faces]
    for (i, f) in enumerate(faces)
        c = sum(corners[v] for v in f) / length(f)
        alignment(k) = abs(dot(normalize((corners[f[k]] + corners[f[mod1(k + 1, length(f))]]) / 2 - c), axis))
        k = argmax(alignment, eachindex(f))
        faces[i] = circshift(f, 1 - k)
    end
    return faces
end

"""
    _align_mates(corners, faces)

Cyclically rotate face vertex lists so that *mated* faces — faces related by a pure
translation, which is what lets a solid tile space without turning — start on the same edge.

This is what makes a bond between mated faces a pure translation rather than a turn. The
relative orientation of a bond is `sᵢ · Δ · sⱼ⁻¹`, where `sᵢ`, `sⱼ` are the two binding site
frames and `Δ` is [`standard_rotation`](@ref); since `Δ` is a π rotation about the frames'
local z, and the local z of a site points at its face's *first* edge, aligning the first
edges of a mated pair gives `sⱼ = sᵢ·Δ` and hence a relative orientation of the identity.

The payoff is that transport around any closed loop of mated bonds composes to the identity,
so ring closures always succeed and a space-filling solid assembles into its lattice. Faces
without a translation mate — every face of a tetrahedron, say — are left alone; such solids
do not tile by translation and have no orientation-free convention to find.

Idempotent, and a no-op for solids with no mated faces.
"""
function _align_mates(corners::Vector{SVector{3,F}}, faces::Vector{Vector{Int}}) where {F}
    faces = [copy(f) for f in faces]
    n = length(faces)
    atol = sqrt(eps(F)) * maximum(norm, corners)

    centroid(f) = sum(corners[v] for v in f) / length(f)
    midpoint(f, k) = (corners[f[k]] + corners[f[mod1(k + 1, length(f))]]) / 2

    paired = falses(n)
    for i in 1:n
        paired[i] && continue
        ci = centroid(faces[i])
        # A mate carries the same corner set, shifted by the translation between centroids.
        j = findfirst(1:n) do k
            k == i && return false
            paired[k] && return false
            length(faces[k]) == length(faces[i]) || return false
            t = centroid(faces[k]) - ci
            all(v -> any(w -> isapprox(corners[v] + t, corners[w]; atol), faces[k]), faces[i])
        end
        isnothing(j) && continue

        target = midpoint(faces[i], 1) + (centroid(faces[j]) - ci)
        k = findfirst(l -> isapprox(midpoint(faces[j], l), target; atol), eachindex(faces[j]))
        isnothing(k) && continue

        faces[j] = circshift(faces[j], 1 - k)
        paired[i] = paired[j] = true
    end
    return faces
end

"""
    _derive_faces(corners)

Find the faces of the convex hull of `corners` by testing every corner triple for a
supporting plane, then winding each face counter-clockwise about its outward normal.
"""
function _derive_faces(corners::Vector{SVector{3,F}}) where {F}
    n = length(corners)
    center = sum(corners) / n
    scale = maximum(norm(c - center) for c in corners)
    atol = sqrt(eps(F)) * scale

    faces = Vector{Int}[]
    planes = Tuple{SVector{3,F},F}[]
    for i in 1:(n - 2), j in (i + 1):(n - 1), k in (j + 1):n
        nrm = cross(corners[j] - corners[i], corners[k] - corners[i])
        norm(nrm) < atol && continue
        nrm = normalize(nrm)
        d = dot(nrm, corners[i])
        # Point the normal away from the interior.
        dot(nrm, center) > d && ((nrm, d) = (-nrm, -d))
        # Supporting plane of the hull?
        all(c -> dot(nrm, c) <= d + atol, corners) || continue
        any(pl -> isapprox(pl[1], nrm; atol) && isapprox(pl[2], d; atol), planes) && continue

        push!(planes, (nrm, d))
        push!(faces, _wind_ccw(corners, findall(c -> abs(dot(nrm, c) - d) <= atol, corners), nrm))
    end
    return faces
end

function _wind_ccw(corners::Vector{SVector{3,F}}, idxs::Vector{Int}, nrm::SVector{3,F}) where {F}
    c = sum(corners[i] for i in idxs) / length(idxs)
    u = normalize(corners[first(idxs)] - c)
    v = cross(nrm, u)
    return sort(idxs; by=i -> atan(dot(corners[i] - c, v), dot(corners[i] - c, u)))
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

Return the midpoint of the `k`th edge of the `i`th face of `p`, i.e. the edge running from
corner `k` to corner `k+1` of that face. The `k = 1` midpoint is the twist reference of the
face's binding site.
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
    edgelength(p::Polyhedron)

Return the length of the shortest edge of `p`.
"""
function edgelength(p::Polyhedron)
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
function rotationgroup(p::Polyhedron{F}) where {F}
    cs = corners(p)
    c0 = sum(cs) / length(cs)
    cs = [c - c0 for c in cs]
    atol = sqrt(eps(F)) * maximum(norm, cs)

    # An orthonormal frame attached to dart k of face i: along the edge, along the outward
    # normal, and their cross product.
    function dartframe(i, k)
        f = facevertices(p, i)
        e1 = normalize(cs[f[mod1(k + 1, length(f))]] - cs[f[k]])
        e2 = facenormal(p, i)
        return hcat(e1, e2, cross(e1, e2))
    end

    Mref = dartframe(1, 1)
    group = RotMatrix3{F}[]
    for i in 1:nfaces(p), k in 1:facedegree(p, i)
        R = dartframe(i, k) * transpose(Mref)
        all(c -> any(c2 -> isapprox(R * c, c2; atol), cs), cs) || continue
        push!(group, RotMatrix3{F}(R))
    end
    return group
end

"""
    geometriclabels(p::Polyhedron)

Return one label per face, grouping faces that are equivalent under [`rotationgroup`](@ref).

Passing these to [`dartencoding`](@ref) yields the solid's true rotation group, whereas
`labels=fill(1, nfaces(p))` yields the *combinatorial* symmetry of the face lattice, which
can be larger: `Pyramid(3)` is combinatorially a tetrahedron and would report 12 instead of
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

################################################################################
# Graph encodings
################################################################################

"""
    dartencoding(p::Polyhedron; labels=1:nfaces(p))

Return `(g, ranges)`, the dart encoding of `p` and the graph vertices belonging to each face.

A *dart* is one corner of one face, so a face of degree `k` owns `k` darts and the graph has
`2 * nedges(p)` vertices in total. The encoding consists of

1. a directed `k`-cycle through each face's own darts, following the face's counter-clockwise
   winding, and
2. a bidirectional edge joining the two darts that sit on each shared polyhedron edge.

The automorphism group of the result is the solid's proper rotation group — the directed face
cycles are what exclude reflections, exactly as the directed cycle does for a polygon in 2D.
`labels` assigns one symmetry label per face, inherited by that face's darts: all labels
distinct gives a symmetry number of 1, all labels equal gives the full rotation group, and
merging some faces gives the subgroup preserving that labelling.

Only `faces(p)` is used; the corner positions play no role.
"""
function dartencoding(p::Polyhedron; labels=1:nfaces(p))
    fs = faces(p)
    length(labels) == length(fs) ||
        throw(ArgumentError("expected $(length(fs)) labels, one per face, got $(length(labels))"))

    ranges = Vector{UnitRange{Int}}(undef, length(fs))
    o = 0
    for (i, f) in enumerate(fs)
        ranges[i] = (o + 1):(o + length(f))
        o += length(f)
    end

    vertex_labels = Cint[labels[i] for (i, f) in enumerate(fs) for _ in f]
    g = NautyDiGraph(o; vertex_labels)

    # Dart k of face i sits on the directed edge f[k] -> f[k+1]. Its partner is the dart of
    # the adjacent face carrying the reversed edge; _check_winding guarantees it exists.
    dart_of_edge = Dict{Tuple{Int,Int},Int}()
    for (i, f) in enumerate(fs), k in eachindex(f)
        dart_of_edge[(f[k], f[mod1(k + 1, length(f))])] = first(ranges[i]) + k - 1
    end

    for (i, f) in enumerate(fs)
        d0 = first(ranges[i])
        for k in eachindex(f)
            d = d0 + k - 1
            add_edge!(g, d, d0 + mod(k, length(f)))
            partner = dart_of_edge[(f[mod1(k + 1, length(f))], f[k])]
            if d < partner
                add_edge!(g, d, partner)
                add_edge!(g, partner, d)
            end
        end
    end
    return g, ranges
end

"""
    cycleencoding(nsites; labels=1:nsites)

Return `(g, ranges)`, the sparse encoding of a particle with `nsites` binding sites: a single
directed cycle carrying one vertex per site.

This is the 2D polygon encoding, and it is also valid in 3D whenever all `labels` are
distinct — no rotation can then preserve the labelling, so the symmetry number is 1 whatever
structure the graph has internally, and it is far cheaper than [`dartencoding`](@ref) (6
vertices instead of 24 for a cube). It is *not* valid when labels are shared: an all-equal
6-cycle has a symmetry number of 6, where a cube needs 24.

One vertex per site for every `nsites >= 1`. Two sites give a 2-cycle, which is a pair of
opposite arcs — that used to be barred, because bonds were told apart from a particle's own
edges by being bidirectional, and a 2-cycle inside a particle would have been misread as a
bond. Bonds are now recognised by joining different particles, so the workaround of padding
fewer than three sites out to a 4-cycle is no longer needed. It was also wrong: two
equivalent sites padded to a 4-cycle report a symmetry number of 4, where the particle's
symmetry is 2.
"""
function cycleencoding(nsites::Integer; labels=1:nsites)
    length(labels) == nsites ||
        throw(ArgumentError("expected $nsites labels, one per site, got $(length(labels))"))
    nsites < 1 && throw(ArgumentError("a particle needs at least one binding site"))

    g = NautyDiGraph(cycle_digraph(nsites); vertex_labels=collect(Cint, labels))
    return g, [i:i for i in 1:nsites]
end

"""
    Polyhedron(sym::Symbol, n=0; a=1.0)

Return a solid realizing the proper rotation group named by `sym`, with edge length `a`:

| `sym` | group   | solid           | order |
|:------|:--------|:----------------|:------|
| `:C`  | `C_n`   | `Pyramid(n)`    | `n`   |
| `:D`  | `D_n`   | `Prism(n)`      | `2n`  |
| `:T`  | `T`     | `Tetrahedron()` | 12    |
| `:O`  | `O`     | `Cube()`        | 24    |
| `:I`  | `I`     | `Dodecahedron()`| 60    |

`:C` and `:D` require `n`. The alternative realizations `Octahedron()`, `Icosahedron()` and
`Antiprism(n)` have the same groups and are constructed by name.
"""
function Polyhedron(sym::Symbol, n::Integer=0; a=1.0)
    if sym in (:C, :D)
        n >= 3 || throw(ArgumentError("$sym requires n >= 3, got n=$n"))
        return sym === :C ? Pyramid(n, a) : Prism(n, a)
    end
    n == 0 || throw(ArgumentError("$sym takes no n"))
    sym === :T && return Tetrahedron(a)
    sym === :O && return Cube(a)
    sym === :I && return Dodecahedron(a)
    return throw(ArgumentError("unknown rotation group $sym, expected one of :C, :D, :T, :O, :I"))
end

"""
    sitegraph(sym::Symbol, n=0; labels)

Shorthand for `dartencoding(Polyhedron(sym, n); labels)`: the graph encoding of the proper
rotation group named by `sym`. See [`Polyhedron`](@ref) for the solid each symbol resolves to.
"""
function sitegraph(sym::Symbol, n::Integer=0; kwargs...)
    return dartencoding(Polyhedron(sym, n); kwargs...)
end

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
