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
    PolyhedronParticleSpecies(p::Polyhedron; colors=1:nfaces(p), labels=colors, encoding=:auto)

Build a particle from the polyhedron `p`, with one binding site at the centroid of each face.

`colors` assigns interaction colors to the binding sites, `labels` the symmetry labels that
determine graph isomorphism, exactly as for [`PolygonParticleSpecies`](@ref). Faces sharing a
label are treated as equivalent, which sets the symmetry number: pass
`labels=geometriclabels(p)` to recover the solid's full rotation group, or leave the default
for all faces distinct.

With `encoding=:auto` the graph encoding is chosen from the labels. All labels distinct means
no rotation can preserve them, so the cheap [`cycleencoding`](@ref) is provably equivalent and
is used; any repeated label needs the full [`dartencoding`](@ref), which is what carries the
rotation group. For a cube that is 6 vertices versus 24. Pass `encoding=:dart` or
`encoding=:cycle` to override.

Binding site frames follow the package convention that the outward normal is the site's local
x axis. The local z axis is fixed to point from the face centroid at the midpoint of the
face's first edge, which is what makes the twist of a face-to-face bond well defined: since
the standard bond rotation preserves local z, two bonded faces meet with their first edges
coincident.

!!! note "One registry per face pair"
    Two faces bond in a single relative twist (see [`standard_offset`](@ref)), not in all `k`
    of them. Which twist that is follows from the site frames, so it is a modelling choice,
    and `Polyhedron` makes the one that lets a solid assemble into its own tiling:

    - faces related by a **translation** are given the same local z, so the bond between them
      is a pure translation and loops of such bonds compose to the identity;
    - faces with no translation mate get their local z along the solid's **principal axis**,
      so their bonds are π rotations about one shared axis and loops of even length compose
      to the identity.

    Between them these cover every space-filling solid in the library — cubes enumerate as
    polycubes (A000162), triangular prisms as polyiamonds (A000577), hexagonal prisms as
    polyhexes (A000228). A solid that tiles only by rotations about *several* axes, and a
    face pair that is neither of the above, still get one twist out of `k`; to use another,
    add a species whose site poses are rotated accordingly.
"""
function PolyhedronParticleSpecies(
    p::Polyhedron{F}; colors=1:nfaces(p), labels=colors, encoding::Symbol=:auto
) where {F}
    n = nfaces(p)
    length(colors) == n ||
        throw(ArgumentError("expected $n colors, one per face, got $(length(colors))"))
    encoding in (:auto, :dart, :cycle) ||
        throw(ArgumentError("encoding must be :auto, :dart or :cycle, got :$encoding"))

    usecycle = encoding === :cycle || (encoding === :auto && allunique(labels))
    g, ranges = usecycle ? cycleencoding(n; labels) : dartencoding(p; labels)

    rmin = inradius(p)
    rmax = bounding_radius(p)
    tol = sqrt(eps(F)) * rmax

    P = Pose{3,F,RotMatrix3{F}}
    sites = Vector{BindingSite{P,F}}(undef, n)
    for i in 1:n
        x = facecentroid(p, i)
        ex = facenormal(p, i)
        ez = normalize(edgemidpoint(p, i, 1) - x)
        psi = RotMatrix3{F}(hcat(ex, cross(ez, ex), ez))
        sites[i] = BindingSite(P(x, psi), colors[i], ranges[i], tol, tol / rmin)
    end

    return _check_encoding(PolyhedronParticleSpecies{F,eltype(sites)}(
        g, sites, p, facenormals(p), _edgedirections(p), rmin, rmax, tol
    ))
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
    length(colors) != nsites(p) && throw(ArgumentError("incorrect number of colors"))
    for k in eachindex(p.sites)
        s = p.sites[k]
        p.sites[k] = BindingSite(s.pose, colors[k], s.vertices, s.touching_tolerance, s.alignment_tolerance)
    end
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
    # Cheap bounds first: the bounding spheres cannot reach, or the inscribed spheres must
    # already intersect.
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
