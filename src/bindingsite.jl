"""
    BindingSite{P<:Pose}

A `BindingSite` describes an anchor point at which particles can be attached together.

The binding site color, together with an interaction matrix, determines whether two
binding sites may bind. One binding site may be represented by multiple contiguous
graph vertices and particle/polyform graph representation.
"""
struct BindingSite{P<:Pose, F<:Real}
    pose::P
    color::Int
    vertices::UnitRange{Int}
    touching_tolerance::F
    alignment_tolerance::F
end
function BindingSite(pose::P, color::Integer, vertices::UnitRange{<:Integer},
                     touching_tolerance::Real, alignment_tolerance::Real) where {P<:Pose}
    F = eltype(P)
    return BindingSite{P,F}(pose, color, vertices, convert(F, touching_tolerance), convert(F, alignment_tolerance))
end

Base.:(==)(b1::BindingSite, b2::BindingSite) = b1.vertices == b2.vertices && b1.color == b2.color
Base.hash(b::BindingSite, h::UInt) = hash(b.color, hash(b.vertices, h))
Base.isless(b1::BindingSite, b2::BindingSite) = isless(b1.vertices, b2.vertices)
Base.:*(p, site::BindingSite) = typeof(site)(p * site.pose, site.color, site.vertices, site.touching_tolerance, site.alignment_tolerance)
Base.:*(site::BindingSite, p) = typeof(site)(site.pose * p, site.color, site.vertices, site.touching_tolerance, site.alignment_tolerance)
@inline shift_vertices(site::BindingSite, v::Integer) = typeof(site)(site.pose, site.color, site.vertices .+ v, site.touching_tolerance, site.alignment_tolerance)
@inline shift_color(site::BindingSite, c::Integer) = typeof(site)(site.pose, site.color + c, site.vertices, site.touching_tolerance, site.alignment_tolerance)
posetype(::BindingSite{P}) where P = P
posetype(::Type{<:BindingSite{P}}) where P = P

"""
    standard_offset(b::BindingSite{<:Pose{D,F}}) where {D,F}
    standard_offset(p::Pose{D,F})
    standard_rotation(::Type{F}, ::Val{D})

Calculate the default offset between binding site `b` and an attached partner binding site. We employ the
convention that attached binding sites are "facing each other": their poses are related
via a 180 degree rotation around their (shared) z-axes. `standard_rotation` is that rotation on its own.
"""
@inline standard_rotation(::Type{F}, ::Val{2}) where {F} = Angle2d{F}(π)
@inline standard_rotation(::Type{F}, ::Val{3}) where {F} = RotXYZ{F}(0, 0, π)
@inline standard_offset(b::BindingSite{<:Pose{D,F}}) where {D,F} = b.pose * standard_rotation(F, Val(D))
@inline standard_offset(p::Pose{D,F}) where {D,F} = p * standard_rotation(F, Val(D))

"""
    contact_pairing(vs1::UnitRange, vs2::UnitRange)

Return the pairs of graph vertices to join when two bonding binding sites occupy the vertex
ranges `vs1` and `vs2`.

Gluing two faces together reverses their orientation: seen from outside the first particle,
the second particle's face runs the other way round. So the vertices are matched
*counter-rotating*, with the first vertex of each range as the fixed point — that is the
vertex pair whose polyhedron edges coincide, by the convention that a binding site's local
z axis points at the midpoint of its face's first edge.

Sites of *different* sizes are joined only where their vertices land at the same angle about
the bond axis: `gcd(k1, k2)` pairs, stepping `k1 ÷ gcd` through one range and `k2 ÷ gcd`
backwards through the other. That is what preserves the bond's residual symmetry. A site
with `k` vertices is invariant under turns by `2π/k`, so the bond is invariant under the
turns common to both — `C_k1 ∩ C_k2 = C_gcd(k1,k2)`. Binding a 6-fold site to a 3-fold one
keeps 3-fold symmetry, and the three links are what let nauty see it.

For sites of one or two vertices this is the identity pairing, so 2D species and patchy
particles are unaffected. This is the single place where the bond convention is defined;
a species needing a different registry would change it here.
"""
@inline function contact_pairing(vs1::UnitRange{Int}, vs2::UnitRange{Int})
    k1, k2 = length(vs1), length(vs2)
    g = gcd(k1, k2)
    step1, step2 = k1 ÷ g, k2 ÷ g
    return (vs1[1 + m * step1] => vs2[mod1(1 - m * step2, k2)] for m in 0:(g - 1))
end

"""
    color(b::BindingSite)

Return the binding site's color.

If the binding site comes from an `BindingRules`, its interactions with other
binding sites are determined by `interactionmatrix(sys)[:, color(b)]`.
"""
color(b::BindingSite) = b.color

"""
    istouching(b1::BindingSite, b2::BindingSite)

Check whether the translational components of the binding sites' poses
are approximately equal.
"""
function istouching(b1::BindingSite, b2::BindingSite)
    return isapprox(b1.pose.x, b2.pose.x;
                    atol=b1.touching_tolerance + b2.touching_tolerance, rtol=0)
end

"""
    isaligned(b1::BindingSite, b2::BindingSite)

Check whether the orientation components of the binding sites' poses
differ by a standard offset.
"""
function isaligned(b1::BindingSite, b2::BindingSite)
    return isapprox(b1.pose.psi, standard_offset(b2).psi;
                    atol=b1.alignment_tolerance + b2.alignment_tolerance, rtol = 0)
end

"""
    isincontact(b1::BindingSite, b2::BindingSite)

Check whether the binding sites are touching and aligned.
"""
function isincontact(b1::BindingSite, b2::BindingSite)
    return istouching(b1, b2) && isaligned(b1, b2)
end

function Base.show(io::Core.IO, b::BindingSite)
    print(io, "BindingSite[c=$(b.color), vs=$(b.vertices)]")
end
function Base.show(io::Core.IO, ::MIME"text/plain", b::BindingSite)
    println(io, "$(dimension(b.pose))-dimensional BindingSite:")
    print(io, " - color: \t")
    println(io, b.color)
    print(io, " - vertices:\t$(b.vertices)\n")
    print(io, " - pose: \t$(b.pose)")
end