"""
    BindingSite{P<:Pose}

A `BindingSite` describes an anchor point at which particles can be attached together.

The binding site color, together with an interaction matrix, determines whether two
binding sites may bind. One binding site may be represented by multiple contiguous
graph vertices and particle/polyform graph representation. The relative poses of two
attached binding sites is given by `standard_offset`.
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

"""
    standard_offset(b::BindingSite{<:Pose{D,F}}) where {D,F}

Calculate the default offset between two attached binding sites. We employ the 
convention that attached binding sites are "facing each other": their poses are related 
via a 180 degree rotation around their (shared) z-axes.
"""
@inline standard_offset(b::BindingSite{<:Pose{2,F}}) where {F} = b.pose * Angle2d{F}(π)
@inline standard_offset(b::BindingSite{<:Pose{3,F}}) where {F} =  b.pose * RotXYZ{F}(0, 0, π)

"""
    color(b::BindingSite)

Return the binding site's color.

If the binding site comes from an `AssemblySystem`, its interactions with other
binding sites are determined by `interactionmatrix(sys)[:, color(b)]`. A color
of `0` indicates that the site is inert, i.e. it does not bind to any other site.
"""
color(b::BindingSite) = b.color

"""
    istouching(b1::BindingSite, b2::BindingSite; kwargs...)

Check whether the translational components of the binding sites' poses
are approximately equal. 

Absolute (`atol`) and relative (`rtol`) tolerances should be supplied
as keyword arguments.
"""
function istouching(b1::BindingSite, b2::BindingSite)
    return isapprox(b1.pose.x, b2.pose.x;
                    atol = b1.touching_tolerance + b2.touching_tolerance, rtol = 0)
end

"""
    isaligned(b1::BindingSite, b2::BindingSite; kwargs...)

Check whether the orientation components of the binding sites' poses
differ by a standard offset.

Absolute (`atol`) and relative (`rtol`) tolerances should be supplied
as keyword arguments.
"""
function isaligned(b1::BindingSite, b2::BindingSite)
    return isapprox(b1.pose.psi, standard_offset(b2).psi;
                    atol = b1.alignment_tolerance + b2.alignment_tolerance, rtol = 0)
end

"""
    isincontact(b1::BindingSite, b2::BindingSite; kwargs...)

Check whether the binding sites' are touching and aligned.

Absolute (`atol`) and relative (`rtol`) tolerances should be supplied
as keyword arguments.
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