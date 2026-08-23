"""
    BindingSite{P<:Pose}

A `BindingSite` describes an anchor point at which particles are attached together.
The binding site color, together with an interaction matrix, determines whether two binding sites may bind.

`gauge` is the order of the site's own rotational symmetry about its outward normal, disregarding the rest
of the particle.
`stab` size of the stabilizer of this site in the rotational symmetry group of the containing particle.
For example, consider a square side of a 3-prism (i.e. an extruded triangle). The side has 4-fold symmetry
(`gauge == 4`), but only two of these rotations are symmetries of the entire particle (`stab == 2`).

A binding site's pose carries three things. Its position and its outward normal (the pose's local x axis)
are physical. The remaining freedom in the frame, the *twist*, is generally not: 
a site with a `gauge`-fold symmetry is unchanged by turns of `2π/gauge` about its normal,
so `psi`, `psi·Rx(2π/gauge)`, ... all describe the same site.

`locking` determines if the binding site fixes the twist between it and its binding partners. A locking
site holds its partner fixed in the relative orientation defined by its frame (see [`standard_offset`](@ref)), 
up to symmetry transformations of the particle, so the number of the possible twists is equal to `stab`. 
If `locking==false`, the site allows all twists that are compatible with the site symmetry alone, disregarding
the symmetry of the particle, so the number of twists is equal to `gauge`.
The two choices of locking coincide whenever the stabilizer of the site is equal to its individual symmetry group.
"""
struct BindingSite{P<:Pose,F<:Real}
    pose::P
    color::Int
    vertices::UnitRange{Int}
    touching_tolerance::F
    alignment_tolerance::F
    gauge::Int
    stab::Int
    locking::Bool
end
function BindingSite(
    pose::P,
    color::Integer,
    vertices::UnitRange{<:Integer},
    touching_tolerance::Real,
    alignment_tolerance::Real,
    gauge::Integer=1,
    stab::Integer=1,
    locking::Bool=true,
) where {P<:Pose}
    F = eltype(P)
    return BindingSite{P,F}(
        pose, color, vertices, convert(F, touching_tolerance), convert(F, alignment_tolerance), gauge, stab, locking
    )
end

@inline BindingSite(
    site::BindingSite;
    pose=site.pose,
    color=site.color,
    vertices=site.vertices,
    gauge=site.gauge,
    stab=site.stab,
    locking=site.locking,
) = typeof(site)(pose, color, vertices, site.touching_tolerance, site.alignment_tolerance, gauge, stab, locking)

Base.:(==)(b1::BindingSite, b2::BindingSite) = b1.vertices == b2.vertices && b1.color == b2.color
Base.hash(b::BindingSite, h::UInt) = hash(b.color, hash(b.vertices, h))
Base.isless(b1::BindingSite, b2::BindingSite) = isless(b1.vertices, b2.vertices)
Base.:*(p, site::BindingSite) = BindingSite(site; pose=p * site.pose)
Base.:*(site::BindingSite, p) = BindingSite(site; pose=site.pose * p)
@inline shift_vertices(site::BindingSite, v::Integer) = BindingSite(site; vertices=site.vertices .+ v)
@inline shift_color(site::BindingSite, c::Integer) = BindingSite(site; color=site.color + c)
@inline setcolor(site::BindingSite, c::Integer) = BindingSite(site; color=c)
@inline setstab(site::BindingSite, s::Integer) = BindingSite(site, stab=s)
posetype(::BindingSite{P}) where {P} = P
posetype(::Type{<:BindingSite{P}}) where {P} = P

@inline standard_rotation(::Type{F}, ::Val{2}, args...) where {F} = Angle2d{F}(π)
@inline standard_rotation(::Type{F}, ::Val{3}, r::Integer=0, L::Integer=1) where {F} = RotXYZ{F}(2F(π) * r / L, 0, π)

"""
    standard_offset(b::BindingSite, r=0, L=1)
    standard_offset(p::Pose, r=0, L=1)

Calculate the offset between binding site `b` and a partner attached in phase `r` of `L`.
We employ the convention that attached binding sites are "facing each other": their poses are
related via a 180 degree rotation around their (shared) z-axes.

A bond in 3D additionally admits a *twist* about the bond (x-)axis.
This twist is generally fixed by the two sites' frames, except that the frame of a symmetric particle
can only be determined up to symmetry transformations.
If there are `L` admissable twists, then any twists angle that is a multiple of `2π/L` is valid.
"""
@inline standard_offset(b::BindingSite{<:Pose{D,F}}, r::Integer=0, L::Integer=1) where {D,F} =
    b.pose * standard_rotation(F, Val(D), r, L)
@inline standard_offset(p::Pose{D,F}, r::Integer=0, L::Integer=1) where {D,F} = p * standard_rotation(F, Val(D), r, L)

"""
    twistfreedom(b::BindingSite)

Return how many possible twists binding site `b` allows. Equal to `stab` if the site is locking,
and equal to `gauge` otherwise. See [`BindingSite`](@ref) and [`bondperiod`](@ref).
"""
@inline twistfreedom(b::BindingSite) = b.locking ? b.stab : b.gauge

"""
    bondperiod(b1::BindingSite, b2::BindingSite)

Return `L`, the number of twists about the bond axis that a bond between `b1` and `b2` admits.
Equal to the least common multiple of the two sites' [`twistfreedom`](@ref)s. See
[`standard_offset`](@ref).
"""
@inline bondperiod(b1::BindingSite, b2::BindingSite) = lcm(twistfreedom(b1), twistfreedom(b2))

"""
    nphases(bbody::BindingSite, battach::BindingSite)

Return how many of the [`bondperiod`](@ref) twists give *distinct* structures when binding site `battch`
to a binding site on a polyform, `bbbody`.

This count excludes twists that are differing by a symmetry of the particle being attached, of which
there are exactly `b2.stab`.
Note the asymmetry: only the *attached* binding site's stabilizer is quotiented out. The symmetries of the
polyform carrying `bbody` can merge structures too, but these are not knowable from one site,
and need to be caught through canonization.
"""
@inline nphases(bbody::BindingSite, battach::BindingSite) = bondperiod(bbody, battach) ÷ battach.stab

"""
    contact_pairing(vs1::UnitRange, vs2::UnitRange, r=0, L=1)

Return the pairs of graph vertices that should be joined when two bonding binding sites occupying the vertex
ranges `vs1` and `vs2` meet in phase `r` of `L` (see [`standard_offset`](@ref)).

Imagine the vertices corresponds to corners of a regular n-gon attached to the site, number the vertices of each site
starting from zero and denote `kᵢ = length(vsᵢ)`. Vertex `a` of site 1 then sits at angle `2πa/k₁` in site 1's frame,
and vertex `b` of site 2  (whose frame is turned by `2πr/L` and then flipped, so its cyclic order runs the other way 
around) sits at `2πr/L - 2πb/k₂`. Writing `K = lcm(k₁, k₂)`, `s₁ = K/k₁`, `s₂ = K/k₂`, the coincidences are the
solutions of

    a*s₁ + b*s₂ = t   (mod K),      t = r*K/L

`t` is always an integer: a site's gauge `gᵢ` divides its vertex count, so `L = lcm(g₁, g₂)`
divides `K`. Since `gcd(s₁, s₂) = 1`, the congruence has exactly `gcd(k₁, k₂)` solutions for
every `t`, and distinct phases give disjoint solution sets, so the pairing count never
depends on the phase, and the graph records which phase a bond is in.
"""
@inline function contact_pairing(vs1::UnitRange{<:Integer}, vs2::UnitRange{<:Integer}, r::Integer=0, L::Integer=1)
    k1, k2 = length(vs1), length(vs2)
    G = gcd(k1, k2)
    K = k1 * k2 ÷ G
    a0, b0 = _phase_shift(k1, k2, K, r, L)
    return (vs1[1 + mod(a0 + j * (k1 ÷ G), k1)] => vs2[1 + mod(b0 - j * (k2 ÷ G), k2)] for j in 0:(G - 1))
end

# One solution of a*s1 + b*s2 = r*K/L (mod K); the rest follow by stepping the two ranges in
# opposite directions. Zero-based, and `(0, 0)` for the base phase.
@inline function _phase_shift(k1::Int, k2::Int, K::Int, r::Integer, L::Integer)
    r == 0 && return 0, 0
    K % L == 0 || throw(
        ArgumentError(
            "a bond in $L phases cannot be encoded by sites of $k1 and $k2 graph vertices: " *
            "$L must divide lcm($k1, $k2) = $K. Give the sites more vertices, i.e. a graph built " *
            "by `dartencoding` rather than `cycleencoding`, or declare a smaller gauge.",
        ),
    )
    s1, s2 = K ÷ k1, K ÷ k2
    t = r * K ÷ L
    # u*s1 = 1 (mod s2), so a0 = t*u solves the congruence mod s2 and hence mod K once b0
    # absorbs the remainder. gcd(s1, s2) = 1 always, so the inverse exists.
    a0 = s2 == 1 ? 0 : mod(t * invmod(s1, s2), k1)
    return a0, mod((t - a0 * s1) ÷ s2, k2)
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
    return isapprox(b1.pose.x, b2.pose.x; atol=b1.touching_tolerance + b2.touching_tolerance, rtol=0)
end

"""
    phase(b1::BindingSite, b2::BindingSite)

Return which of the bond's phases the two sites are in, or `nothing` if their
orientations are not related by any of them.

This function is symmetric in its arguments.
"""
function phase(b1::BindingSite, b2::BindingSite)
    atol = b1.alignment_tolerance + b2.alignment_tolerance
    L = bondperiod(b1, b2)
    for r in 0:(L - 1)
        isapprox(b1.pose.psi, standard_offset(b2, r, L).psi; atol, rtol=0) && return r
    end
    return nothing
end

"""
    isaligned(b1::BindingSite, b2::BindingSite)

Check whether the orientation components of the binding sites' poses
differ by a standard offset, in any phase the bond admits.
"""
isaligned(b1::BindingSite, b2::BindingSite) = !isnothing(phase(b1, b2))

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