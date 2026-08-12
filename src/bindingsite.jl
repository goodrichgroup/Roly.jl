"""
    BindingSite{P<:Pose}

A `BindingSite` describes an anchor point at which particles can be attached together.

The binding site color, together with an interaction matrix, determines whether two
binding sites may bind. One binding site may be represented by multiple contiguous
graph vertices and particle/polyform graph representation.

`gauge` is the order of the site's *own* rotational symmetry about its outward normal, and
`stab` how much of that the whole particle keeps; see the note below. Both are fixed once the
particle's labelling is, so they are stored rather than recomputed per attachment attempt.

!!! note "What a binding site's orientation means"
    A site's pose carries three things, and they do not have the same status. Its position and
    its outward normal — the pose's local x axis — are physical. The remaining freedom in the
    frame, the *twist reference*, is not: a site with a `gauge`-fold symmetry is unchanged by
    turns of `2π/gauge` about its normal, so `psi`, `psi·Rx(2π/gauge)`, … all describe the same
    site. A square polyhedron face is such a site; a patch is one when the model resolves it
    into several graph vertices.

    Two counts meet here and are worth keeping apart. `gauge` is the *site's* symmetry. How
    much of it the whole particle keeps — the site's stabiliser in the particle's rotation
    group, `stab` — is a different, smaller number. A triangular prism is only 2-fold about a
    square side face, yet the face is still a square: `gauge` 4, `stab` 2. They are used for
    different things. `gauge` is what [`site_symmetry`](@ref) needs, since asking whether a
    rotation carries face *i* onto face *j* must not depend on which edge was called first, and
    the stabiliser would reject genuine symmetries — a rotation can carry one face's first dart
    onto another face's second. `stab` is what [`nregistrations`](@ref) needs, to skip the
    attachments that differ only by a symmetry of the particle being attached.

    `locking` picks between them for the purpose of deciding what a bond here means. A locking
    site — the default — holds its partner in the one relative orientation its frame names, and
    only its particle's own symmetry can loosen that, so its say in the bond is `stab`. A
    rotation-free site declines to fix the twist at all and admits every orientation its face
    geometrically permits, so its say is `gauge`. The two coincide whenever a face is no more
    symmetric than the body around it, and always in 2D, where a site has no turn about its
    in-plane normal to argue about.

    The twist reference becomes observable only through [`standard_offset`](@ref), which of the
    orientations two sites could meet in recognises exactly one — fixed by the two sites' twist
    references *relative to each other*. The individual choices mean nothing; the relative
    choices across a particle's sites are its bond registry, and decide how orientation is
    transported along a chain of bonds, hence which rings close. That is a modelling choice,
    not a fact about the shape, and no symmetry check can validate it. Nor is there a canonical
    choice to be had: for the Platonic solids the rotation group acts freely and transitively
    on darts, so the orbit of any one dart is *every* dart and no symmetric assignment of one
    per face exists. Every convention breaks the symmetry; the only question is whether it
    breaks it the way the intended lattice needs.
"""
struct BindingSite{P<:Pose, F<:Real}
    pose::P
    color::Int
    vertices::UnitRange{Int}
    touching_tolerance::F
    alignment_tolerance::F
    gauge::Int
    stab::Int
    locking::Bool
end
function BindingSite(pose::P, color::Integer, vertices::UnitRange{<:Integer},
                     touching_tolerance::Real, alignment_tolerance::Real,
                     gauge::Integer=1, stab::Integer=1, locking::Bool=true) where {P<:Pose}
    F = eltype(P)
    return BindingSite{P,F}(pose, color, vertices, convert(F, touching_tolerance),
                            convert(F, alignment_tolerance), gauge, stab, locking)
end

# Rebuild a site changing only the three fields anything ever changes, so the symmetry counts
# and tolerances travel along untouched.
@inline _respecify(site::BindingSite, pose, color, vertices) =
    typeof(site)(pose, color, vertices, site.touching_tolerance, site.alignment_tolerance,
                 site.gauge, site.stab, site.locking)

Base.:(==)(b1::BindingSite, b2::BindingSite) = b1.vertices == b2.vertices && b1.color == b2.color
Base.hash(b::BindingSite, h::UInt) = hash(b.color, hash(b.vertices, h))
Base.isless(b1::BindingSite, b2::BindingSite) = isless(b1.vertices, b2.vertices)
Base.:*(p, site::BindingSite) = _respecify(site, p * site.pose, site.color, site.vertices)
Base.:*(site::BindingSite, p) = _respecify(site, site.pose * p, site.color, site.vertices)
@inline shift_vertices(site::BindingSite, v::Integer) = _respecify(site, site.pose, site.color, site.vertices .+ v)
@inline shift_color(site::BindingSite, c::Integer) = _respecify(site, site.pose, site.color + c, site.vertices)
@inline setcolor(site::BindingSite, c::Integer) = _respecify(site, site.pose, c, site.vertices)
@inline setstab(site::BindingSite, s::Integer) =
    typeof(site)(site.pose, site.color, site.vertices, site.touching_tolerance,
                 site.alignment_tolerance, site.gauge, s, site.locking)
posetype(::BindingSite{P}) where P = P
posetype(::Type{<:BindingSite{P}}) where P = P

# The face-to-face rotation on its own, without a pose to apply it to; see `standard_offset`.
@inline standard_rotation(::Type{F}, ::Val{2}, r::Integer=0, L::Integer=1) where {F} = Angle2d{F}(π)
@inline standard_rotation(::Type{F}, ::Val{3}, r::Integer=0, L::Integer=1) where {F} =
    RotXYZ{F}(2F(π) * r / L, 0, π)

"""
    standard_offset(b::BindingSite{<:Pose{D,F}}, r=0, L=1) where {D,F}
    standard_offset(p::Pose{D,F}, r=0, L=1)

Calculate the offset between binding site `b` and a partner attached in registration `r` of `L`.
We employ the convention that attached binding sites are "facing each other": their poses are
related via a 180 degree rotation around their (shared) z-axes. `standard_rotation` is that
rotation on its own.

A bond additionally admits a *twist* about the bond axis — the site's local x, its outward
normal. Which twist is meant is fixed by the two sites' frames, and would be a single number,
except that a symmetric particle has no single frame to offer: when a rotation carries the
particle onto itself while turning a site about its normal, the two frames it relates describe
that particle equally well and neither can claim the bond. So the admissible twists are what
survives regauging either site, and form the `L` multiples of `2π/L` that
[`bondperiod`](@ref) counts. `r = 0` is the frames' own choice, and is the only registration
when neither site has anything to be ambiguous about.
"""
@inline standard_offset(b::BindingSite{<:Pose{D,F}}, r::Integer=0, L::Integer=1) where {D,F} =
    b.pose * standard_rotation(F, Val(D), r, L)
@inline standard_offset(p::Pose{D,F}, r::Integer=0, L::Integer=1) where {D,F} =
    p * standard_rotation(F, Val(D), r, L)

"""
    twistfreedom(b::BindingSite)

Return how much say `b` has in the twist of a bond it takes part in: `stab` if it locks,
`gauge` if it is rotation-free. See [`BindingSite`](@ref) and [`bondperiod`](@ref).
"""
@inline twistfreedom(b::BindingSite) = b.locking ? b.stab : b.gauge

"""
    bondperiod(b1::BindingSite, b2::BindingSite)

Return `L`, the number of twists about the bond axis that a bond between `b1` and `b2` admits:
the least common multiple of the two sites' [`twistfreedom`](@ref)s. See
[`standard_offset`](@ref).

The lcm rather than either site's own count, because the admissible set has to be closed under
both sites' turns at once, and the rotations of order `q₁` and `q₂` generate a cyclic group of
order `lcm(q₁, q₂)`. So one site being rotation-free is enough to open a bond up, whatever the
other one says.
"""
@inline bondperiod(b1::BindingSite, b2::BindingSite) = lcm(twistfreedom(b1), twistfreedom(b2))

"""
    nregistrations(b1::BindingSite, b2::BindingSite)

Return how many of the [`bondperiod`](@ref) twists give *distinct* structures when `b2`'s
particle is attached to `b1`, so that `0:(nregistrations(b1, b2) - 1)` enumerates them.

The ones skipped are those differing by a symmetry of the particle being attached — `b2.stab`
of them, turning about the bond axis. They put the same body in the same place with only its
sites permuted, so they give nothing new. They form a subgroup of the `L` twists, since `stab`
divides `b2.gauge` and equals it or `b2.stab` either way, hence divides `L`; that is what makes
`0:(L ÷ stab - 1)` a transversal.

Note the asymmetry: only the *attached* particle's stabiliser is quotiented out. The
polyform's own symmetries would merge structures too, but they are not knowable from one site,
and canonical form catches them. By the same token this is an optimisation rather than a
correctness requirement — the duplicates would be merged anyway — but it is where the saving
is.
"""
@inline nregistrations(b1::BindingSite, b2::BindingSite) = bondperiod(b1, b2) ÷ b2.stab

"""
    contact_pairing(vs1::UnitRange, vs2::UnitRange, r=0, L=1)

Return the pairs of graph vertices to join when two bonding binding sites occupy the vertex
ranges `vs1` and `vs2` and meet in registration `r` of `L` (see [`standard_offset`](@ref)).

Vertices are joined exactly where they land at the same angle about the bond axis. Number the
vertices of each site from zero; vertex `a` of site 1 sits at azimuth `2πa/k₁` in site 1's
frame, and vertex `b` of site 2 — whose frame is turned by `2πr/L` and then flipped, so its
cyclic order runs the other way round — sits at `2πr/L - 2πb/k₂`. Writing `K = lcm(k₁, k₂)`,
`s₁ = K/k₁`, `s₂ = K/k₂`, the coincidences are the solutions of

    a·s₁ + b·s₂ ≡ t   (mod K),      t = r·K/L

`t` is always an integer: a site's gauge divides its vertex count, so `L = lcm(g₁, g₂)`
divides `K`. Since `gcd(s₁, s₂) = 1`, the congruence has exactly `gcd(k₁, k₂)` solutions for
every `t`, and distinct registrations give disjoint solution sets — so the pairing count never
depends on the registration, and the graph records which registration a bond is in.

Two consequences are worth naming. The bond keeps its residual symmetry: a site with `k`
vertices is invariant under turns by `2π/k`, so the bond is invariant under the turns common
to both, `C_k₁ ∩ C_k₂ = C_gcd(k₁,k₂)`, and the `gcd` links are what let nauty see it. And at
`r = 0` this is the counter-rotating pairing with the two first vertices as the fixed point,
so 2D species and patchy particles, which have one registration, are unaffected.

This is the single place where the bond convention is defined; a species needing a different
registry would change it here.
"""
@inline function contact_pairing(vs1::UnitRange{Int}, vs2::UnitRange{Int},
                                 r::Integer=0, L::Integer=1)
    k1, k2 = length(vs1), length(vs2)
    G = gcd(k1, k2)
    K = k1 * k2 ÷ G
    a0, b0 = _registration_shift(k1, k2, K, r, L)
    return (vs1[1 + mod(a0 + j * (k1 ÷ G), k1)] => vs2[1 + mod(b0 - j * (k2 ÷ G), k2)]
            for j in 0:(G - 1))
end

# One solution of a*s1 + b*s2 = r*K/L (mod K); the rest follow by stepping the two ranges in
# opposite directions. Zero-based, and `(0, 0)` for the base registration.
@inline function _registration_shift(k1::Int, k2::Int, K::Int, r::Integer, L::Integer)
    r == 0 && return 0, 0
    K % L == 0 || throw(ArgumentError(
        "a bond in $L registrations cannot be encoded by sites of $k1 and $k2 graph vertices: " *
        "$L must divide lcm($k1, $k2) = $K. Give the sites more vertices, i.e. a graph built " *
        "by `dartencoding` rather than `cycleencoding`, or declare a smaller gauge."
    ))
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
    return isapprox(b1.pose.x, b2.pose.x;
                    atol=b1.touching_tolerance + b2.touching_tolerance, rtol=0)
end

"""
    registration(b1::BindingSite, b2::BindingSite)

Return which of the bond's registrations the two sites are in, or `nothing` if their
orientations are not related by any of them.

This is what lets a bond be recognised however it arose. A ring closure is not placed by
[`raise!`](@ref) — it is discovered, when a site of the new particle happens to land on one
already in the polyform — so the registration it closes in has to be read back off the
geometry, and the same value then decides the vertex pairing via [`contact_pairing`](@ref).
Symmetric in its arguments: `psi₂ = psi₁·Rx(2πr/L)·Δ` gives back the same `r` read the other
way round, since `Δ` is a π rotation that inverts the twist and squares to the identity.
"""
function registration(b1::BindingSite, b2::BindingSite)
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
differ by a standard offset, in any registration the bond admits.
"""
isaligned(b1::BindingSite, b2::BindingSite) = !isnothing(registration(b1, b2))

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