"""
    Particle{P<:Pose}

A `Particle` describes an indivisible subunit that makes up a `Polyform`.

Every particle is a member of a `ParticleSpecies` (identified by `speciesindex` within
an `BindingRules`) and is located in space at a specific `Pose`.
"""
struct Particle{P<:Pose}
    pose::P
    leadingvertex::Int
    speciesindex::Int
end

"""
    Particle(sys::BindingRules, speciesindex::Integer, pose=nothing; leadingvertex::Integer)

Create a particle of species `speciesindex` within `sys` at a given `pose`.
"""
function Particle(sys::BindingRules, speciesindex::Integer, pose=nothing; leadingvertex::Integer)
    if isnothing(pose)
        pose = one(posetype(sys))
    end
    return Particle(pose, leadingvertex, speciesindex)
end

Base.:*(p, part::Particle) = typeof(part)(p * part.pose, part.leadingvertex, part.speciesindex)
Base.:*(part::Particle, p) = typeof(part)(part.pose * p, part.leadingvertex, part.speciesindex)
Base.:+(p, part::Particle) = typeof(part)(p + part.pose, part.leadingvertex, part.speciesindex)
Base.:+(part::Particle, p) = typeof(part)(part.pose + p, part.leadingvertex, part.speciesindex)

@inline graphvertices(p::Particle, sys::BindingRules) =
    (1:nv(graphrep(species(sys, p.speciesindex)))) .+ (p.leadingvertex - 1)

"""
    leadingvertex(p::Particle)

Return the leading vertex of particle `p`.
"""
leadingvertex(p::Particle) = p.leadingvertex

"""
    shift_leadingvertex(p::Particle, v::Integer)

Return a copy of `p` whose block of graph vertices starts `v` further along.

Used when removing a particle compacts the original vertex numbering; compare
[`shift_vertices`](@ref) for binding sites.
"""
@inline shift_leadingvertex(p::Particle, v::Integer) = typeof(p)(p.pose, p.leadingvertex + v, p.speciesindex)

"""
    speciesindex(p::Particle)

Return the species index of particle `p` within its `BindingRules`.
"""
speciesindex(p::Particle) = p.speciesindex

"""
    nsites(p::Particle, sys::BindingRules)

Return the number of binding sites of particle `p`.
"""
nsites(p::Particle, sys::BindingRules) = nsites(species(sys, p.speciesindex))

"""
    bindingsites(p::Particle, sys::BindingRules, i::Integer)

Return the `i`th binding site of particle `p`.
"""
function bindingsites(p::Particle, sys::BindingRules, i::Integer)
    return p.pose * shift_vertices(bindingsites(species(sys, p.speciesindex), i), leadingvertex(p) - 1)
end

"""
    bindingsites(p::Particle, sys::BindingRules)

Return an iterator over the binding sites of particle `p`.
"""
function bindingsites(p::Particle, sys::BindingRules)
    return (bindingsites(p, sys, i) for i in 1:nsites(p, sys))
end

"""
    could_contact(p1::Particle, p2::Particle, sys::BindingRules)

Return `true` if the particles could potentially be in contact.
"""
function could_contact(p1::Particle, p2::Particle, sys::BindingRules)
    return could_contact(species(sys, p1.speciesindex) => p1.pose, species(sys, p2.speciesindex) => p2.pose)
end

"""
    overlap(p1::Particle, p2::Particle, sys::BindingRules)

Return `true` if the particles are overlapping.
"""
function overlap(p1::Particle, p2::Particle, sys::BindingRules)
    spcs1, spcs2 = species(sys, p1.speciesindex), species(sys, p2.speciesindex)

    # Binding rules can override the overlap check
    if sys._onlattice
        atol = sqrt(eps(numtype(spcs1))) * (bounding_radius(spcs1) + bounding_radius(spcs2))
        return norm(p1.pose.x - p2.pose.x) < atol
    end
    return overlap(spcs1 => p1.pose, spcs2 => p2.pose)
end

"""
    isconvex(p::Particle, sys::BindingRules)

Return true if the particle has a convex shape, which enables minor optimizations when
checking for overlaps.
"""
isconvex(p::Particle, sys::BindingRules) = isconvex(species(sys, p.speciesindex))

function Base.show(io::Core.IO, p::Particle)
    print(io, "Particle[si=$(p.speciesindex), lv=$(p.leadingvertex)]")
end
function Base.show(io::Core.IO, ::MIME"text/plain", p::Particle)
    print(io, "Particle:\n")
    print(io, " - speciesindex: $(p.speciesindex)\n")
    print(io, " - leadingvertex: $(p.leadingvertex)\n")
    print(io, " - pose: $(p.pose)")
end