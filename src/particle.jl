"""
    Particle{P<:Pose}

A `Particle` describes an indivisible subunit that makes up a `Polyform`.

Every particle is a member of a `ParticleSpecies` (identified by `species_index` within
an `BindingRules`) and is located in space at a specific `Pose`.
"""
struct Particle{P<:Pose}
    pose::P
    leading_vertex::Int
    species_index::Int
end

"""
    Particle(sys::BindingRules, species_index::Integer, pose=nothing; leading_vertex::Integer)

Create a particle of species `species_index` within `sys` at a given `pose`.
"""
function Particle(sys::BindingRules, species_index::Integer, pose=nothing; leading_vertex::Integer)
    if isnothing(pose)
        pose = one(posetype(sys))
    end
    return Particle(pose, leading_vertex, species_index)
end

Base.:*(p, part::Particle) = typeof(part)(p * part.pose, part.leading_vertex, part.species_index)
Base.:*(part::Particle, p) = typeof(part)(part.pose * p, part.leading_vertex, part.species_index)
Base.:+(p, part::Particle) = typeof(part)(p + part.pose, part.leading_vertex, part.species_index)
Base.:+(part::Particle, p) = typeof(part)(part.pose + p, part.leading_vertex, part.species_index)

@inline graphvertices(p::Particle, sys::BindingRules) =
    (1:nv(graphrep(species(sys, p.species_index)))) .+ (p.leading_vertex - 1)

"""
    leading_vertex(p::Particle)

Return the leading vertex of particle `p`.
"""
leading_vertex(p::Particle) = p.leading_vertex

"""
    shift_leadingvertex(p::Particle, v::Integer)

Return a copy of `p` whose block of graph vertices starts `v` further along.

Used when removing a particle compacts the original vertex numbering; compare
[`shift_vertices`](@ref) for binding sites.
"""
@inline shift_leadingvertex(p::Particle, v::Integer) =
    typeof(p)(p.pose, p.leading_vertex + v, p.species_index)

"""
    species_index(p::Particle)

Return the species index of particle `p` within its `BindingRules`.
"""
species_index(p::Particle) = p.species_index

"""
    nsites(p::Particle, sys::BindingRules)

Return the number of binding sites of particle `p`.
"""
nsites(p::Particle, sys::BindingRules) = nsites(species(sys, p.species_index))

"""
    bindingsites(p::Particle, sys::BindingRules, i::Integer)

Return the `i`th binding site of particle `p`.
"""
function bindingsites(p::Particle, sys::BindingRules, i::Integer)
    return p.pose * shift_vertices(bindingsites(species(sys, p.species_index), i), leading_vertex(p) - 1)
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
    return could_contact(species(sys, p1.species_index) => p1.pose,
                         species(sys, p2.species_index) => p2.pose)
end

"""
    overlap(p1::Particle, p2::Particle, sys::BindingRules)

Return `true` if the particles are overlapping.

When every species of `sys` tiles space and every bond carries a tile onto a tile — see
[`_tiles`](@ref) — the particles of any assembly sit on cells of one tiling, distinct cells have
disjoint interiors, and overlapping therefore means occupying the same cell. Coincident centres
answer it, and no geometry runs. Every other system takes the species' own `overlap`.
"""
function overlap(p1::Particle, p2::Particle, sys::BindingRules)
    spcs1, spcs2 = species(sys, p1.species_index), species(sys, p2.species_index)
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
isconvex(p::Particle, sys::BindingRules) = isconvex(species(sys, p.species_index))

function Base.show(io::Core.IO, p::Particle)
    print(io, "Particle[si=$(p.species_index), lv=$(p.leading_vertex)]")
end
function Base.show(io::Core.IO, ::MIME"text/plain", p::Particle)
    print(io, "Particle:\n")
    print(io, " - species_index: $(p.species_index)\n")
    print(io, " - leading_vertex: $(p.leading_vertex)\n")
    print(io, " - pose: $(p.pose)")
end