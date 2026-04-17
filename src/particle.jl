"""
    Particle{P<:Pose,SPC<:ParticleSpecies}

A `Particle` describes an indivisible subunit that makes up a `Polyform`.

Every particle is a member of a `ParticleSpecies` and is located in space
at a specific `Pose`.
"""
struct Particle{P<:Pose,SPC<:ParticleSpecies}
    pose::P
    leading_vertex::Int
    species::SPC
end

"""
    Particle(ps::ParticleSpecies, pose=nothing; leading_vertex::Integer)

Create a particle of species `ps` at a given `pose`.
"""
function Particle(ps::ParticleSpecies{D,F}, pose=nothing; leading_vertex::Integer) where {D,F}
    if isnothing(pose)
        pose = Pose{D,F}()
    end
    return Particle(pose, leading_vertex, ps)
end

# TODO: this should probably be a function on a polyform?
graphvertices(p::Particle) = (1:nv(graphrep(species(p)))) .+ (p.leading_vertex - 1)

"""
    leading_vertex(p::Particle)

Return the leading vertex of particle `p`.
"""
leading_vertex(p::Particle) = p.leading_vertex

"""
    species(p::Particle)

Return the particle species of particle `p`.
"""
species(p::Particle) = p.species

"""
    nsites(p::Particle)

Return the number of binding sites of particle `p`.
"""
nsites(p::Particle) = nsites(species(p))

"""
    bindingsites(p::Particle, i::Integer)

Return the `i`th binding site of particle `p`.
"""
function bindingsites(p::Particle, i::Integer) 
    return shift_vertices(bindingsites(species(p), i), leading_vertex(p) - 1) * p.pose
end

"""
    bindingsites(p::Particle)

Return an iterator over the binding sites of particle `p`.
"""
function bindingsites(p::Particle)
    return (bindingsites(p, i) for i in 1:nsites(p))
end

"""
    could_contact(p1::Particle, p2::Particle)

Return `true` if the particles could potentially be in contact.
"""
function could_contact(p1::Particle, p2::Particle)
    return could_contact(p1.species => p1.pose, p2.species => p2.pose)
end

"""
    overlap(p1::Particle, p2::Particle)

Return `true` if the particles are overlapping.
"""
function overlap(p1::Particle, p2::Particle)
    return overlap(p1.species => p1.pose, p2.species => p2.pose)
end

"""
    isconvex(::Particle)

Return true if the particle has a convex shape, which enables minor optimizations when checking
for overlaps.
"""
isconvex(p::Particle) = isconvex(p.species)

function Base.show(io::Core.IO, p::Particle{P}) where {P}
    print(io, "Particle[cs=$([color(b) for b in bindingsites(p)])]\n")
end
function Base.show(io::Core.IO, ::MIME"text/plain", p::Particle)
    print(io, "$(dimension(P))-dimensional Particle:\n")
    print(io, " - colors: $([color(b) for b in bindingsites(p)])")
    print(io, " - vertices:\t$(graphvertices(p))")
end


function render!(ax, p::Particle; kwargs...)
    return render!(ax, species(p), p.pose; kwargs...)
end
function render(p::Particle; kwargs...)
    return render(species(p), p.pose; kwargs...)
end

