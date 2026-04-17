struct Particle{P,SPC<:ParticleSpecies}
    pose::P
    leading_vertex::Int
    species::SPC
end
function Particle(ps::ParticleSpecies, pose=nothing; leading_vertex::Integer)
    if isnothing(pose)
        pose = Pose{dimension(ps)}()
    end
    return Particle{typeof(pose),typeof(ps)}(pose, leading_vertex, ps)
end
function Particle(sys::AssemblySystem, i::Integer, pose=nothing; leading_vertex::Integer)
    ps = species(sys, i)
    return Particle(ps, pose; leading_vertex)
end

function Base.show(io::Core.IO, p::Particle{P}) where {P}
    print(io, "$(dimension(P))-dimensional Particle:\n")
    print(io, " - species: \t")
    println(io, p.species)
    print(io, " - vertices:\t$(graphvertices(p))")
end

graphvertices(p::Particle) = (1:nv(graphrep(species(p)))) .+ (p.leading_vertex - 1)
leading_vertex(p::Particle) = p.leading_vertex
species(p::Particle) = p.species
nsites(p::Particle) = nsites(species(p))

function bindingsite(p::Particle, i::Integer) 
    return shift_vertices(bindingsite(species(p), i), leading_vertex(p) - 1) * p.pose
end
function bindingsites(p::Particle)
    return (bindingsite(p, i) for i in 1:nsites(p))
end

function could_contact(p1::Particle, p2::Particle)
    return could_contact(p1.species => p1.pose, p2.species => p2.pose)
end
function overlap(p1::Particle, p2::Particle)
    return overlap(p1.species => p1.pose, p2.species => p2.pose)
end
isconvex(p::Particle) = isconvex(p.species)

function render!(ax, p::Particle; kwargs...)
    return render!(ax, species(p), p.pose; kwargs...)
end
function render(p::Particle; kwargs...)
    return render(species(p), p.pose; kwargs...)
end