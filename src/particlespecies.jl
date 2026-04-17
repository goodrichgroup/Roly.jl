struct BindingSite{P}
    pose::P
    color::Int
    vertices::UnitRange{Int}
end
function Base.show(io::Core.IO, b::BindingSite)
    print(io, "$(dimension(b.pose))-dimensional BindingSite:\n")
    print(io, " - color: \t")
    println(io, b.color)
    print(io, " - vertices:\t$(b.vertices)\n")
    #TODO make this general
    print(io, " - x: $(b.pose.x)\n")
    print(io, " - ψ: $(rotation_angle(b.pose.ψ) / π)π")
end

Base.:*(p::Pose, site::BindingSite) = typeof(site)(p * site.pose, site.color, site.vertices)
Base.:*(site::BindingSite, p::Pose) = typeof(site)(site.pose * p, site.color, site.vertices)
shift_vertices(site::BindingSite, v::Integer) = typeof(site)(site.pose, site.color, site.vertices .+ v)
shift_color(site::BindingSite, c::Integer) = typeof(site)(site.pose, site.color + c, site.vertices)
@inline orientation_offset(b::BindingSite{<:Pose{2,F}}) where {F} = b.pose * Angle2d{F}(π)
@inline orientation_offset(b::BindingSite{<:Pose{3,F}}) where {F} =  b.pose * RotXYZ{F}(0, 0, π)

color(b::BindingSite) = b.color
function istouching(b1::BindingSite, b2::BindingSite; kwargs...)
    return isapprox(b1.pose.x, b2.pose.x; kwargs...)
end
function isaligned(b1::BindingSite, b2::BindingSite; kwargs...)
    return isapprox(b1.pose.ψ, orientation_offset(b2).ψ; kwargs...)
end
function isincontact(b1::BindingSite, b2::BindingSite; kwargs...)
    return isapprox(b1.pose, b2.pose; kwargs...)
end

abstract type ParticleSpecies end
const SpeciesAndPose{SPC} = Union{Pair{SPC,<:Pose},Tuple{SPC,<:Pose}} where {SPC<:ParticleSpecies}

"""
    could_contact(p1::SpeciesAndPose, p2::SpeciesAndPose)

Perform a quick check to determine whether two particle species at specific poses could
potentially be in contact. If `could_contact` returns `false`, it is assumed that contact
is impossible and no further contact checks will be performed.
"""
function graphrep(p::ParticleSpecies) end
function dimension(::ParticleSpecies) end
function nsites(p::ParticleSpecies) end
function bindingsite(p::ParticleSpecies, i::Integer) end
function bindingsites(p::ParticleSpecies)
    return (bindingsite(p, i) for i in 1:nsites(p))
end
function setcolor!(p::ParticleSpecies, c::Integer) end

function symmetrynumber(p::ParticleSpecies)
    _, autg = nauty(graphrep(p))
    return convert(Int, autg.n)
end

"""
    isconvex(::ParticleSpecies)

Return true if the particle species has a convex shape, which enables minor optimizations when checking
for overlaps.
"""
function isconvex(::ParticleSpecies) 
    return false
end

function could_contact(p1::SpeciesAndPose, p2::SpeciesAndPose) end
function overlap(p1::SpeciesAndPose, p2::SpeciesAndPose) end

"""
    equivalent_site_indices(spcs::ParticleSpecies)

Return the list of equivalent binding sites, that is equivalence classes of binding sites under the 
symmetry of the particle species.
"""
function equivalent_site_indices(spcs::ParticleSpecies)
    symmetrynumber(spcs) == 1 && return [[i] for i in 1:nsites(spcs)]
    throw(ErrorException("building blocks with symmetry are not yet implemented"))
end