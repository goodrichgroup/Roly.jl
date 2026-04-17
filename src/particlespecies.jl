abstract type ParticleSpecies{D,F} end
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