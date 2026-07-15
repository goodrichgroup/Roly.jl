"""
    ParticleSpecies{D,B<:BindingSite}

Abstract supertype of all particle species. `D` is the spatial dimension and `B` is the
concrete `BindingSite` type.
"""
abstract type ParticleSpecies{D,B<:BindingSite} end

"""
    numtype(::ParticleSpecies)

Return the numeric element type, derived from the `BindingSite` type parameter.
"""
numtype(::ParticleSpecies{D,B}) where {D,B} = eltype(posetype(B))
numtype(::Type{<:ParticleSpecies{D,B}}) where {D,B} = eltype(posetype(B))

"""
    posetype(::ParticleSpecies)

Return the concrete `Pose` type used by the particle species.
"""
posetype(::ParticleSpecies{D,B}) where {D,B} = posetype(B)
posetype(::Type{<:ParticleSpecies{D,B}}) where {D,B} = posetype(B)

const SpeciesAndPose{SPC} = Pair{SPC,<:Pose} where {SPC<:ParticleSpecies}

"""
    bindingsites(p::ParticleSpecies, i::Integer)

Return the `i`th binding site of particle species `p`.
"""
function bindingsites end

"""
    bindingsites(p::ParticleSpecies)

Return an iterator over all binding sites of a particle of species `p`.
"""
function bindingsites(p::ParticleSpecies)
    return (bindingsites(p, i) for i in 1:nsites(p))
end

"""
    dimension(p::ParticleSpecies)

Return the spatial dimension a particle of species `p` lives in.
"""
dimension(::ParticleSpecies{D}) where {D} = D

"""
    graphrep(p::ParticleSpecies)

Return the graph representation of the particle species `p`.
"""
function graphrep end

"""
    nsites(p::ParticleSpecies)

Return the number of binding sites of a particle of species `p`.
"""
function nsites end

"""
    setcolors!(p::ParticleSpecies, colors::AbstractVector{<:Integer})

Assign colors to the binding sites of particle species `p`.
"""
function setcolors! end

"""
    isconvex(::ParticleSpecies)

Return true if the particle species has a convex shape, which enables minor optimizations when checking
for overlaps.
"""
function isconvex(::ParticleSpecies)
    return false
end

"""
    bounding_radius(p::ParticleSpecies)

Return the radius of a bounding sphere centred at the particle's pose origin.
"""
function bounding_radius end

"""
    could_contact(p1::SpeciesAndPose, p2::SpeciesAndPose)

Return `true` if the particle species at their respective poses could potentially be in contact.

This should be a quick and efficient pre-check; a value of `true` does not necessarily mean that
contact has occured. If `false`, however, it is assumed that contact is impossible and no
further contact checks will be performed.
"""
function could_contact(p1::SpeciesAndPose, p2::SpeciesAndPose; kwargs...)
    (s1, pose1), (s2, pose2) = p1, p2
    return norm(pose1.x - pose2.x) < bounding_radius(s1) + bounding_radius(s2)
end

"""
    overlap(p1::SpeciesAndPose, p2::SpeciesAndPose)

Return `true` if the particle species at their respective poses are overlapping.
"""
function overlap end

"""
    symmetrynumber(p::ParticleSpecies)

Return the symmetry number of the particle species `p`.

The symmetry number is equal to the size of the automorphism group
of `graphrep(p)`.
"""
function symmetrynumber(p::ParticleSpecies)
    _, autg = nauty(graphrep(p))
    return convert(Int, autg.n)
end

"""
    equivalent_site_indices(spcs::ParticleSpecies)

Return a list of equivalent binding sites, that is equivalence classes of binding sites under the 
symmetry of the particle species.
"""
function equivalent_site_indices(spcs::ParticleSpecies)
    symmetrynumber(spcs) == 1 && return [[i] for i in 1:nsites(spcs)]
    throw(ErrorException("building blocks with symmetry are not yet implemented"))
end