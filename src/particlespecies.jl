"""
    ParticleSpecies{D,F,P<:Pose{D,F}}

Abstract supertype of all particle species. `D` is the spatial dimension, `F` is the
numeric type, and `P` is the concrete `Pose` type used by the species.
"""
abstract type ParticleSpecies{D,F,P<:Pose{D,F}} end

"""
    numtype(::ParticleSpecies{D,F,P}) -> Type

Return the numeric element type `F`.
"""
numtype(::ParticleSpecies{D,F}) where {D,F} = F
numtype(::Type{<:ParticleSpecies{D,F}}) where {D,F} = F

"""
    posetype(::ParticleSpecies{D,F,P}) -> Type

Return the concrete `Pose` type `P` used by the particle species.
"""
posetype(::ParticleSpecies{D,F,P}) where {D,F,P} = P
posetype(::Type{<:ParticleSpecies{D,F,P}}) where {D,F,P} = P

const SpeciesAndPose{SPC} = Union{Pair{SPC,<:Pose},Tuple{SPC,<:Pose}} where {SPC<:ParticleSpecies}

"""
    bindingsites(p::ParticleSpecies, i::Integer)

Return the `i`th binding site of particle species `p`.
"""
function bindingsites(::ParticleSpecies, ::Integer) end

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
function dimension(::ParticleSpecies) end

"""
    graphrep(p::ParticleSpecies)

Return the graph representation of the particle species `p`.
"""
function graphrep(::ParticleSpecies) end

"""
    nsites(p::ParticleSpecies)

Return the number of binding sites of a particle of species `p`.
"""
function nsites(::ParticleSpecies) end

"""
    setcolor!(p::ParticleSpecies, c::Integer)

Set the base color of particle species `p` to `c`.

The colors of the binding sites will then start at `c`.
"""
function setcolor!(::ParticleSpecies, ::Integer) end

"""
    isconvex(::ParticleSpecies)

Return true if the particle species has a convex shape, which enables minor optimizations when checking
for overlaps.
"""
function isconvex(::ParticleSpecies) 
    return false
end

"""
    could_contact(p1::SpeciesAndPose, p2::SpeciesAndPose)

Return `true` if the particle species at their respective poses could potentially be in contact.

This should be a quick and efficient pre-check; a value of `true` does not necessarily mean that
contact has occured. If `false`, however, it is assumed that contact is impossible and no
further contact checks will be performed.
"""
function could_contact(::SpeciesAndPose, ::SpeciesAndPose) end

"""
    overlap(p1::SpeciesAndPose, p2::SpeciesAndPose)

Return `true` if the particle species at their respective poses are overlapping.
"""
function overlap(::SpeciesAndPose, ::SpeciesAndPose) end

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

Return the list of equivalent binding sites, that is equivalence classes of binding sites under the 
symmetry of the particle species.
"""
function equivalent_site_indices(spcs::ParticleSpecies)
    symmetrynumber(spcs) == 1 && return [[i] for i in 1:nsites(spcs)]
    throw(ErrorException("building blocks with symmetry are not yet implemented"))
end