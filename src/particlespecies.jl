"""
    ParticleSpecies{D,B<:BindingSite}

Abstract supertype of all particle species in `D` spatial dimensions.
"""
abstract type ParticleSpecies{D,B<:BindingSite} end

"""
    SpeciesAndPose{SPC}

A species paired with a placement, `species => pose`.

Use this type when writing methods for [`overlap`](@ref) and [`could_contact`](@ref).
"""
const SpeciesAndPose{SPC} = Pair{SPC,<:Pose} where {SPC<:ParticleSpecies}

"""
    _tilingcell(ps::ParticleSpecies)

Identify the tiling `ps` sits on, or `nothing` if it does not tile.

This function is currently considered internal and it is not recommended to define custom methods.
"""
_tilingcell(::ParticleSpecies) = nothing

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

"""
    bindingsite(p::ParticleSpecies, i::Integer)

Return the `i`th binding site of particle species `p`.
"""
function bindingsite end

"""
    bindingsites(p::ParticleSpecies)

Return an iterator over all binding sites of a particle of species `p`.
"""
function bindingsites(p::ParticleSpecies)
    return (bindingsite(p, i) for i in 1:nsites(p))
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
function setcolors!(p::ParticleSpecies, colors::AbstractVector{<:Integer})
    _recolor!(p, p.sites, colors)
    return nothing
end

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

Return the radius of a bounding sphere centerd at the particle's pose origin.
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
    # DO NOT `canonize` here: a species' graph must stay in construction order, since
    # `BindingSite.vertices` indexes it directly and `bindingsite(::Particle, ...)` reaches those
    # vertices by shifting the range by a leading vertex. Reordering here would silently
    # decouple the two.
    _, autg = nauty(graphrep(p))
    return convert(Int, autg.n)
end

# The type of `ps` by its full name, e.g. `"Roly.PolygonParticleSpecies"`, which `fromdict` resolves
# back to the type through the loaded modules.
_typetag(ps::ParticleSpecies) = join((fullname(parentmodule(typeof(ps)))..., nameof(typeof(ps))), ".")

function _typefromtag(tag::AbstractString)
    names = Symbol.(split(tag, "."))
    loaded = Base.loaded_modules_array()
    i = findfirst(m -> nameof(m) === first(names), loaded)
    T = first(names) === :Main ? Main : isnothing(i) ? nothing : loaded[i]
    isnothing(T) && throw(ArgumentError("The species type $tag cannot be read: its package is not loaded."))
    for name in names[2:end]
        (T isa Module && isdefined(T, name)) || throw(ArgumentError("$tag is not a particle species type."))
        T = getfield(T, name)
    end
    (T isa Type && T <: ParticleSpecies) || throw(ArgumentError("$tag is not a particle species type."))
    return T
end

"""
    todict(ps::ParticleSpecies)

Serialize `ps` into a `Dict{String,Any}` of strings, numbers and arrays.

This default records the type of `ps`, its number of sites and their colors, which is not enough
to rebuild it. A species type that can be rebuilt from its description defines its own method.
"""
function todict(ps::ParticleSpecies)
    return Dict{String,Any}(
        "type" => _typetag(ps), "nsites" => nsites(ps), "colors" => [color(bindingsite(ps, i)) for i in 1:nsites(ps)]
    )
end

"""
    fromdict(::Type{ParticleSpecies}, d)

Return the particle species that `d`, written by [`todict`](@ref), describes.

`d` names its species type in its entry `type`, and reading it calls `fromdict(T, d)` with that type
`T`. A species type opts in to being read by defining that method.
"""
fromdict(::Type{ParticleSpecies}, d) = fromdict(_typefromtag(d["type"]), d)

"""
    fromdict(::Type{T}, d) where {T<:ParticleSpecies}

Throw, since the species type `T` defines no method to be read back with.
"""
function fromdict(::Type{T}, d) where {T<:ParticleSpecies}
    throw(ArgumentError("$T cannot be read back, since it defines no `Roly.fromdict(::Type{<:$(nameof(T))}, d)`."))
end
