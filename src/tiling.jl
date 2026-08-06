# Convert a tiling contact (original-space vertex ranges) to a bond index.
_contact_bondindex(poly, (vs1, vs2)) = bondindex(poly, tocanon(poly, first(vs1)), tocanon(poly, first(vs2)))

# ── LatticeIter ───────────────────────────────────────────────────────────────

"""
    LatticeIter(poly::Polyform; nreps=1)

An iterator over all valid unit-cell tilings of `poly` by pure translation.
Each element is a tuple

    (latvecs, contacts)

where `latvecs::Vector{SVector{D,F}}` are the lattice vectors and
`contacts::Vector{Vector{NTuple{2,UnitRange{Int}}}}` are the inter-cell contacts
for each lattice direction.

A tiling is valid when the shifted copies of `poly` tile the plane without overlap
and with `n_opensites ÷ 2` valid inter-cell bonds (one per open binding site pair).

# Example
```julia
iter = LatticeIter(poly)
for (latvecs, contacts) in iter
    ...
end
```
"""
struct LatticeIter{D,P<:Polyform,F}
    polyform::P
    nreps::Int
    potential_latvecs::Vector{SVector{D,F}}
    n_opensites::Int
    function LatticeIter(poly::Polyform; nreps=1)
        F = numtype(poly)
        D = dimension(bindingrules(poly))
        latvecs = _potential_lattice_vectors(poly)
        n_open = length(collect_open_bindingsites(poly))
        return new{D,typeof(poly),F}(poly, nreps, latvecs, n_open)
    end
end

Base.IteratorSize(::Type{<:LatticeIter}) = Base.SizeUnknown()

function Base.eltype(::Type{LatticeIter{D,P,F}}) where {D,P,F}
    return Tuple{Vector{SVector{D,F}},Vector{Vector{NTuple{2,UnitRange{Int}}}}}
end

function Base.show(io::IO, iter::LatticeIter{D}) where {D}
    nvecs = length(iter.potential_latvecs)
    return print(
        io, "LatticeIter{$D}[$nvecs candidate vector$(nvecs == 1 ? "" : "s"), ", "$(iter.n_opensites) open sites]"
    )
end

struct LatticeIterState{PT}
    vec_idxs::Vector{Int}
    particles::Vector{Vector{PT}}
    contacts::Vector{Vector{NTuple{2,UnitRange{Int}}}}
end

function Base.iterate(iter::LatticeIter)
    PT = eltype(iter.polyform.particles)
    state = LatticeIterState{PT}([1], Vector{PT}[], Vector{NTuple{2,UnitRange{Int}}}[])
    return Base.iterate(iter, state)
end

function Base.iterate(iter::LatticeIter, state::LatticeIterState)
    iter.n_opensites > 1 || return nothing
    !isempty(iter.potential_latvecs) || return nothing

    poly = iter.polyform
    sys = bindingrules(poly)
    latvecs = iter.potential_latvecs
    nvecs = length(latvecs)
    nreps = iter.nreps
    n_opensites = iter.n_opensites
    particles = poly.particles
    PT = eltype(particles)

    vec_idxs, tile_particles, tile_contacts = state.vec_idxs, state.particles, state.contacts

    while true
        latest_contacts = NTuple{2,UnitRange{Int}}[]
        latest_particles = PT[]

        k = length(vec_idxs)
        for m in Iterators.product((0:nreps for _ in 1:(k - 1))..., 1:nreps)
            all(iszero, m) && continue
            shift = sum(mi * vi for (mi, vi) in zip(m, latvecs[vec_idxs]))
            shifted_particles = particles .+ Ref(shift)

            for shifted_particle in shifted_particles
                for tile_ps in tile_particles
                    has_overlap, contacts = _overlap_and_contacts(tile_ps, shifted_particle, sys)
                    has_overlap && @goto upward_traverse
                    # reject if any bonds are not allowed by the binding rules
                    any(isnothing(_contact_bondindex(poly, c)) for c in contacts) && @goto upward_traverse
                end
                has_overlap, contacts = overlap_and_contacts(poly, shifted_particle)
                has_overlap && @goto upward_traverse
                # reject if any bonds are not allowed by the binding rules
                any(isnothing(_contact_bondindex(poly, c)) for c in contacts) && @goto upward_traverse
                append!(latest_contacts, contacts)
            end
            append!(latest_particles, shifted_particles)
        end

        push!(tile_contacts, latest_contacts)
        push!(tile_particles, latest_particles)
        push!(vec_idxs, vec_idxs[end])
        sum(length, tile_contacts; init=0) == n_opensites ÷ 2 || continue
        return (latvecs[vec_idxs[1:(end - 1)]], copy.(tile_contacts)), state

        @label upward_traverse
        while vec_idxs[end] == nvecs || sum(length, tile_contacts; init=0) == n_opensites ÷ 2
            pop!(vec_idxs)
            isempty(vec_idxs) && return nothing
            pop!(tile_contacts)
            pop!(tile_particles)
        end
        vec_idxs[end] += 1
    end
end

# Return all candidate lattice vectors: translation vectors between pairs of open
# binding sites that could tile under a uniform shift.
function _potential_lattice_vectors(poly::Polyform{D}) where {D}
    F = numtype(poly)
    sys = bindingrules(poly)
    pairs = sort!(collect_compatible_pairs(poly))
    open_sites = sort!(collect_open_bindingsites(poly))
    latvecs = SVector{D,F}[]
    visited = Set{typeof(first(open_sites))}()
    for (site, siteloc) in pairs
        c = siteloc2color(sys, siteloc)
        for aligned_site in filter(s -> color(s) == c && isaligned(s, site), open_sites)
            aligned_site ∈ visited && continue
            push!(latvecs, site.pose.x - aligned_site.pose.x)
        end
        push!(visited, site)
    end
    return latvecs
end

"""
    istranslationtile(poly::Polyform; nreps=1)

Return `true` if `poly` tiles the plane periodically by pure translation.
"""
istranslationtile(poly::Polyform; kwargs...) = !isempty(LatticeIter(poly; kwargs...))

"""
    cantile(polys; maxtilesize, kwargs...)

Check whether any meta-structure assembled from polyforms `polys` (up to size `maxtilesize`)
tiles the plane by translation. Returns `true` if a tiling is found, `false` if the
search is exhausted without finding one, or `missing` if the search reached
`maxtilesize` without a definitive answer.
"""
function cantile(polys::AbstractVector{<:Polyform}; maxtilesize, kwargs...)
    tile_found = Ref(false)
    f(s, args...) = istranslationtile(s; kwargs...) ? (tile_found[] = true; BREAK) : ACCEPT
    _, largest_size = polyenum(f, BindingRules(polys); maxsize=maxtilesize)
    tile_found[] && return true
    largest_size < maxtilesize && return false
    return missing
end

"""
    cantile(poly::Polyform; kwargs...)

Check whether `poly` can tile the plane. Returns `true` if a tiling is found, `false` if the
search is exhausted without finding one, or `missing` if the search reached
`maxtilesize` without a definitive answer.
"""
cantile(poly::Polyform; kwargs...) = cantile([poly]; kwargs...)

