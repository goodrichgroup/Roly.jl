function cantile(poly::AbstractVector{<:Polyform}; maxtilesize, kwargs...)
    tile_found = Ref(false)
    f(s, args...) = isunitcell(s; kwargs...) ? (tile_found[] = true; BREAK) : ACCEPT
    _, largest_size = polyenum(f, AssemblySystem(poly); maxsize=maxtilesize)
    
    tile_found[] && return true
    largest_size < maxtilesize && return false
    return missing
end
cantile(poly::Polyform; kwargs...) = cantile([poly]; kwargs...)

function isunitcell(poly::Polyform; kwargs...)
    return !isnothing(tile_latvecs(poly; kwargs...))
end

function tile_latvecs(poly::Polyform; kwargs...)
    out = latvecs_and_contacts(poly; kwargs...)
    isnothing(out) && return nothing
    return out[1]
end

function tile_contacts(poly::Polyform; kwargs...)
    out = latvecs_and_contacts(poly; kwargs...)
    isnothing(out) && return nothing
    return out[2]
end

function tile_bonds(poly::Polyform; kwargs...)
    # tcs = tile_contacts(poly; kwargs...)
    # isnothing(tcs) && return nothing

    tbonds = Vector{Int}[]
    for (_, tcs) in LatticeIter(poly)
        bonds = Int[]
        for contacts in tcs
            for contact in contacts
                idx = bondindex(poly, contact)
                push!(bonds, idx)
            end
        end
        push!(tbonds, bonds)
    end
    return tbonds
end

function bondindex(poly, contact)
    sys = assemblysystem(poly)
    labs = labels(graphrep(poly))
    
    vs1, vs2 = contact
    v1, v2 = first(vs1), first(vs2)

    v1 = tocanon(poly, v1)
    v2 = tocanon(poly, v2)

    c1, c2 = label2color(sys, labs[v1]), label2color(sys, labs[v2])
    return findfirst(==(minmax(c1, c2)), bonded_colors(sys))
end

struct LatticeIter{D,P,F}
    polyform::P
    nreps::Int
    potential_latvecs::Vector{SVector{D,F}}
    n_opensites::Int
    function LatticeIter(poly::Polyform; nreps=1)
        # TODO FLOATTYPE
        return new{dimension(assemblysystem(poly)),typeof(poly),Float64}(poly, nreps, potential_lattice_vectors(poly), length(Roly.collect_open_bindingsites(poly)))
    end
end

Base.IteratorSize(::Type{<:LatticeIter}) = Base.SizeUnknown()
Base.eltype(::Type{<:LatticeIter{D,Any,F}}) where {D,F} = Tuple{Vector{SVector{D,F}},Vector{Vector{NTuple{2,UnitRange{Int}}}}}

struct LatticeIterState{PT}
    vec_idxs::Vector{Int}
    particles::Vector{Vector{PT}}
    contacts::Vector{Vector{NTuple{2,UnitRange{Int}}}}
end

function Base.iterate(latiter::LatticeIter)
    state = LatticeIterState{eltype(latiter.polyform.particles)}([1], [], [])
    return Base.iterate(latiter, state)
end

function Base.iterate(latiter::LatticeIter, state::LatticeIterState)
    poly = latiter.polyform
    sys = assemblysystem(poly)
    nreps = latiter.nreps

    latvecs = latiter.potential_latvecs
    nvecs = length(latvecs)

    n_opensites = latiter.n_opensites

    n_opensites > 1 || return nothing
    nvecs > 0 || return nothing

    particles = poly.particles
    
    vec_idxs, tile_particles, tile_contacts = state.vec_idxs, state.particles, state.contacts

    while true
        latest_contacts = []
        latest_particles = []

        k = length(vec_idxs)
        for m in Iterators.product((0:nreps for _ in 1:k-1)..., 1:nreps)
            all(iszero, m) && continue

            shifted_particles = particles .+ Ref(sum(m .* latvecs[vec_idxs])) # TODO

            for shifted_particle in shifted_particles
                for tile_ps in tile_particles
                    overlap, contacts = Roly._overlap_and_contacts(tile_ps, shifted_particle, sys)
                    overlap && @goto upward_traverse
                    for contact in contacts
                        isnothing(bondindex(poly, contact)) && @goto upward_traverse
                    end
                end

                overlap, contacts = Roly.overlap_and_contacts(poly, shifted_particle)
                overlap && @goto upward_traverse
                for contact in contacts
                    isnothing(bondindex(poly, contact)) && @goto upward_traverse
                end
                append!(latest_contacts, contacts)
            end

            append!(latest_particles, shifted_particles)
        end
   
        push!(tile_contacts, latest_contacts)
        push!(tile_particles, latest_particles)
        push!(vec_idxs, vec_idxs[end])
        return (latvecs[vec_idxs[1:end-1]], copy.(tile_contacts)), state

        @label upward_traverse
        while vec_idxs[end] == nvecs || sum(length, tile_contacts; init=0) == n_opensites ÷ 2
            pop!(vec_idxs)
            isempty(vec_idxs) && return nothing

            pop!(tile_contacts)
            pop!(tile_particles)
        end
        vec_idxs[end] += 1
    end
    error()
end

function latvecs_and_contacts(poly::Polyform; nreps=1)
    particles = poly.particles
    render(poly)

    #TODO this should be in potential_lattice_vectors
    n_opensites = length(Roly.collect_open_bindingsites(poly))
    n_opensites > 1 || return nothing

    latvecs = potential_lattice_vectors(poly)
    isempty(latvecs) && return nothing
    nvecs = length(latvecs)
    
    tile_particles = Vector{eltype(particles)}[]
    tile_contacts = Vector{NTuple{2,UnitRange{Int}}}[]
    vec_idxs = [1]

    while sum(length, tile_contacts; init=0) < n_opensites ÷ 2
        latest_contacts = []
        latest_particles = []

        k = length(vec_idxs)
        for m in Iterators.product((0:nreps for _ in 1:k-1)..., 1:nreps)
            all(iszero, m) && continue

            shifted_particles = particles .+ Ref(sum(m .* latvecs[vec_idxs])) # TODO

            for shifted_particle in shifted_particles
                for tile_ps in tile_particles
                    overlap, _ = Roly._overlap_and_contacts(tile_ps, shifted_particle, assemblysystem(poly))
                    overlap && @goto upward_traverse
                end

                overlap, contacts = Roly.overlap_and_contacts(poly, shifted_particle)
                overlap && @goto upward_traverse
                append!(latest_contacts, contacts)
            end

            append!(latest_particles, shifted_particles)
        end
   
        push!(tile_contacts, latest_contacts)
        push!(tile_particles, latest_particles)
        push!(vec_idxs, vec_idxs[end])
        continue

        @label upward_traverse
        while vec_idxs[end] == nvecs
            pop!(vec_idxs)
            isempty(vec_idxs) && return nothing

            pop!(tile_contacts)
            pop!(tile_particles)
        end
        vec_idxs[end] += 1
    end
    return latvecs[vec_idxs[1:end-1]], tile_contacts
end


function potential_lattice_vectors(poly::Polyform{D}) where {D}
    sys = assemblysystem(poly)

    # TODO
    pairs = sort!(Roly.collect_compatible_pairs(poly))
    open_sites = sort!(Roly.collect_open_bindingsites(poly))

    latvecs = SVector{D,Float64}[] #TODO

    visited_sites = Set() # TODO

    for (site, siteloc) in pairs
        c = Roly.siteloc2color(sys, siteloc)

        compatible_sites = filter(s->color(s) == c, open_sites)
        aligned_sites = filter(s->Roly.isaligned(s, site), compatible_sites)

        for aligned_site in aligned_sites
            aligned_site ∈ visited_sites && continue
            v = site.pose.x - aligned_site.pose.x
            push!(latvecs, v)
        end
        push!(visited_sites, site)
    end
    return latvecs
end