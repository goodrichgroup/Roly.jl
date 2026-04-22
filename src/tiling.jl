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
    return !isnothing(lattice_vectors(poly; kwargs...))
end

function lattice_vectors(poly::Polyform; nreps=1)
    particles = poly.particles

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
        if vec_idxs[end] == nvecs
            pop!(vec_idxs)
            isempty(vec_idxs) && return nothing

            pop!(tile_contacts)
            pop!(tile_particles)
        end
        vec_idxs[end] += 1
    end
    return latvecs[vec_idxs[1:end-1]]
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