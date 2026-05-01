function cantile(poly::AbstractVector{<:Polyform}; maxtilesize, return_tile=false, kwargs...)
    metasys = AssemblySystem(poly)
    tile_found = Ref(false)
    tile = Ref(Polyform(metasys))

    f(s, args...) = isunitcell(s; kwargs...) ? (tile_found[] = true; tile[] = copy(s); BREAK) : ACCEPT
    _, largest_size = polyenum(f, metasys; maxsize=maxtilesize)
    
    tile_found[] && return return_tile ? tile[] : true
    largest_size < maxtilesize && return false
    return missing
end
cantile(poly::Polyform; kwargs...) = cantile([poly]; kwargs...)
cantile(poly, bondslots; kwargs...) = cantile(poly; bondslots, kwargs...)

#TODO we also need to return true if it can self-close (if it is a loop)

function isunitcell(poly::Polyform; kwargs...)
    return !isnothing(lattice_vectors(poly; kwargs...))
end

function lattice_vectors(poly::Polyform; bondslots=nothing, nreps=1)
    particles = poly.particles
    sys = assemblysystem(poly)

    slots = if !isnothing(bondslots)
        ns = nspecies(sys)
        nb_orig = length(_original_bonded_colors(sys))
        slots = zeros(Int, nb_orig)
        for (i, n) in pairs(composition(poly)[1:ns])
            slots[bondslots[i]] .+= n
        end
        slots
    else
        nothing
    end

    #TODO this should be in potential_lattice_vectors
    n_opensites = isnothing(bondslots) ? length(Roly.collect_open_bindingsites(poly)) : min(2*sum(slots), length(Roly.collect_open_bindingsites(poly)))
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
                    overlap, _ = Roly._overlap_and_contacts(tile_ps, shifted_particle, sys)
                    overlap && @goto upward_traverse
                end

                overlap, contacts = Roly.overlap_and_contacts(poly, shifted_particle)
                overlap && @goto upward_traverse
                append!(latest_contacts, contacts)
            end

            append!(latest_particles, shifted_particles)
        end

        # TODO check if contacts are interacting!!
        if !isnothing(slots)
            labs = labels(graphrep(poly))
            orig_bondcolors = _original_bonded_colors(sys)
            bonds = zero(slots)
            combined_contacts = copy(tile_contacts)
            push!(combined_contacts, latest_contacts)
            for (vs1, vs2) in Iterators.flatten(combined_contacts)
                v1, v2 = first(vs1), first(vs2)
                c1 = _original_color(sys, label2color(sys, labs[v1]))
                c2 = _original_color(sys, label2color(sys, labs[v2]))
                bidx = findfirst(==(minmax(c1, c2)), orig_bondcolors)
                isnothing(bidx) && @goto upward_traverse
                bonds[bidx] += 1
            end
            any(<(0), slots - bonds) && @goto upward_traverse
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