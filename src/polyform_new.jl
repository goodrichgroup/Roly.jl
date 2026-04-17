"""
    Polyform{D,T,F} where {D,T,F}

`Polyform` is a struct containing all information required to instantiate a
polyform, which here is defined as an aggregate of building blocks that are connected
(bound) together at various binding sites. In a self-assembly context, a polyform is
often referred to a "self-assembled structure".
"""
mutable struct Polyform{D,P<:Particle,S<:AssemblySystem,G<:AbstractNautyGraph}
    graphrep::G
    sigma::Int
    canonvs::Vector{Int}
    particles::Dict{Int,P}
    assemblysystem::S
end
function Polyform(sys::AssemblySystem{D}) where {D}
    P = posetype(sys)
    g = NautyDiGraph(0)
    ps = species(sys, 1)

    sigma = 1
    canonvs = []
    return Polyform{D,Particle{P,typeof(ps)},typeof(sys),typeof(g)}(g, sigma, canonvs, Dict(), sys)
end
function Polyform(sys::AssemblySystem{D}, i::Integer) where {D}
    P = posetype(sys)
    ps = species(sys, i)
    g = copy(graphrep(ps))
    canonize!(g)
    return Polyform{D,Particle{P,typeof(ps)},typeof(sys),typeof(g)}(g, symmetrynumber(ps), copy(vertices(g)),
                                                                    Dict(1=>Particle(ps; leading_vertex=1)), sys)
end

Base.copy(p::Polyform) = typeof(p)(copy(p.graphrep), p.sigma, copy(p.canonvs), copy(p.particles), p.assemblysystem)
function Base.copy!(dst::Polyform, src::Polyform)
    copy!(dst.graphrep, src.graphrep)
    dst.sigma = src.sigma
    copy!(dst.canonvs, src.canonvs)
    copy!(dst.particles, src.particles)
    dst.assemblysystem = src.assemblysystem
    return dst
end

function Base.show(io::Core.IO, p::Polyform{D}) where {D}
    print(io, "Polyform{$D}[n=$(nparticles(p))]")
end
Base.show(io::Core.IO, ::Type{Polyform{D}}) where {D} = print(io, "Polyform{$D}")
Base.:(==)(p1::Polyform, p2::Polyform) = graphrep(p1) == graphrep(p2)

@inline particle(p::Polyform, v::Integer) = get(p.particles, v, nothing)
@inline nparticles(p::Polyform) = length(p.particles)
@inline nsites(p::Polyform) = sum(nsites, values(p.particles))
@inline graphrep(p::Polyform) = p.graphrep
@inline symmetrynumber(p::Polyform) = p.sigma
@inline assemblysystem(p::Polyform) = p.assemblysystem
@inline canonical_vertices(p::Polyform) = p.canonvs
is_leadingvertex(p::Polyform, v::Integer) = v in keys(p.particles)

@inline interior_edges(p::Polyform) = _filter_edges_by_species(p, Val(true))
@inline exterior_edges(p::Polyform) = _filter_edges_by_species(p, Val(false))
function _filter_edges_by_species(p::Polyform, ::Val{same_species}) where same_species
    sys = assemblysystem(p)
    labs = labels(graphrep(p))

    filtered_edges = Iterators.filter(edges(graphrep(p))) do (; src, dst)
        src_spcs = label2species(sys, labs[src])
        dst_spcs = label2species(sys, labs[dst])
        return same_species ? src_spcs == dst_spcs : src_spcs != dst_spcs
    end
    return filtered_edges
end

"""
    bindingsite(p::Polyform, i::Integer)

Return the `i`th binding site of polyform `p`. The order of binding sites is
given by the isomorphism class of the graph encoding of `p` and is deterministic.
"""
function bindingsite(p::Polyform, i::Integer)
    k = 0
    for v in canonical_vertices(p)
        prtcl = particle(p, v)
        isnothing(prtcl) && continue

        for j in 1:nsites(prtcl)
            # TODO: skip bound or non-interacting sites
            k += 1
            if k == i
                return bindingsite(prtcl, j)
            end
        end
    end
    return nothing
end
function bindingsites(p::Polyform)
    return (bindingsite(p, i) for i in 1:nsites(p))
end

function open_bindingsites(p::Polyform)
    return Iterators.takewhile(!isnothing, (bindingsite(p, k) for k in 1:nsites(p) 
           if !isbound_site(p, k) && !isinert(assemblysystem(p), color(bindingsite(p, k)))))
end
function open_bindingsite(p::Polyform, i::Integer)
    #TODO simplify this by assinging color=0 to all inert sites
    return take_nth(open_bindingsites(p), i)
end
nopensites(p::Polyform) = sum(x->1, open_bindingsites(p); init=0)

function isbound_vertex(p::Polyform, v::Integer)
    vspcs = label2species(assemblysystem(p), label(graphrep(p), v))
    # calling `outneighbors` on a NautyGraph would allocate, this is a workaround for that
    neighbor_bitvec = NautyGraphs.adjrow(graphrep(p), v)
    return any(label2species(assemblysystem(p), label(graphrep(p), n)) != vspcs 
               for n in eachindex(neighbor_bitvec) if neighbor_bitvec[n])
end

function isbound_site(p::Polyform, site_index::Integer)
    return isbound_vertex(p, invperm(canonical_vertices(p))[first(bindingsite(p, site_index).vertices)])
end

function overlap_and_contacts(poly::Polyform, part::Particle; allow_noninteracting=false, allow_misaligned=false, kwargs...)
    intmat = interactionmatrix(assemblysystem(poly))

    contacts = NTuple{2,UnitRange{Int}}[]
    for polypart in values(poly.particles)
        could_contact(polypart, part) || continue
        overlap(polypart, part) && return true, nothing

        # TODO: we could only iterate over open binding sites
        for b1 in bindingsites(polypart), b2 in bindingsites(part)
            istouching(b1, b2; kwargs...) || continue

            interacting = intmat[color(b1), color(b2)]
            if !allow_noninteracting && !interacting
                return true, nothing
            end
            if !allow_misaligned && !isaligned(b1, b2; kwargs...)
                return true, nothing
            end
            push!(contacts, (b1.vertices, b2.vertices))
        end
    end
    return false, contacts 
end

function select_jth_compatible_siteloc(poly::Polyform, j::Integer)
    i = 0
    for s in 1:nopensites(poly)
        site = open_bindingsite(poly, s)
        for siteloc in compatible_sitelocs(assemblysystem(poly), color(site))
            i += 1
            if i == j
                return site, siteloc
            end
        end
    end
    return nothing, nothing
end

function raise!(poly::Polyform, j::Integer; kwargs...)
    #TODO: handle n==0
    site, siteloc = select_jth_compatible_siteloc(poly::Polyform, j::Integer)
    isnothing(site) && return nothing
    species_index, site_index = siteloc

    particle_species = species(assemblysystem(poly), species_index)
    leading_vertex = nv(graphrep(poly)) + 1
    attached_particle = attach(particle_species => site_index, site; leading_vertex)

    overlap, contacting_vertices = overlap_and_contacts(poly, attached_particle; kwargs...)
    overlap && return missing

    poly.particles[leading_vertex] = attached_particle

    g_attach = graphrep(species(attached_particle))
    add_vertices!(graphrep(poly); vertex_labels=labels(g_attach))
    for (;src, dst) in edges(g_attach)
        add_edge!(graphrep(poly), src + leading_vertex - 1, dst + leading_vertex - 1)
    end

    for (vs1, vs2) in contacting_vertices
        for (v1, v2) in zip(vs1, vs2)
            v1 = invperm(canonical_vertices(poly))[v1]
            add_edge!(graphrep(poly), v1, v2)
            add_edge!(graphrep(poly), v2, v1)
        end
    end

    append!(poly.canonvs, graphvertices(attached_particle))
    perm, autg = nauty(graphrep(poly); canonize=true)
    poly.canonvs .= @view poly.canonvs[perm]
    poly.sigma = convert(Int, autg.n)

    return poly
end

function lower!(poly::Polyform)
    n = nparticles(poly)
    if n == 0 
        return nothing
    elseif n == 1
        pop!(poly.particles, only(keys(poly.particles)))
        resize!(poly.canonvs, 0)
        rem_vertices!(graphrep(poly), vertices(graphrep(poly)))
        poly.sigma = 1
        return poly
    end

    # println("LOWER")
    for v in Iterators.reverse(canonical_vertices(poly))
        # @show canonical_vertices(poly)

        # find the leading vertex corresponding to v
        for k in v:-1:1
            is_leadingvertex(poly, k) || continue
            v = k
            break
        end
        # @show v

        part = particle(poly, v)
        vs0 = graphvertices(part)
        # vs = canonical_vertices(poly)[vs0]
        vs = sort!(invperm(canonical_vertices(poly))[vs0])
        # @show vs

        # @show graphrep(poly).graphset
        # for p in values(poly.particles)
        #     println(leading_vertex(p) => p.pose)
        # end
        are_cutvertices(graphrep(poly), vs) && continue
        # println("##########")
        # @show vs0
        # @show vs

        pop!(poly.particles, leading_vertex(part))
        rem_vertices!(graphrep(poly), vs)
        deleteat!(poly.canonvs, vs)

        vertex_shift(v) = sum(x -> x <= v, vs0)
        @views poly.canonvs .-= vertex_shift.(poly.canonvs)

        perm, autg = nauty(graphrep(poly); canonize=true)
        poly.canonvs .= @view perm[poly.canonvs]
        poly.sigma = convert(Int, autg.n)
        return poly
    end
    error()
end

function attach(species_and_site::Pair{<:ParticleSpecies, <:Integer}, at::BindingSite; leading_vertex)
    particle_species, site_index = species_and_site
    site = bindingsite(particle_species, site_index)
    # We want the pose of bindingsite(particle, site_index) to be equal to the pose of at + the standard_offset
    particle_pose = inv(site.pose) * orientation_offset(at)
    return Particle(particle_species, particle_pose; leading_vertex)
end

function render!(ax, p::Polyform; kwargs...)
    out = for part in values(p.particles)
        render!(ax, part; kwargs...)
    end
    return out
end
function render(p::Polyform; hidedecorations=true, kwargs...)
    f = Figure()
    ax = Axis(f[1, 1], aspect=DataAspect())
    hidedecorations && hidedecorations!(ax)
    render!(ax, p; kwargs...)
    return f
end