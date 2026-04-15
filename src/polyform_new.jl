mutable struct PolyformGraphRep{G<:AbstractNautyGraph}
    const g::G
    sigma::Int
    canonvs::Vector{Int}
    function PolyformGraphRep(g::AbstractNautyGraph)
        g = copy(g)
        canonize!(g)
        return new{typeof(g)}(g, 0, collect(vertices(g)))
    end
    function PolyformGraphRep(args...; kwargs...)
        g = NautyDiGraph(args...; kwargs...)
        canonize!(g)
        new{NautyDiGraph}(g, 0, collect(vertices(g)))
    end
end

NautyGraphs.labels(pgr::PolyformGraphRep) = NautyGraphs.labels(pgr.g)
NautyGraphs.edges(pgr::PolyformGraphRep) = NautyGraphs.edges(pgr.g)
NautyGraphs.vertices(pgr::PolyformGraphRep) = NautyGraphs.vertices(pgr.g)
NautyGraphs.nv(pgr::PolyformGraphRep) = NautyGraphs.nv(pgr.g)
NautyGraphs.ne(pgr::PolyformGraphRep) = NautyGraphs.ne(pgr.g)

"""
    canonical_vertices(pgr::PolyformGraphRep)

Return the list of graphrep vertices in canonical order.
"""
function canonical_vertices(pgr::PolyformGraphRep)
    return pgr.canonvs
end


"""
    Polyform{D,T,F} where {D,T,F}

`Polyform` is a struct containing all information required to instantiate a
polyform, which here is defined as an aggregate of building blocks that are connected
(bound) together at various binding sites. In a self-assembly context, a polyform is
often referred to a "self-assembled structure".
"""
struct Polyform{D,P<:Particle,S<:AssemblySystem,PGR<:PolyformGraphRep}
    graphrep::PGR
    particles::Dict{Int,P}
    assemblysystem::S
end
function Polyform(sys::AssemblySystem{D}) where {D}
    P = posetype(sys)
    pgr = PolyformGraphRep(0)
    return Polyform{D,Particle{P,typeof(sys)},typeof(sys),typeof(pgr)}(pgr, Dict(), sys)
end
function Polyform(sys::AssemblySystem{D}, i::Integer) where {D}
    P = posetype(sys)
    ps = species(sys, i)
    pgr = PolyformGraphRep(graphrep(ps))
    return Polyform{D,Particle{P,typeof(ps)},typeof(sys),typeof(pgr)}(pgr, Dict(1=>Particle(ps; leading_vertex=1)), sys)
end

function Base.show(io::Core.IO, p::Polyform{D}) where {D}
    print(io, "Polyform{$D}[n=$(nparticles(p))]")
end
Base.show(io::Core.IO, ::Type{Polyform{D}}) where {D} = print(io, "Polyform{$D}")

@inline particle(p::Polyform, v::Integer) = get(p.particles, v, nothing)
@inline nparticles(p::Polyform) = length(p.particles)
@inline nsites(p::Polyform) = sum(nsites, values(p.particles))
@inline graphrep(p::Polyform) = p.graphrep
@inline symmetrynumber(p::Polyform) = graphrep(p).sigma
@inline assemblysystem(p::Polyform) = p.assemblysystem

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
    for v in canonical_vertices(graphrep(p))
        prtcl = particle(p, v)
        isnothing(prtcl) && continue

        for j in 1:nsites(prtcl)
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

# function contacting_sites(p1::SpeciesAndPose{<:PolygonParticleSpecies}, p2::SpeciesAndPose{<:PolygonParticleSpecies}; buffer=0.1)
#     spcs1, pose1 = p1
#     spcs2, pose2 = p2

#     for b1 in bindingsites(spcs1), b2 in bindingsites(spcs2)
#         p1 = pose1 * b1.pose
#         p2 = pose2 * b2.pose
#         x_overlap = isapprox(p1.x, p2.x; kwargs...)
#         if x_overlap
#             interacting = intmat[si, siteloc2index(, j)]
#             if !allow_noninteracting && !interacting
#                 return false, contacts
#             end

#             ψ_overlap = isapprox(ψs1[j], ψs2[j]; kwargs...)
#             if !allow_misaligned && !ψ_overlap
#                 return false, contacts
#             end
#             push!(contacts, (i, j))
#         end
#     end
#     return 
# end

function overlap_and_contacts(poly::Polyform, part::Particle; allow_noninteracting=false, allow_misaligned=false, kwargs...)
    intmat = interactionmatrix(assemblysystem(poly))
    return poly, part

    contacts = NTuple{2,Int}[]
    for polypart in values(poly.particles)
        could_contact(polypart, part) || continue
        # overlap(polypart, part) && return true, nothing

        for b1 in bindingsites(polypart), b2 in bindingsites(part)
            istouching(b1, b2; kwargs...) || continue

            interacting = intmat[color(b1), color(b2)]
            if !allow_noninteracting && !interacting
                return true, nothing
            end
            if !allow_misaligned && !isaligned(b1, b2; kwargs...)
                return true, nothing
            end
            push!(contacts, (i, j))
        end
    end
    return false, contacts 
end

function raise!(p::Polyform, j::Integer; kwargs...)
    # attempt to connect the `k`th possible neighbor to the `s`th binding site
    n = nsites(p)
    s, k = mod1(j, n), cld(j, n)

    site = bindingsite(p, s)
    leading_vertex = nv(graphrep(p)) + 1

    siteloc = take_nth(compatible_sitelocs(assemblysystem(p), color(site)), k)
    isnothing(siteloc) && return missing
    species_index, site_index = siteloc

    particle_species = species(assemblysystem(p), species_index)
    attached_particle = attach(particle_species => site_index, site; leading_vertex)

    overlap, contacting_sites = overlap_and_contacts(p, attached_particle; kwargs...)
    return overlap, contacting_sites
    overlap && return missing

        # add edges and push particle to dict
    return nothing
end
function lower!(p::Polyform)
    n = nparticles(p)
    n == 0 && return false

    for v in canonical_vertices(graphrep(p))
        vs = graphvertices(particle(p, v))
        #TODO: I think vertices must allocate. Maybe we can perform the check without allocating it first?
        are_cutvertices(p.anatomy, vs) && continue

        # delete vertices and particle
    end
end

function attach(species_and_site::Pair{<:ParticleSpecies, <:Integer}, at::BindingSite; leading_vertex)
    particle_species, site_index = species_and_site
    site = bindingsite(particle_species, site_index)
    # We want the pose of bindingsite(particle, site_index) to be equal to the pose of at + the standard_offset
    particle_pose = orientation_offset(at) * inv(site.pose)
    return Particle(particle_species, particle_pose; leading_vertex)
end