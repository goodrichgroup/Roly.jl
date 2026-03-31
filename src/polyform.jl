"""
    Polyform{D,T,F} where {D,T,F}

`Polyform` is a struct containing all information required to instantiate a
polyform, which here is defined as an aggregate of building blocks that are connected
(bound) together at various binding sites. In a self-assembly context, a polyform is
often referred to a "self-assembled structure".
"""

mutable struct PolyformAnatomy{G<:AbstractNautyGraph}
    const anatomy::G            # always keep this canonized (?)
    σ::Int
    function PolyformAnatomy(anatomy::AbstractNautyGraph)
        anatomy = copy(anatomy)
        canonize!(anatomy)
        return new{typeof(anatomy)}(anatomy, 0)
    end
    function PolyformAnatomy(args...; kwargs...)
        anatomy = NautyDiGraph(args...; kwargs...)
        canonize!(anatomy)
        new{NautyDiGraph}(anatomy, 0)
    end
end

Base.show(io::Core.IO, ::Type{PolyformAnatomy{G}}) where {G} = print(io, "PolyformAnatomy[$G]")

Base.copy(pa::PolyformAnatomy) = PolyformAnatomy(copy(pa.anatomy), pa.σ)
function Base.copy!(dest::PolyformAnatomy{G}, src::PolyformAnatomy{G}) where {G}
    copy!(dest.anatomy, src.anatomy)
    dest.σ = src.σ
    return dest
end

function canonid(pa::PolyformAnatomy)
    return canonical_id(pa.anatomy)
end
Base.hash(pa::PolyformAnatomy, h::UInt=zero(Uint)) = hash(pa.anatomy, h)
Base.:(==)(p::PolyformAnatomy, h::PolyformAnatomy) = p.anatomy ≃ h.anatomy
Graphs.vertices(pa::PolyformAnatomy) = vertices(pa.anatomy)
Graphs.edges(pa::PolyformAnatomy) = edges(pa.anatomy)
sitelabels(pa::PolyformAnatomy) = NautyGraphs.labels(pa.anatomy)

function symmetrynumber(pa::PolyformAnatomy) 
    if iszero(pa.σ)
        _, autgrp = nauty(pa)
        pa.σ = autgrp.n
    end
    return pa.σ
end

struct Polyform{D,F<:AbstractFloat,R<:Rotation{D},PA<:PolyformAnatomy,AS<:AssemblySystem}
    anatomy::PA
    xs::Vector{Point{D,F}}
    ψs::Vector{R}
    assemblysystem::AS
    _bonds::Vector{Int}        # bonds[v] is the vertex v is bound to (or zero if unbound)
    _nparticles::Int
end
function Polyform(sys::AssemblySystem{D}) where {D}
    P = positiontype(sys)
    R = orientationtype(sys)
    A = PolyformAnatomy(0)
    return Polyform{D,numtype(sys),R,typeof(A),typeof(sys)}(A, P[], R[], sys, Int[], 0)
end
function Polyform(sys::AssemblySystem{D}, buildingblock::BuildingBlock{D}) where {D}
    R = orientationtype(sys)
    A = PolyformAnatomy(anatomy(buildingblock))
    return Polyform{D,numtype(sys),R,typeof(A),typeof(sys)}(A, copy(sites(buildingblock).xs), copy(sites(buildingblock).ψs), sys, zeros(Int, nvertices(buildingblock)), 1)
end
Polyform(sys::AssemblySystem, i::Integer) = Polyform(sys, buildingblocks(sys)[i])

anatomy(p::Polyform) = p.anatomy
sitelabels(p::Polyform) = sitelabels(anatomy(p))
symmetrynumber(p::Polyform) = anatomy(p).σ
assemblysystem(p::Polyform) = p.assemblysystem
positions(p::Polyform) = p.xs
orientations(p::Polyform) = p.ψs

Graphs.vertices(p::Polyform) = vertices(anatomy(p))

function interior_edges(p::Polyform)
    return _filter_edges_by_species(p, Val(true))
end
function exterior_edges(p::Polyform)
    return _filter_edges_by_species(p, Val(false))
end
function _filter_edges_by_species(p::Polyform, ::Val{same_species}) where same_species
    sys = assemblysystem(p)
    labs = sitelabels(p)
    es = edges(anatomy(p))
    interior_edges = Iterators.filter(es) do (; src, dst)
        src_spcs = label2species(sys, labs[src])
        dst_spcs = label2species(sys, labs[dst])
        return same_species ? src_spcs == dst_spcs : src_spcs != dst_spcs
    end
    return interior_edges
end

nparticles(p::Polyform) = p._nparticles
nvertices(p::Polyform) = nv(p.anatomy)
function Base.show(io::Core.IO, p::Polyform{D,T,F}) where {D,T,F}
    print(io, "Polyform{$D,$T,$F}[n=$(nparticles(p))]")
end
Base.show(io::Core.IO, ::Type{Polyform{D,T,F}}) where {D,T,F} = print(io, "Polyform{$D,$T,$F}")
dimension(::Polyform{D}) where {D} = D
dimension(::Type{Polyform{D}}) where {D} = D

# function canonize!(p::Polyform{D,T}) where {D,T}
#     canonperm, autg = nauty(p.anatomy; canonize=true)
#     p.canonical_order .= canonperm
#     p.σ = convert(T, autg.n)
#     return
# end

# function canonid(p::Polyform)
#     return canonid(anatomy(p))
# end
# is_isomorphic(p::Polyform, h::Polyform) = p ≃ h
# ≃(p::Polyform, h::Polyform) = NautyGraphs.is_isomorphic(p.anatomy, h.anatomy)


# function Base.copy(p::Polyform)
#     return Polyform(copy(p.anatomy), copy(p.bond_partners), copy(p.canonical_order), copy(p.species), copy(p.xs), copy(p.ψs), p.σ)
# end
# function Base.copy!(dest::Polyform{D,T,F}, src::Polyform{D,T,F}) where {D,T,F}
#     copy!(dest.anatomy, src.anatomy)
#     copy!(dest.bond_partners, src.bond_partners)
#     copy!(dest.canonical_order, src.canonical_order)
#     copy!(dest.species, src.species)
#     copy!(dest.xs, src.xs)
#     copy!(dest.ψs, src.ψs)
#     dest.σ = src.σ
#     return dest
# end


# function particle2vertex(p::Polyform, assembly_system, particle::Integer, site::Integer)
#     @assert particle <= nparticles(p)
#     @assert site <= nsites(assembly_system.geometries[p.species[particle]])
#     vps = (vertices_per_site(assembly_system.geometries[i]) for i in @view p.species[1:particle])
#     return 1 + sum(sum(nvs) for nvs in Iterators.take(vps, particle-1); init=0) + sum(@view last(vps)[1:site-1])
# end
# function particle2vertex(p::Polyform, assembly_system, particle::Integer)
#     @assert particle <= nparticles(p)
#     vps = (vertices_per_site(assembly_system.geometries[i]) for i in @view p.species[1:particle])
#     vertices = sum(sum(nvs) for nvs in Iterators.take(vps, particle-1); init=0) .+ cumsum((1, last(vps)[1:end-1]...))
#     return vertices
# end

# function particle2multivertex(p::Polyform, assembly_system, particle::Integer, site::Integer)
#     @assert particle <= nparticles(p)
#     @assert site <= nsites(assembly_system.geometries[p.species[particle]])
#     vps = (vertices_per_site(assembly_system.geometries[i]) for i in @view p.species[1:particle])
#     start_vertex = 1 + sum(sum(nvs) for nvs in Iterators.take(vps, particle-1); init=0) + sum(@view last(vps)[1:site-1])
#     end_vertex = start_vertex + last(vps)[site] - 1
#     return start_vertex:end_vertex
# end
# function particle2multivertex(p::Polyform, assembly_system, particle::Integer)
#     @assert particle <= nparticles(p)
#     vps = (vertices_per_site(assembly_system.geometries[i]) for i in @view p.species[1:particle])
#     start_vertex = 1 + sum(sum(nvs) for nvs in Iterators.take(vps, particle-1); init=0)
#     end_vertex = start_vertex + sum(last(vps)) - 1
#     return start_vertex:end_vertex
# end
# function vertex2particle(p::Polyform, assembly_system, vertex::Integer)
#     @assert vertex <= nvertices(p)
#     vps = (vertices_per_site(assembly_system.geometries[i]) for i in p.species)
#     particle, site = 0, 0
#     for vp in vps
#         particle += 1
#         for v in vp
#             site += 1
#             vertex -= v
#             if vertex <= 0
#                 return particle, site
#             end
#         end
#         site = 0
#     end
#     return
# end
# function vertex2label(p::Polyform, assembly_system, vertex::Integer)
#     particle, site = vertex2particle(p, assembly_system, vertex)
#     spcs = species(p)[particle]
#     return spcssite2label(spcs, site, assembly_system)
# end