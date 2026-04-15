struct BindingSiteCollection{D,F,R<:Rotation{D,F}} #TODO This can be completely static
    xs::Vector{Point{D,F}} # make both of these staticarrays
    ψs::Vector{R}
end
function PolygonBindingSites(n::Integer, a::F=1.0) where {F<:Real}
    # TODO: this should have a constructor that takes polygon vertices and automatically creates sites on each polygon side
    r_in = convert(F, 0.5a * cot(π / n))
    ψs = [Angle2d{F}(-F(π) * (1/2 + 2/n * i)) for i in 0:(n-1)]
    xs = [Point{2,F}(pol2cart(r_in, rotation_angle(ψ))) for ψ in ψs]
    return BindingSiteCollection{2,F,Angle2d{F}}(xs, ψs)
end

@inline standard_offset(::R) where {R<:Rotation{2}} = R(π)
@inline standard_offset(::Rotation{3,F}) where {F} = RotXYZ{F}(0, 0, π)
@inline positiontype(::BindingSiteCollection{D,F}) where {D,F} = Point{D,F}
@inline orientationtype(::BindingSiteCollection{D,F,R}) where {D,F,R} = R
@inline positiontype(::Type{<:BindingSiteCollection{D,F}}) where {D,F} = Point{D,F}
@inline orientationtype(::Type{<:BindingSiteCollection{D,F,R}}) where {D,F,R} = R
@inline dimension(::BindingSiteCollection{D}) where {D} = D
@inline nsites(bs::BindingSiteCollection) = length(bs.xs)

Base.show(io::Core.IO, bs::BindingSiteCollection) = print(io, "$(dimension(bs))d BindigSiteCollection with $(nsites(bs)) sites")

function Base.copy(bs::BindingSiteCollection)
    return BindingSiteCollection(copy(bs.xs), copy(bs.ψs))
end

struct BuildingBlock{D,F,G<:AbstractGeometry{F},S<:BindingSiteCollection{D,F}}
    anatomy::NautyDiGraph
    sites::S
    geometry::G
    σ::Int
    _siteof::Vector{Int} # sites[v] is the binding site that vertex v corresponds to; is generally 1-1 with vertices(anatomy) unless higher symmetry is present #TODO: can prob optimize away in most cases
    #TODO: enfore site numbering from 1 to nsites
    #TODO: enforce that the location of binding sites is between rmin and rmax as defined by the geometry, or better yet, add make rmin, rmax fields of this object
end
function PolygonBuildingBlock(n::Integer, a::Real=1.0; labels=1:n)
    geom = PolygonGeometry(n, a)
    sites = PolygonBindingSites(n, a)
    siteof = collect(1:n)
    anatomy = NautyDiGraph(cycle_digraph(n); vertex_labels=labels)
    _, autgrp = nauty(anatomy)
    σ = autgrp.n
    return BuildingBlock{2,typeof(a),typeof(geom),typeof(sites)}(anatomy, sites, geom, σ, siteof)
end
Base.show(io::Core.IO, bb::BuildingBlock) = print(io, "$(dimension(bb))d BuildingBlock with $(nsites(bb)) sites and $(nvertices(bb)) vertices")

@inline anatomy(bb::BuildingBlock) = bb.anatomy
@inline dimension(::BuildingBlock{D}) where {D} = D
@inline sitevertices(bb::BuildingBlock) = NautyGraphs.vertices(anatomy(bb))
@inline sitevertices(bb::BuildingBlock, site) = findall(==(site), bb._siteof)
@inline nvertices(bb::BuildingBlock, site) = count(==(site), bb._siteof)
@inline nvertices(bb::BuildingBlock) = nv(anatomy(bb))
@inline nlabels(bb::BuildingBlock) = length(unique(labels(anatomy(bb))))
@inline nsites(bb::BuildingBlock) = nsites(bb.sites)
@inline sites(bb::BuildingBlock) = bb.sites
@inline site(bb::BuildingBlock, v) = sites(bb)[bb._siteof[v]]
@inline symmetrynumber(bb::BuildingBlock) = bb.σ
@inline numtype(::BuildingBlock{D,F}) where {D,F} = F
@inline numtype(::Type{<:BuildingBlock{D,F}}) where {D,F} = F
@inline positiontype(::BuildingBlock{D,F,G,S}) where {D,F,G,S} = positiontype(S)
@inline positiontype(::Type{<:BuildingBlock{D,F,G,S}}) where {D,F,G,S} = positiontype(S)
@inline orientationtype(::BuildingBlock{D,F,G,S}) where {D,F,G,S} = orientationtype(S)
@inline orientationtype(::Type{<:BuildingBlock{D,F,G,S}}) where {D,F,G,S} = orientationtype(S)

function Base.copy(bb::BB) where {BB}
    return BB(copy(bb.anatomy), copy(bb.sites), copy(bb.geometry), bb.σ, copy(bb._siteof))
end



const UnitTriangle = PolygonBuildingBlock(3, Float32(1))
const UnitSquare = PolygonBuildingBlock(4, Float32(1))
const UnitHexagon = PolygonBuildingBlock(6, Float32(1))