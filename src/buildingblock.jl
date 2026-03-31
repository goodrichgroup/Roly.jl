struct BindingSiteCollection{D,F,R<:Rotation{D,F}}
    xs::Vector{Point{D,F}}
    ψs::Vector{R}
end
function PolygonBindingSites(n::Integer, a::F=1.0) where {F<:Real}
    # TODO: this should have a constructor that takes polygon vertices and automatically creates sites on each polygon side
    r_in = convert(F, 0.5a * cot(π / n))
    ψs = [Angle2d{F}(-F(π) * (1/2 + 2/n * i)) for i in 0:(n-1)]
    xs = [Point{2,F}(pol2cart(r_in, rotation_angle(ψ))) for ψ in ψs]
    return BindingSiteCollection{2,F,Angle2d{F}}(xs, ψs)
end

standard_offset(::R) where {R<:Rotation{2}} = R(π)
standard_offset(::Rotation{3,F}) where {F} = RotZYZ{F}(π, 0, 0)

positiontype(::BindingSiteCollection{D,F}) where {D,F} = Point{D,F}
orientationtype(::BindingSiteCollection{D,F,R}) where {D,F,R} = R
positiontype(::Type{<:BindingSiteCollection{D,F}}) where {D,F} = Point{D,F}
orientationtype(::Type{<:BindingSiteCollection{D,F,R}}) where {D,F,R} = R

Base.show(io::Core.IO, bs::BindingSiteCollection) = print(io, "$(dimension(bs))d BindigSiteCollection with $(nsites(bs)) sites")

dimension(::BindingSiteCollection{D}) where {D} = D
nsites(bs::BindingSiteCollection) = length(bs.xs)

function Base.copy(bs::BindingSiteCollection)
    return BindingSiteCollection(copy(bs.xs), copy(bs.ψs))
end

struct BuildingBlock{D,F,G<:AbstractGeometry{F},S<:BindingSiteCollection{D,F}}
    anatomy::NautyDiGraph
    sites::S
    siteof::Vector{Int} # sites[v] is the binding site vertex v corresponds to; is generally 1-1 with vertices(anatomy) unless higher symmetry is present #TODO: can prob optimize away in most cases
    geometry::G
    #TODO: enfore site numbering from 1 to nsites
    #TODO: enforce that the location of binding sites is between rmin and rmax as defined by the geometry, or better yet, add make rmin, rmax fields of this object
end
function PolygonBuildingBlock(n::Integer, a::Real=1.0; labels=1:n)
    geom = PolygonGeometry(n, a)
    sites = PolygonBindingSites(n, a)
    siteof = collect(1:n)
    anatomy = NautyDiGraph(cycle_digraph(n); vertex_labels=labels)
    return BuildingBlock{2,typeof(a),typeof(geom),typeof(sites)}(anatomy, sites, siteof, geom)
end
Base.show(io::Core.IO, bb::BuildingBlock) = print(io, "$(dimension(bb))d BuildingBlock with $(nsites(bb)) sites and $(nvertices(bb)) vertices")

anatomy(bb::BuildingBlock) = bb.anatomy
dimension(::BuildingBlock{D}) where {D} = D
vertices(bb::BuildingBlock) = NautyGraphs.vertices(anatomy(bb))
vertices(bb::BuildingBlock, site) = findall(==(site), bb.siteof)
nvertices(bb::BuildingBlock, site) = count(==(site), bb.siteof)
nvertices(bb::BuildingBlock) = nv(anatomy(bb))

nlabels(bb::BuildingBlock) = length(unique(labels(anatomy(bb))))
nsites(bb::BuildingBlock) = nsites(bb.sites)
sites(bb::BuildingBlock) = bb.sites
site(bb::BuildingBlock, v) = sites(bb)[bb.siteof[v]]

numtype(::BuildingBlock{D,F}) where {D,F} = F
numtype(::Type{<:BuildingBlock{D,F}}) where {D,F} = F
positiontype(::BuildingBlock{D,F,G,S}) where {D,F,G,S} = positiontype(S)
positiontype(::Type{<:BuildingBlock{D,F,G,S}}) where {D,F,G,S} = positiontype(S)
orientationtype(::BuildingBlock{D,F,G,S}) where {D,F,G,S} = orientationtype(S)
orientationtype(::Type{<:BuildingBlock{D,F,G,S}}) where {D,F,G,S} = orientationtype(S)

function Base.copy(bb::BB) where {BB}
    return BB(copy(bb.anatomy), copy(bb.sites), copy(bb.siteof), copy(bb.geometry))
end

const UnitTriangle = PolygonBuildingBlock(3, Float32(1))
const UnitSquare = PolygonBuildingBlock(4, Float32(1))
const UnitHexagon = PolygonBuildingBlock(6, Float32(1))