struct BindingSiteCollection{D,F}
    xs::Vector{Point{D,F}}        # Displacement vectors from the center to the binding sites
    θs_ref::Vector{Rotation{D,F}} # Rotation necessary to rotate side i into reference position #TODO change this into the relative orientation (its inverse)
    Δθs::Vector{Rotation{D,F}}    # Rotation to be applied to the added particle (assumed to be in reference orientation) once attached to side i
end
function PolygonBindingSites(n::Integer, a::F=1.0) where {F<:Real}
    # TODO: this should have a constructor that takes polygon vertices and automatically creates sites on each polygon side
    r_in = convert(F, 0.5a * cot(π / n))
    θs_ref = [Angle2d{F}(2/n * i) for i in 0:(n-1)]
    Δθs = [Angle2d{F}(1 - 2/n*i) for i in 0:(n-1)]
    rs = r_in * ones(F, n)
    xs = [SVector{2,F}(pol2cart(r, -rotation_angle(θ) - 1/2)) for (r, θ) in zip(rs, θs_ref)]
    return BindingSiteCollection{2,F}(xs, θs_ref, Δθs)
end

Base.show(io::Core.IO, bs::BindingSiteCollection) = print(io, "$(dimension(bs))d BindigSiteCollection with $(nsites(bs)) sites")

dimension(::BindingSiteCollection{D}) where {D} = D
nsites(bs::BindingSiteCollection) = length(bs.xs)

function Base.copy(bs::BindingSiteCollection)
    return BindingSiteCollection(copy(bs.xs), copy(bs.θs_ref), copy(bs.Δθs))
end

struct BuildingBlock{D,G<:AbstractGeometry,B<:BindingSiteCollection}
    anatomy::NautyDiGraph
    sites::B
    siteof::Vector{Int} # sites[v] is the binding site vertex v corresponds to; is generally 1-1 with vertices(anatomy) unless higher symmetry is present #TODO: can prob optimize away in most cases
    geometry::G
    #TODO: enfore site numbering from 1 to nsites
end
function PolygonBuildingBlock(n::Integer, a::Real=1.0; labels=1:n)
    geom = PolygonGeometry(n, a)
    sites = PolygonBindingSites(n, a)
    siteof = collect(1:n)
    anatomy = NautyDiGraph(cycle_digraph(n); vertex_labels=labels)
    return BuildingBlock{2,typeof(geom),typeof(sites)}(anatomy, sites, siteof, geom)
end
Base.show(io::Core.IO, bb::BuildingBlock) = print(io, "$(dimension(bb))d BuildingBlock with $(nsites(bb)) sites and $(nvertices(bb)) vertices")

dimension(::BuildingBlock{D}) where {D} = D
vertices(bb::BuildingBlock, site) = findall(==(site), bb.siteof)
nvertices(bb::BuildingBlock, site) = count(==(site), bb.siteof)
nvertices(bb::BuildingBlock) = nv(bb.anatomy)

nsites(bb::BuildingBlock) = nsites(bb.sites)
site(bb::BuildingBlock, v) = bb.siteof[v]
sites(bb::BuildingBlock) = 1:nsites(bb)

function Base.copy(bb::BB) where {BB}
    return BB(copy(bb.anatomy), copy(bb.sites), copy(bb.siteof), copy(bb.geometry))
end
