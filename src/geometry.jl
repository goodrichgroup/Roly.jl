abstract type AbstractGeometry{F<:AbstractFloat} end

struct PolygonGeometry{F<:AbstractFloat} <: AbstractGeometry{F}
    corners::Vector{Point{2,F}}             # Corners of the Polygon used for overlap checking
    Rmin::F
    Rmax::F
end
function PolygonGeometry(n::Integer, a::F) where {F}
    Rmin = convert(F, 0.5a * cot(π / n))
    Rmax = convert(F, 0.5a * csc(π / n))
    θs = [2/n * i for i in 0:(n-1)]
    corners = [convert.(F, Rmax * [cos(π * (θ + 1/n)), sin(π * (θ + 1/n))]) for θ in θs]
    return PolygonGeometry{F}(corners, Rmin, Rmax)
end
Base.show(io::Core.IO, geom::PolygonGeometry) = print(io, "$(dimension(geom))d PolygonGeometry with $(ncorners(geom)) vertices")

dimension(::PolygonGeometry) = 2
ncorners(geom::PolygonGeometry) = length(geom.corners)

function Base.copy(geom::PolygonGeometry)
    return PolygonGeometry(copy(geom.corners), geom.Rmin, geom.Rmax)
end


# struct PolyhedronGeometry{F<:AbstractFloat} <: AbstractGeometry{F}
#     xs::Vector{Point{3,F}}                # Displacement vectors from the center to the binding sites
#     θs_ref::Vector{Quaternion{F}}         # Rotation necessary to rotate side i into reference position   
#     Δθs::Vector{Quaternion{F}}            # Rotation to be applied to the added particle (assumed to be in reference orientation) once attached to side i 
#     corners::Vector{Point{3,F}}
#     anatomy::NautyDiGraph                 # Oriented graph of the dual polyhedron associated with the bb symmetry group
#     site_vertices::Vector{T}
#     R_min::F
#     R_max::F

#     function PolyhedronGeometry{T}(shape::Symbol, a::F) where {T,F}
#         if shape == :cube
#             xs = a * [Point{3,F}(1, 0, 0),
#                       Point{3,F}(0, 0, -1),
#                       Point{3,F}(0, 1, 0),
#                       Point{3,F}(0, 0, 1),
#                       Point{3,F}(0, -1, 0),
#                       Point{3,F}(-1, 0, 0)]

#             corners = vec([a * SVector(i, j, k) for i in F[-1, 1], j in F[-1, 1], k in F[-1, 1]])
#             θs_ref = [Quaternion{F}(1, 0, 0, 0),
#                       inv(√2) * Quaternion{F}(1, 0, -1, 0),
#                       inv(√2) * Quaternion{F}(1, 0, 0, -1),
#                       inv(√2) * Quaternion{F}(1, 0, 1, 0),
#                       inv(√2) * Quaternion{F}(1, 0, 0, 1),
#                       Quaternion{F}(0, 0, 0, 1)]
#             Δθs = [Quaternion{F}(0, 0, 0, 1),
#                    Quaternion{F}(0, 0, 0, 1) * (inv(√2) * Quaternion{F}(1, 0, -1, 0)),
#                    inv(√2) * Quaternion{F}(1, 0, 0, -1),
#                    Quaternion{F}(0, 0, 0, 1) * (inv(√2) * Quaternion{F}(1, 0, 1, 0)),
#                    inv(√2) * Quaternion{F}(1, 0, 0, 1),
#                    Quaternion{F}(1, 0, 0, 0)]

#             g0 = sparse([0 1 0 0;
#                          0 0 1 0;
#                          0 0 0 1;
#                          1 0 0 0])
#             G = blockdiag(g0, g0, g0, g0, g0, g0)
#             edges = [1, 7], [2, 20], [3, 13], [4, 10], [5, 21], [6, 17], 
#             [8, 9], [11, 16], [12, 22], [14, 19], [15, 23], [18, 24]
#             for e in edges
#                 i, j = e
#                 G[i, j] = G[j, i] = 1
#             end
#             anatomy = NautyDiGraph(G)
#             site_vertices = fill(T(4), 6)

#             R_min = minimum(norm.(xs))
#             R_max = maximum(norm.(corners))
#         else
#             error("Shape $shape not supported.")
#         end
        
#         return new{T,F}(xs, θs_ref, Δθs, corners, anatomy, site_vertices, R_min, R_max)
#     end
# end

# const UnitCubeGeometry = PolyhedronGeometry{DefInt}(:cube, DefFloat(1.))
# #TODO Implement basic 3d shapes: Platonic solids and polygon extrusions

# dimension(::PolyhedronGeometry) = 3

function attachment_offset(site_i::Integer, site_j::Integer,
                           geom_i::AbstractGeometry,
                           geom_j::AbstractGeometry)
    x_i, Δθ = geom_i.xs[site_i], geom_i.Δθs[site_i]
    x_j, θ_ref = geom_j.xs[site_j], geom_j.θs_ref[site_j]

    Δψ = Δθ * θ_ref
    Δx = x_i - rotate(x_j, Δψ)
    return Δx, Δψ
end

# TODO: this function is susceptible to roundoff errors, be very careful with thresh
function face_pairs(Δx::Point{D,F}, ψi::RotationOperator{F}, ψj::RotationOperator{F}, geom_i::AbstractGeometry, geom_j::AbstractGeometry, convex=false, thresh=1e-2) where {D,F}
    pairs = Tuple{Int,Int}[]

    for (i, z_i) in enumerate(geom_i.xs), (j, z_j) in enumerate(geom_j.xs)
        z_i = rotate(z_i, ψi)
        z_j = rotate(z_j, ψj)

        matching = maximum(abs, z_i - z_j - Δx) < F(thresh)
        if matching
            push!(pairs, (i, j))
            if convex
                break
            end
        end
    end
    return pairs
end

function overlap(Δx, ψ_i, ψ_j, geom_i, geom_j; buffer=0.1)
    # Assumes lattice!!
    return norm(Δx) < 2geom_i.R_min - buffer
    # return __check_overlap(geom_i.mesh, transform(geom_j.mesh, Δx, Δψ))
end

function contact_status(Δx::Point, ψi::RotationOperator, ψj::RotationOperator, geom_i::AbstractGeometry, geom_j::AbstractGeometry)
    if norm(Δx) > geom_i.R_max + geom_j.R_max
        return true, Tuple{Int,Int}[]
    end

    if overlap(Δx, ψi, ψj, geom_i, geom_j)
        return false, Tuple{Int,Int}[]
    end

    pairs = face_pairs(Δx, ψi, ψj, geom_i, geom_j)
    return true, pairs
end
