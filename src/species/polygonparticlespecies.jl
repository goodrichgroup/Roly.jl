struct PolygonParticleSpecies{F,B<:BindingSite} <: ParticleSpecies{2,F,Pose{2,F,Angle2d{F}}}
    g::NautyDiGraph
    sites::Vector{B}
    corners::Vector{SVector{2,F}}
    rmin::F
    rmax::F
    isregular::Bool
    skin::F
end
function PolygonParticleSpecies(n::Integer, a::F=1.0; colors=1:n, labels=colors) where {F<:Real}
    # TODO: this should have a constructor that takes polygon vertices and automatically creates sites on each polygon side
    r_in = convert(F, 0.5a * cot(π / n))
    r_out = convert(F, 0.5a * csc(π / n))

    tol = sqrt(eps(F)) * r_out
    sites = BindingSite{Pose{2,F,Angle2d{F}},F}[]
    for (i, (c, l)) in enumerate(zip(colors, labels))
        ψ = Angle2d{F}(-F(π) * (1 / 2 + 2 / n * (i - 1)))
        x = SVector{2,F}(pol2cart(r_in, rotation_angle(ψ)))
        push!(sites, BindingSite(Pose(x, ψ), c, i:i, tol, tol / r_in))
    end

    g = NautyDiGraph(cycle_digraph(n); vertex_labels=labels)
    corners = [
        SVector{2,F}(r_out * cos(-π / 2 - (2k - 1) * π / n), r_out * sin(-π / 2 - (2k - 1) * π / n)) for k in 1:n
    ]
    return PolygonParticleSpecies{F,eltype(sites)}(g, sites, corners, r_in, r_out, true, tol)
end

function Base.show(io::Core.IO, ps::PolygonParticleSpecies)
    return print(io, "$(dimension(ps))d PolygonParticleSpecies with $(nsites(ps)) sites")
end

function Base.copy(pps::PolygonParticleSpecies)
    return PolygonParticleSpecies(
        copy(pps.g), copy(pps.sites), copy(pps.corners), pps.rmin, pps.rmax, pps.isregular, pps.skin
    )
end

dimension(::PolygonParticleSpecies) = 2

graphrep(p::PolygonParticleSpecies) = p.g
nsites(p::PolygonParticleSpecies) = length(p.sites)
bindingsites(p::PolygonParticleSpecies, i::Integer) = p.sites[i]
function setcolors!(p::PolygonParticleSpecies, colors::AbstractVector{<:Integer})
    length(colors) != nsites(p) && throw(ArgumentError("incorrect number of colors"))
    for k in eachindex(p.sites)
        s = p.sites[k]
        p.sites[k] = BindingSite(s.pose, colors[k], s.vertices, s.touching_tolerance, s.alignment_tolerance)
    end
    return nothing
end

function can_skip_overlap_check(p1::PolygonParticleSpecies, p2::PolygonParticleSpecies)
    p1.isregular && p2.isregular || return false
    n = nsites(p1)
    n in (3, 4, 6) || return false
    return nsites(p2) == n
end

isconvex(::PolygonParticleSpecies) = true

function _SAT_overlap(p1::SpeciesAndPose{<:PolygonParticleSpecies}, p2::SpeciesAndPose{<:PolygonParticleSpecies})
    spcs1, pose1 = p1
    spcs2, pose2 = p2
    isconvex(spcs1) && isconvex(spcs2) || error("SAT overlap check requires convex polygons")

    # The SAT algorithm checks for intersections of convex polygons by exploiting a simple theorem:
    # For any nonintersecting convex shapes, there exists an axis in space on which the projections
    # of the shapes are nonoverlapping. Therefore, we project all polygon corners on axes given by
    # the surface normals of the polygons. As soon as we identify an axis where projections do not overlap
    # we know that the polygons are not intersecting.
    skin_sum = spcs1.skin + spcs2.skin
    corners1, corners2 = spcs1.corners, spcs2.corners
    for (corners, pose) in ((corners1, pose1), (corners2, pose2))
        n = length(corners)
        for i in 1:n
            edge = pose.psi * (corners[mod1(i + 1, n)] - corners[i])
            normal = SVector(-edge[2], edge[1])
            lo1, hi1 = extrema(dot(normal, pose1 * c) for c in corners1)
            lo2, hi2 = extrema(dot(normal, pose2 * c) for c in corners2)
            (hi2 < lo1 - skin_sum || hi1 < lo2 - skin_sum) && return false
        end
    end
    return true
end

function could_contact(
    p1::SpeciesAndPose{<:PolygonParticleSpecies}, p2::SpeciesAndPose{<:PolygonParticleSpecies}; kwargs...
)
    spcs1, pose1 = p1
    spcs2, pose2 = p2
    return norm(pose1.x - pose2.x) < spcs1.rmax + spcs2.rmax
end
function overlap(p1::SpeciesAndPose{<:PolygonParticleSpecies}, p2::SpeciesAndPose{<:PolygonParticleSpecies}; kwargs...)
    spcs1, pose1 = p1
    spcs2, pose2 = p2
    skin_sum = spcs1.skin + spcs2.skin
    if can_skip_overlap_check(spcs1, spcs2)
        return norm(pose1.x - pose2.x) < (spcs1.rmin + spcs2.rmin) - skin_sum
    else
        return _SAT_overlap(p1, p2)
    end
end

const UnitTriangle = PolygonParticleSpecies(3)
const UnitSquare = PolygonParticleSpecies(4)
const UnitHexagon = PolygonParticleSpecies(6)
