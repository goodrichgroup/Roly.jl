"""
    PolygonParticleSpecies{F}

A 2D convex regular polygon with `n` edges and one binding site per edge.
"""
struct PolygonParticleSpecies{F,B<:BindingSite} <: ParticleSpecies{2,B}
    g::NautyDiGraph
    sites::Vector{B}
    corners::Vector{SVector{2,F}}
    rmin::F
    rmax::F
    isregular::Bool
    skin::F
end

"""
    PolygonParticleSpecies(n, a=1.0; colors=1:n, labels=colors)

Construct a regular `n`-gon with edge length `a`. Each edge carries one binding site.

`colors` assigns interaction colors to the binding sites (used by the interaction matrix to
specify which sites bond). `labels` assigns symmetry labels to the graph vertices that determine
graph isomorphism: sites with equal labels are treated as equivalent under the particle's
symmetry group, which affects both isomorphism detection and the symmetry number. By default
`labels=colors`, so distinct colors imply distinct labels.
To treat all sites as equivalent (e.g. for enumerating unlabeled polyforms), pass
`labels=fill(1, n)` explicitly — note this also sets `symmetrynumber` to `n`.
"""
function PolygonParticleSpecies(n::Integer, a::F=1.0; colors=1:n, labels=colors) where {F<:Real}
    r_in = convert(F, 0.5a * cot(π / n))
    r_out = convert(F, 0.5a * csc(π / n))

    tol = sqrt(eps(F)) * r_out
    g, ranges = cycleencoding(n; labels)
    sites = BindingSite{Pose{2,F,Angle2d{F}},F}[]
    for (i, c) in enumerate(colors)
        ψ = Angle2d{F}(-F(π) * (1 / 2 + 2 / n * (i - 1)))
        x = SVector{2,F}(pol2cart(r_in, rotation_angle(ψ)))
        push!(sites, BindingSite(Pose(x, ψ), c, ranges[i], tol, tol / r_in, 1))
    end

    corners = [
        SVector{2,F}(r_out * cos(-π / 2 - (2k - 1) * π / n), r_out * sin(-π / 2 - (2k - 1) * π / n)) for k in 1:n
    ]
    return _check_encoding(
        PolygonParticleSpecies{F,eltype(sites)}(g, sites, corners, r_in, r_out, true, tol)
    )
end

function Base.show(io::Core.IO, ps::PolygonParticleSpecies)
    return print(io, "$(dimension(ps))d PolygonParticleSpecies with $(nsites(ps)) sites")
end

function Base.copy(pps::PolygonParticleSpecies)
    return PolygonParticleSpecies(
        copy(pps.g), copy(pps.sites), copy(pps.corners), pps.rmin, pps.rmax, pps.isregular, pps.skin
    )
end


graphrep(p::PolygonParticleSpecies) = p.g
nsites(p::PolygonParticleSpecies) = length(p.sites)
bindingsites(p::PolygonParticleSpecies, i::Integer) = p.sites[i]
function setcolors!(p::PolygonParticleSpecies, colors::AbstractVector{<:Integer})
    length(colors) != nsites(p) && throw(ArgumentError("incorrect number of colors"))
    for k in eachindex(p.sites)
        s = p.sites[k]
        p.sites[k] = BindingSite(s.pose, colors[k], s.vertices, s.touching_tolerance, s.alignment_tolerance, s.gauge)
    end
    return nothing
end

function can_skip_overlap_check(p1::SpeciesAndPose{<:PolygonParticleSpecies}, p2::SpeciesAndPose{<:PolygonParticleSpecies})
    spcs1, pose1 = p1
    spcs2, pose2 = p2

    norm(pose1.x - pose2.x) >= spcs1.rmax + spcs2.rmax && return true
    spcs1.isregular && spcs2.isregular || return false

    n = nsites(spcs1)
    n in (3, 4, 6) && nsites(spcs2) == n || return false

    sym_angle = 2π / n
    θ_rel = mod(rotation_angle(pose2.psi) - rotation_angle(pose1.psi), sym_angle)
    tol = (spcs1.skin + spcs2.skin) / spcs1.rmin
    face_to_face = mod(π, sym_angle)
    return θ_rel < tol || θ_rel > sym_angle - tol || abs(θ_rel - face_to_face) < tol
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
            (hi2 < lo1 + skin_sum || hi1 < lo2 + skin_sum) && return false
        end
    end
    return true
end


function overlap(p1::SpeciesAndPose{<:PolygonParticleSpecies}, p2::SpeciesAndPose{<:PolygonParticleSpecies}; kwargs...)
    can_skip_overlap_check(p1, p2) || return _SAT_overlap(p1, p2)
    spcs1, pose1 = p1
    spcs2, pose2 = p2
    return norm(pose1.x - pose2.x) < (spcs1.rmin + spcs2.rmin) - (spcs1.skin + spcs2.skin)
end

bounding_radius(ps::PolygonParticleSpecies) = ps.rmax

"""
    UnitTriangle

An equilateral triangle with unit-length edges and one binding site per edge.
"""
const UnitTriangle = PolygonParticleSpecies(3)

"""
    UnitSquare

A square with unit-length edges and one binding site per edge.
"""
const UnitSquare = PolygonParticleSpecies(4)

"""
    UnitHexagon

A regular hexagon with unit-length edges and one binding site per edge.
"""
const UnitHexagon = PolygonParticleSpecies(6)
