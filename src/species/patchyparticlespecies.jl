"""
    PatchyParticleSpecies{D,F,B}

A `D`-dimensional sphere of radius `r` with binding sites on its surface.
"""
struct PatchyParticleSpecies{D,F,B<:BindingSite} <: ParticleSpecies{D,B}
    g::NautyDiGraph
    sites::Vector{B}
    r::F
    twists::Vector{F}
    skin::F
end

"""
    PatchyParticleSpecies(g::NautyDiGraph, r, patch_positions, patch_twists=zeros(...);
                          colors, vertices, sitesyms, locking, labels)

Construct a `D`-dimensional sphere of radius `r` with binding sites on its surface.

General constructor, taking the particle's symmetry as the graph encoding `g`.

  - `patch_positions`, `patch_twists`: where each patch sits and how its frame is turned about
    its own normal, see [`normal_pose`](@ref).
  - `colors=1:npatches`: the interaction color of each patch.
  - `vertices=[i:i for i in 1:npatches]`: which vertices of `g` each patch owns. The ranges must
    be disjoint and cover `g`. The default gives each patch one vertex, the layout
    [`cycleencoding`](@ref) produces.
  - `sitesyms=1`: each patch's own rotational symmetry about its normal, a scalar or one entry
    per patch. A patch owning a single vertex leaves no room for a turn about its normal, so its
    symmetry is 1 whatever the arrangement; a patch that needs its own symmetry needs the
    several vertices [`dartencoding`](@ref) gives it.
  - `locking=true`: whether a patch pins its partner's twist, a scalar or one entry per patch.
  - `labels=nothing`: the symmetry label of each patch. `nothing` derives them from the patch
    geometry and colors with [`siteorbits`](@ref); pass one label per patch to impose a labeling
    of your own, which [`check_encoding`](@ref) then has to agree with.
"""
function PatchyParticleSpecies(
    g::NautyDiGraph,
    r::Real,
    patch_positions,
    patch_twists=zeros(float(typeof(r)), length(patch_positions));
    colors=1:length(patch_positions),
    vertices=[i:i for i in eachindex(patch_positions)],
    sitesyms=1,
    locking=true,
    labels=nothing,
)
    D = length(first(patch_positions))
    F = float(eltype(first(patch_positions)))
    r = F(r)
    n = length(patch_positions)

    length(vertices) == n || throw(ArgumentError("expected $n vertex ranges, one per patch, got $(length(vertices))"))
    # Vertices owned by no patch are allowed: an encoding may carry structural ones, as the
    # two-site cycle does to keep a particle from reading as a bond. What is not allowed is two
    # patches claiming one vertex, or a patch claiming one the graph has not got.
    owned = collect(Iterators.flatten(vertices))
    (allunique(owned) && all(v -> 1 <= v <= nv(g), owned)) || throw(
        ArgumentError(
            "the patches' vertex ranges must be disjoint and lie within the graph's $(nv(g)) vertices",
        ),
    )

    tol = sqrt(eps(F)) * r
    poses = [normal_pose(patch_positions[i], patch_twists[i]) for i in 1:n]
    syms = _expandperface(sitesyms, n, "site symmetries")
    locks = _expandperface(locking, n, "locking flags")

    labs = isnothing(labels) ? siteorbits(poses, syms, collect(colors)) : collect(labels)
    length(labs) == n || throw(ArgumentError("expected $n labels, one per patch, got $(length(labs))"))
    stabs = stabilizerorders(poses, syms, labs)

    vertexlabels = zeros(Cint, nv(g))
    for i in 1:n, v in vertices[i]
        vertexlabels[v] = labs[i]
    end
    setlabels!(g, vertexlabels)

    sites = [BindingSite(poses[i], colors[i], vertices[i], tol, tol / r, syms[i], stabs[i], locks[i]) for i in 1:n]
    return check_encoding(PatchyParticleSpecies{D,F,eltype(sites)}(g, sites, r, F.(collect(patch_twists)), tol))
end

"""
    PatchyDisk(angles, r=1; colors=1:length(angles))

A 2D disk of radius `r` with one binding site placed at each angle in `angles` (radians,
measured counterclockwise from the +x axis).
"""
function PatchyDisk(angles, r=1; colors=1:length(angles))
    F = float(eltype(angles))
    r = F(r)
    n = length(angles)
    tol = sqrt(eps(F)) * r
    positions = [SVector(r * cos(F(phi)), r * sin(F(phi))) for phi in angles]
    poses = [normal_pose(positions[i], F(0)) for i in 1:n]

    # A 2D site has no turn about its in-plane normal, so its sitesym is 1 throughout.
    sitesyms = ones(Int, n)
    labels = siteorbits(poses, sitesyms, collect(colors))
    stabs = stabilizerorders(poses, sitesyms, labels)

    g, ranges = cycleencoding(n; labels)
    sites = [BindingSite(poses[i], colors[i], ranges[i], tol, tol / r, 1, stabs[i]) for i in 1:n]
    return check_encoding(PatchyParticleSpecies{2,F,eltype(sites)}(g, sites, r, zeros(F, n), tol))
end

"""
    PatchySphere(p::Polyhedron, r=1; colors=1:nfaces(p), locking=true, twists=0)

A 3D sphere of radius `r` carrying one patch per face of polyhedron `p`, so that the patches inherit the
polyhedron's rotation group.

See [`PolyhedronParticleSpecies`](@ref) for documentation of the keyword arguments.
"""
function PatchySphere(p::Polyhedron, r::Real=1; colors=1:nfaces(p), locking=true, twists=0)
    _patchysphere(p, r, colors, locking, twists, nothing)
end

"""
    _reseat_radially(pose, r)

Move a face's binding site frame onto the sphere of radius `r`: the site sits at `r` along the
centroid direction with a radial normal, keeping its twist reference by re-orthogonalising the
face's local z. The two frames coincide whenever the centroid direction is the face normal, as
on any Platonic solid.
"""
function _reseat_radially(pose::Pose{3,F}, r::Real) where {F}
    ex = normalize(pose.x)
    ez = pose.psi[:, 3]
    ez = normalize(ez - dot(ez, ex) * ex)
    return Pose{3,F,RotMatrix3{F}}(F(r) * ex, RotMatrix3{F}(hcat(ex, cross(ez, ex), ez)))
end

function _patchysphere(p::Polyhedron{F}, r::Real, colors, locking, twists, usecycle::Union{Nothing,Bool}) where {F}
    n = nfaces(p)
    length(colors) == n || throw(ArgumentError("expected $n colors, one per face, got $(length(colors))"))

    r = F(r)
    tol = sqrt(eps(F)) * r
    patchposes(fs) = [_reseat_radially(q, r) for q in _faceposes(p, fs)]

    g, sites = _facesites(
        p,
        patchposes,
        colors,
        _expandperface(locking, n, "locking flags"),
        _expandperface(twists, n, "twists"),
        usecycle,
        tol,
        tol / r,
    )
    # stored as the general constructor takes them: turns of `normal_pose` about each patch normal
    turns = [_twistangle(normal_pose(b.pose.x).psi, b.pose.psi) for b in sites]
    return check_encoding(PatchyParticleSpecies{3,F,eltype(sites)}(g, sites, r, turns, tol))
end

function Base.show(io::Core.IO, ps::PatchyParticleSpecies{D}) where {D}
    return print(io, "$(D)d PatchyParticleSpecies with $(nsites(ps)) sites")
end

Base.copy(ps::PatchyParticleSpecies) = typeof(ps)(copy(ps.g), copy(ps.sites), ps.r, copy(ps.twists), ps.skin)

graphrep(ps::PatchyParticleSpecies) = ps.g
nsites(ps::PatchyParticleSpecies) = length(ps.sites)
bindingsite(ps::PatchyParticleSpecies, i::Integer) = ps.sites[i]
isconvex(::PatchyParticleSpecies) = true

"""
    todict(ps::PatchyParticleSpecies)

Serialize `ps` into a `Dict{String,Any}` of the arguments of its general constructor
[`PatchyParticleSpecies(g, r, patch_positions, patch_twists; ...)`](@ref PatchyParticleSpecies),
whichever constructor built it: its radius, the position, twist, color, site symmetry, locking and
graph vertices of each patch, and the vertex count and edges of its graph.
"""
function todict(ps::PatchyParticleSpecies)
    g = graphrep(ps)
    return Dict{String,Any}(
        "type" => _typetag(ps),
        "r" => ps.r,
        "positions" => [collect(b.pose.x) for b in ps.sites],
        "twists" => copy(ps.twists),
        "colors" => [color(b) for b in ps.sites],
        "sitesyms" => [b.sitesym for b in ps.sites],
        "locking" => [b.locking for b in ps.sites],
        "vertices" => [[first(b.vertices), last(b.vertices)] for b in ps.sites],
        "graph" => Dict{String,Any}("nvertices" => nv(g), "edges" => [[src(e), dst(e)] for e in edges(g)]),
    )
end

"""
    fromdict(::Type{<:PatchyParticleSpecies}, d)

Return the patchy species that `d`, written by [`todict`](@ref), describes, built with the general
constructor.
"""
function fromdict(::Type{<:PatchyParticleSpecies}, d)
    graph = d["graph"]
    g = NautyDiGraph(graph["nvertices"])
    for (u, v) in graph["edges"]
        add_edge!(g, u, v)
    end
    positions = d["positions"]
    D = length(first(positions))
    return PatchyParticleSpecies(
        g,
        Float64(d["r"]),
        [SVector{D,Float64}(x) for x in positions],
        Float64.(d["twists"]);
        colors=Int.(d["colors"]),
        vertices=[a:b for (a, b) in d["vertices"]],
        sitesyms=Int.(d["sitesyms"]),
        locking=Bool.(d["locking"]),
    )
end

function could_contact(::SpeciesAndPose{<:PatchyParticleSpecies}, ::SpeciesAndPose{<:PatchyParticleSpecies}; kwargs...)
    return true # this check would be identical with overlap, so no need to do it twice
end

function overlap(p1::SpeciesAndPose{<:PatchyParticleSpecies}, p2::SpeciesAndPose{<:PatchyParticleSpecies}; kwargs...)
    spcs1, pose1 = p1
    spcs2, pose2 = p2
    return norm(pose1.x - pose2.x) < (spcs1.r + spcs2.r) - (spcs1.skin + spcs2.skin)
end

bounding_radius(ps::PatchyParticleSpecies) = ps.r