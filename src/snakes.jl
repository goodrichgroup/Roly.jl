"""
    SnakeGraph(rules::BindingRules)

A graph representation of `rules`, used to determine whether the number of polyforms compatible with
`rules` is finite.

A `SnakeGraph` encodes the possible snakes, one-dimensional chain-like polyforms, as a directed graph on
the binding sites of the species in `rules`. Each snake corresponds to a walk on the `SnakeGraph`, making
it possible to systematically evaluate them.

Every snake is a succession of bonded binding sites. The `SnakeGraph` has one "outgoing" node for every
non-inert binding site, and one "incoming" node for every such site in every twist it can be bonded in.
Outgoing nodes are connected to the incoming nodes of the sites they can bind to. Incoming nodes are
connected to the outgoing nodes on the same particle. Every edge is associated with a rigid
transformation `h ∈ SE(d)`, making the `SnakeGraph` a graph with group-valued edge labels.

See also [`isacyclic`](@ref), [`isconfined`](@ref), [`isbounded`](@ref).
"""
struct SnakeGraph{R<:BindingRules,P<:Pose} <: AbstractGraph{Int}
    rules::R
    graph::SimpleDiGraph{Int}
    species::Vector{Int}
    sites::Vector{Int}
    twists::Vector{BondTwist}
    incoming::BitVector
    motions::Dict{Edge{Int},P}
end

"""
    SnakeGraph(rules::BindingRules{D}) where {D}

Return the `SnakeGraph` associated with `rules`.
"""
function SnakeGraph(rules::BindingRules{D}) where {D}
    P = posetype(rules)
    R = rotationtype(rules)
    F = numtype(rules)
    spcs, sites, twists, incoming = Int[], Int[], BondTwist[], Bool[]
    nodeindex = Dict{Tuple{Int,Int,BondTwist,Bool},Int}()
    notwist = BondTwist(0, 1)   # an outgoing node is not yet in any bond
    function node!(spc, site, tw, inc)
        return get!(nodeindex, (spc, site, tw, inc)) do
            push!(spcs, spc)
            push!(sites, site)
            push!(twists, tw)
            push!(incoming, inc)
            return length(spcs)
        end
    end
    motions = Dict{Edge{Int},P}()

    # bond edges
    for spc in 1:nspecies(rules)
        for s in 1:nsites(species(rules, spc))
            site = bindingsite(species(rules, spc), s)
            isinert(rules, color(site)) && continue

            exitnode = node!(spc, s, notwist, false)
            for loc in compatible_speciessitelocs(rules, color(site))
                site2 = bindingsite(rules, loc)
                for tw in distinctbondtwists(site, site2)
                    entrynode = node!(loc.species, loc.site, tw, true)
                    motions[Edge(exitnode, entrynode)] = P(zero(SVector{D,F}), R(standard_rotation(F, Val(D), tw)))
                end
            end
        end
    end

    # internal edges
    for i in eachindex(spcs, sites, twists, incoming)
        incoming[i] || continue
        ps = species(rules, spcs[i])
        entrypose = bindingsite(ps, sites[i]).pose
        for k in 1:nsites(ps)
            k == sites[i] && continue
            exitnode = get(nodeindex, (spcs[i], k, notwist, false), nothing)
            isnothing(exitnode) && continue
            motions[Edge(i, exitnode)] = inv(entrypose) * bindingsite(ps, k).pose
        end
    end

    graph = SimpleDiGraph(length(spcs))
    for e in keys(motions)
        add_edge!(graph, e)
    end
    return SnakeGraph{typeof(rules),P}(rules, graph, spcs, sites, twists, BitVector(incoming), motions)
end

bindingrules(sg::SnakeGraph) = sg.rules

Graphs.nv(sg::SnakeGraph) = nv(sg.graph)
Graphs.ne(sg::SnakeGraph) = ne(sg.graph)
Graphs.vertices(sg::SnakeGraph) = vertices(sg.graph)
Graphs.edges(sg::SnakeGraph) = edges(sg.graph)
Graphs.has_vertex(sg::SnakeGraph, v::Integer) = has_vertex(sg.graph, v)
Graphs.has_edge(sg::SnakeGraph, u::Integer, v::Integer) = has_edge(sg.graph, u, v)
Graphs.inneighbors(sg::SnakeGraph, v::Integer) = inneighbors(sg.graph, v)
Graphs.outneighbors(sg::SnakeGraph, v::Integer) = outneighbors(sg.graph, v)
Graphs.is_directed(::Type{<:SnakeGraph}) = true
Graphs.edgetype(::SnakeGraph) = Edge{Int}
Base.eltype(::SnakeGraph) = Int

function Base.show(io::Core.IO, sg::SnakeGraph)
    nin = count(sg.incoming)
    return print(io, "SnakeGraph[", nv(sg) - nin, " out-nodes, ", nin, " in-nodes, ", ne(sg), " edges]")
end

"""
    isacyclic(rules::BindingRules)
    isacyclic(sg::SnakeGraph)

Return `true` if `rules` only allows finite, treelike polyforms, and `false` if its polyforms may contain
cycles or grow without bound.

`true` means that the [`SnakeGraph`](@ref) of `rules` has no cycle, so no snake passes through the same
binding site twice. The converse does not hold: a cycle can be blocked by overlapping particles, so
`false` does not imply that rings or unbounded growth actually occur.

See also [`isconfined`](@ref), [`isbounded`](@ref).
"""
isacyclic(rules::BindingRules) = isacyclic(SnakeGraph(rules))
isacyclic(sg::SnakeGraph) = !is_cyclic(sg.graph)

# Generate the holonomy groups of all strongly connected components.
# We basically construct the directed cycle basis from an arbitrary vertex,
# and then form the group generators correponding to each basis element.
function _holonomies(sg::SnakeGraph)
    P = posetype(sg.rules)
    components = Tuple{Vector{Int},Vector{P}}[]
    for comp in strongly_connected_components(sg.graph)
        length(comp) == 1 && !has_edge(sg.graph, comp[1], comp[1]) && continue
        member = Set(comp)
        frames = Dict(first(comp) => one(P))
        generators = P[]
        queue = [first(comp)]

        # Implicitly construct a spanning tree, once an edge that closes a loop is found,
        # record it as a group generator
        while !isempty(queue)
            u = popfirst!(queue)
            for v in outneighbors(sg.graph, u)
                v in member || continue
                if haskey(frames, v) # found an edge that closes a loop
                    push!(generators, frames[u] * sg.motions[Edge(u, v)] * inv(frames[v]))
                else
                    frames[v] = frames[u] * sg.motions[Edge(u, v)]
                    push!(queue, v)
                end
            end
        end
        push!(components, (comp, generators))
    end
    return components
end

# Whether the rigid motions `gs` share a fixed point, which they do exactly when the group they
# generate keeps every orbit bounded. The motion `x -> R x + t` fixes `x` exactly when `(R - I) x = -t`,
# and the stacked system has a solution exactly when appending `b` leaves the rank of `A` unchanged.
function _hasfixedpoint(gs::AbstractVector{<:Pose{D,F}}, scale::Real) where {D,F}
    A = reduce(vcat, [SMatrix{D,D,F}(g.psi) - I for g in gs])
    b = reduce(vcat, [-g.x for g in gs]) / scale
    return rank([A b]; atol=sqrt(eps(F))) == rank(A; atol=sqrt(eps(F)))
end

# The strongly connected components whose holonomy group has no fixed point: their closed walks can
# carry a snake arbitrarily far. Every other component keeps the particles of a snake inside a bounded
# region, which only finitely many of them fit into.
function _unconfinedcomponents(sg::SnakeGraph)
    scale = maximum(bounding_radius(ps) for ps in species(sg.rules))
    return [comp for (comp, gens) in _holonomies(sg) if !_hasfixedpoint(gens, scale)]
end

# The nodes from which a snake can still grow without bound: those that reach an unconfined component.
function _unconfinednodes(sg::SnakeGraph, components=_unconfinedcomponents(sg))
    unconfined = falses(nv(sg))
    for comp in components
        unconfined[comp] .= true
    end
    queue = findall(unconfined)
    while !isempty(queue)
        v = pop!(queue)
        for u in inneighbors(sg.graph, v)
            unconfined[u] && continue
            unconfined[u] = true
            push!(queue, u)
        end
    end
    return unconfined
end

"""
    isconfined(rules::BindingRules)
    isconfined(sg::SnakeGraph)

Return `true` if the geometry of the binding sites bounds the size of the polyforms of `rules`,
which implies that the number of polyforms is finite. Return `false` if the way particles bind does not by itself
rule out unbounded growth.

`true` means that in every strongly connected component of [`SnakeGraph`](@ref) of `rules`, the motions along the closed
walks share a fixed point. Every snake then stays within a bounded distance of those points, where only finitely many 
particles fit. The converse does not hold: overlapping particles can stop snakes that the motions alone would let 
escape, so `false` does not imply that polyforms grow without bound.

See also [`isacyclic`](@ref), [`isbounded`](@ref).
"""
isconfined(rules::BindingRules) = isconfined(SnakeGraph(rules))
isconfined(sg::SnakeGraph) = isempty(_unconfinedcomponents(sg))


mutable struct _Snake{Q<:Particle}
    particles::Vector{Q}
    nodes::Vector{Int}
end
Base.copy(s::_Snake) = _Snake(copy(s.particles), copy(s.nodes))
function Base.copy!(dst::_Snake, src::_Snake)
    copy!(dst.particles, src.particles)
    copy!(dst.nodes, src.nodes)
    return dst
end

# Writes the `j`-th extension of `v` into `u`, forming every contact the new particle makes, so a
# snake may cross-link onto itself. Returns `missing` if it does not fit and `nothing` past the last.
# The extensions of a snake depend only on the node its newest particle was entered through, and
# `moves` lists them for every node, as the bond edges `(x, y)` that lead on from it.
function _extendsnake!(u::_Snake, v::_Snake, j::Integer, sg::SnakeGraph, moves::Vector{Vector{Tuple{Int,Int}}})
    nodemoves = moves[v.nodes[end]]
    j > length(nodemoves) && return nothing
    rules = sg.rules
    x, y = nodemoves[j]
    part = v.particles[end]
    ps = species(rules, part.speciesindex)
    mate = bindingsite(species(rules, sg.species[y]), sg.sites[y])
    pose = part.pose * bindingsite(ps, sg.sites[x]).pose * sg.motions[Edge(x, y)] * inv(mate.pose)
    new = typeof(part)(pose, leadingvertex(part) + nv(graphrep(ps)), sg.species[y])
    first(_overlap_and_contacts(v.particles, new, rules)) && return missing
    copy!(u, v)
    push!(u.particles, new)
    push!(u.nodes, y)
    return u
end

# Searches the snakes of up to `maxlength` particles that start from a particle entered through the
# in-node `root`, and returns the first periodic witness found together with the search status. A new
# particle is compared with every earlier one entered through the same node, the latest first.
function _walksnakes(sg::SnakeGraph, unconfined::BitVector, root::Integer, maxlength::Integer)
    rules = sg.rules
    # the bond edges leading on from each node towards unconfined nodes; the possible snake moves
    moves = [[(x, y) for x in outneighbors(sg.graph, v) for y in outneighbors(sg.graph, x) if unconfined[y]]
             for v in vertices(sg)]
    start = _Snake([Particle(rules, sg.species[root]; leadingvertex=1)], [root])
    P, Q = posetype(rules), particletype(rules)
    witness = Ref{Union{Nothing,_PeriodicSnake{P,Q}}}(nothing)
    extend!(u, v, j, _) = _extendsnake!(u, v, j, sg, moves)
    result = _backtrack(extend!, start; maxdepth=maxlength - 1) do s, _
        n = length(s.particles)
        for j in (n - 1):-1:1
            s.nodes[j] == s.nodes[n] || continue
            w = _screwescape(s.particles, j, n, rules)
            isnothing(w) && continue
            witness[] = w
            return BREAK
        end
        return ACCEPT
    end
    return witness[], result.status
end

# Searche the snakes of `sg` with iterative deepening, from length 2 and doubling up to `maxlength`,
# so that the shortest periodic witness is found first. Return a `_PeriodicSnake` and status code.
function _periodicsnake(sg::SnakeGraph, unconfined::BitVector=_unconfinednodes(sg); maxlength::Integer)
    roots = Int[]
    seen = Set{Tuple{Int,Int}}()
    for v in vertices(sg)
        (sg.incoming[v] && unconfined[v]) || continue
        key = (sg.species[v], sg.sites[v])
        key in seen && continue
        push!(seen, key)
        push!(roots, v)
    end
    isempty(roots) && return nothing, :bounded

    # generate snakes up to length L, which is doubled every round until it hits maxlength
    L = min(2, maxlength)
    while true
        exhausted = true
        for root in roots
            w, status = _walksnakes(sg, unconfined, root, L)
            isnothing(w) || return w, :unbounded
            status == MaxDepthReached && (exhausted = false)
        end
        exhausted && return nothing, :bounded
        L == maxlength && return nothing, :undecided
        L = min(2L, maxlength)
    end
end
_periodicsnake(rules::BindingRules; kwargs...) = first(_periodicsnake(SnakeGraph(rules); kwargs...))

struct _PeriodicSnake{P<:Pose,Q<:Particle}
    motion::P
    particles::Vector{Q}
end

# Return the range `i:(j - 1)` of `particles` and a repeatable rigid motion (as a `_PeriodicSnake`).
# If no repeatable motion exists, return nothing.
function _screwescape(particles::AbstractVector{<:Particle}, i::Integer, j::Integer, rules::BindingRules)
    C = particles[i:(j - 1)]
    g = particles[j].pose * inv(particles[i].pose)

    axis, h = _screwclimb(g)
    isnothing(h) && return nothing               # a pure rotation cannot escape its own annulus

    us = (dot(p.pose.x, axis) for p in C)
    rmax = maximum(bounding_radius(species(rules, speciesindex(p))) for p in C)
    extent = maximum(us) - minimum(us) + 2 * rmax
    gm = g
    for _ in 1:ceil(Int, extent / abs(h))
        for p in C
            moved = gm * p
            _overlap_and_contacts(C, moved, rules)[1] && return nothing
        end
        gm = gm * g
    end
    return _PeriodicSnake(g, C)
end

# Return `(axis, climb)` of a rigid motion. Any rigid motion can be viewed as screw along `axis`,
# where `climb` is the distance along the axis the body is translated by. (The normal distance from the axis
# is not returned by this function.)
function _screwclimb(g::Pose{D,F}) where {D,F}
    tol = sqrt(eps(F))
    turning = abs(rem(rotation_angle(g.psi), 2 * F(π), RoundNearest)) > tol
    if !turning
        # with no turn there is no rotation axis, so the axis is just the direction of translation
        n = norm(g.x)
        n > tol || return nothing, nothing        # the identity: no repeat at all
        return g.x / n, n
    end
    D == 3 || return nothing, nothing
    axis = rotation_axis(g.psi)
    climb = dot(g.x, axis)
    return axis, abs(climb) > tol ? climb : nothing
end

"""
    isbounded(rules::BindingRules[; maxlength=4nv(SnakeGraph(rules)])

Return `true` if `rules` admits finitely many polyforms and `false` otherwise. Return `missing` if the problem
cannot be decided within the imposed `maxlength`.

Whether a set of binding rules is bounded is mathematically undecidable in general [1, 2]. This function applies
increasingly strict tests that nevertheless often give a definite answer in practice.

Every infinite polyform contains an infinite snake, a 1d chain of particles, so boundedness can be decided by 
considering only the snakes of `rules`. This is done in three steps:

 1. if the [`SnakeGraph`](@ref) of `rules` has no cycle, every snake if finite and the rules must be bounded.
 2. if, in every strongly connected component of that graph, the motions along the closed walks share a fixed point,
    every snake stays inside a bounded region. Particle overlaps then guarantee that the rules must be bounded.
 3. otherwise, snakes are enumerated up to `maxlength`. If the enumeration finds a snake that can indefinitely
    bind to itself, this proves that the rules are unbounded. Otherwise, if we can exhaustively enumerate all snakes
    the rules must be bounded. If the length-limit starts binding, the problem cannot be decided and a larger 
    `maxlength` is needed, but there is no guarantee that a finite `maxlength` will lead to a settled outcome.

See also [`isacyclic`](@ref), [`isconfined`](@ref).

# References

 1. Kari, J. Infinite snake tiling problems. Developments in Language Theory (eds Ito, M. & Toyama, M.),
    Lecture Notes in Computer Science vol. 2450 (Springer, 2003). https://doi.org/10.1007/3-540-45005-X_6
 2. Adleman, L., Kari, J., Kari, L. & Reishus, D. On the decidability of self-assembly of infinite ribbons.
    Proc. 43rd Annual IEEE Symposium on Foundations of Computer Science 1-8 (2002). https://doi.org/10.1109/sfcs.2002.1181977
"""
function isbounded(rules::BindingRules; maxlength::Union{Nothing,Integer}=nothing)
    sg = SnakeGraph(rules)
    isacyclic(sg) && return true

    unconfined = _unconfinednodes(sg)
    any(unconfined) || return true

    L = something(maxlength, 4 * nv(sg))
    _, verdict = _periodicsnake(sg, unconfined; maxlength=L)

    verdict == :unbounded && return false
    verdict == :bounded && return true
    return missing
end
