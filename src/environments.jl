"""
    ParticleEnvironment

The environment of a particle: every particle within a fixed graph distance of it, as an
isomorphism class of the bond graph with the root particle marked.

  - `graph`: canonized graph representation, with the labels of the root's vertices offset so the
    root is distinguished
  - `rootvertex`: canonical index of one of the root's site vertices
  - `rules`: the `BindingRules` the labels refer to

Identity (`hash`/`==`) is the graph; the radius is a property of the extraction or enumeration
that produced the environment. Every environment is itself a valid polyform.
"""
struct ParticleEnvironment{G<:AbstractNautyGraph,S<:BindingRules}
    graph::G
    rootvertex::Int
    rules::S
end

Base.hash(e::ParticleEnvironment, h::UInt) = hash(e.graph, h)
Base.:(==)(a::ParticleEnvironment, b::ParticleEnvironment) = a.rules === b.rules && a.graph == b.graph
Base.show(io::Core.IO, e::ParticleEnvironment) = print(io, "ParticleEnvironment[nv=$(nv(e.graph))]")

# Root marks are added in multiples of this, so marked labels never collide with site labels.
function _markoffset(rules::BindingRules)
    m = 0
    for i in 1:nspecies(rules)
        m = max(m, maximum(labels(graphrep(species(rules, i)))))
    end
    return m
end

# Mark and canonize: bump the labels of the i-th group of vertices by i*offset. Marking every vertex
# of a root particle (rather than a single one) preserves the particle's internal symmetry, so e.g.
# a threefold-symmetric root does not split one class into three.
#
# Site-resolution marks (unbound vs outside-the-crop) need no encoding: they are determined by graph
# distance from the roots, which any isomorphism of marked graphs preserves.
function _canonmarked(g::AbstractNautyGraph, groups, offset)
    h = copy(g)
    for (i, vs) in enumerate(groups), v in vs
        setlabel!(h, v, labels(h)[v] + i * offset)
    end
    # `canonize!` returns the new -> old permutation (`newlabels[i] == oldlabels[perm[i]]`); the
    # callers hold old indices and want new ones, hence the inversion.
    return h, invperm(collect(Int, canonize!(h)))
end

# Buffers for the particle-hop search, reusable across calls and polyforms.
struct EnvironmentBuffers
    dist::Vector{Int}               # particle -> bond-hops from the roots, -1 if unreached
    queue::Vector{Int}
    vertex2particle::Vector{Int}    # original vertex -> owning particle
    neighbors::Vector{Int}          # neighbor scratch for `bfs!`
end
EnvironmentBuffers() = EnvironmentBuffers(Int[], Int[], Int[], Int[])
function Base.copy(bufs::EnvironmentBuffers)
    return EnvironmentBuffers(copy(bufs.dist), copy(bufs.queue), copy(bufs.vertex2particle),
                              copy(bufs.neighbors))
end

# Particle-hop distances from the `roots` via `bfs!`: vertices of the same particle are free, bonds
# cost one hop. Capped at `maxdepth`, so the cost is the ball, not the whole polyform. Fills and
# returns `bufs.dist`.
function _particledists!(bufs::EnvironmentBuffers, poly::Polyform, roots; maxdepth)
    rules = bindingrules(poly)
    g = graphrep(poly)

    resize!(bufs.dist, nparticles(poly))
    resize!(bufs.queue, nparticles(poly))
    resize!(bufs.vertex2particle, nv(g))
    for (i, part) in enumerate(poly.particles), ov in graphvertices(part, rules)
        bufs.vertex2particle[ov] = i
    end

    fill!(bufs.dist, -1)
    bfs!(bufs.dist, bufs.queue, roots; maxdepth) do p
        # Any edge leaving the particle is a bond; interior edges stay within it.
        empty!(bufs.neighbors)
        for ov in graphvertices(particles(poly, p), rules)
            outneighs = NautyGraphs.adjrow(g, tocanon(poly, ov))
            for w in eachindex(outneighs)
                outneighs[w] || continue
                q = bufs.vertex2particle[toorig(poly, w)]
                q == p || push!(bufs.neighbors, q)
            end
        end
        bufs.neighbors
    end
    return bufs.dist
end

"""
    ParticleEnvironment(poly::Polyform, root::Integer; depth::Integer)

Extract the environment of particle `root` in `poly`: the ball of all particles within graph
distance `depth`.
"""
function ParticleEnvironment(poly::Polyform, root::Integer; depth::Integer, bufs=EnvironmentBuffers())
    rules = bindingrules(poly)
    dist = _particledists!(bufs, poly, (root,); maxdepth=depth)

    verts = Int[]
    offsets = zeros(Int, nparticles(poly))   # particle -> position of its first vertex in `verts`
    for p in 1:nparticles(poly)
        0 <= dist[p] <= depth || continue
        offsets[p] = length(verts) + 1
        for ov in graphvertices(particles(poly, p), rules)
            push!(verts, tocanon(poly, ov))
        end
    end

    h = graphrep(poly)[verts]
    rootgroup = offsets[root]:(offsets[root] + nsites(particles(poly, root), rules) - 1)
    hm, old2new = _canonmarked(h, (rootgroup,), _markoffset(rules))
    return ParticleEnvironment(hm, old2new[offsets[root]], rules)
end

"""
    particleenvironments(poly::Polyform; depth::Integer)

Return the environment of every particle of `poly`, in particle order.
"""
function particleenvironments(poly::Polyform; depth::Integer)
    bufs = EnvironmentBuffers()
    return [ParticleEnvironment(poly, p; depth, bufs) for p in 1:nparticles(poly)]
end

# ------------------------------------------------------------------------------------------------
# Enumeration of all particle environments of a rule set, by reverse search over rooted balls.

mutable struct EnvironmentState{P<:Polyform,G<:AbstractNautyGraph}
    poly::P               # working polyform; the root is always particle 1
    key::G                # root-marked canonical graph, the reverse-search identity
    rootvertex::Int       # canonical index of the root's first site vertex in `key`
end

function EnvironmentState(rules::BindingRules)
    poly = Polyform(rules)
    return EnvironmentState(poly, copy(graphrep(poly)), 0)
end
Base.copy(s::EnvironmentState) = EnvironmentState(copy(s.poly), copy(s.key), s.rootvertex)
function Base.copy!(dst::EnvironmentState, src::EnvironmentState)
    copy!(dst.poly, src.poly)
    copy!(dst.key, src.key)
    dst.rootvertex = src.rootvertex
    return dst
end

_rootgroup(poly::Polyform) =
    [tocanon(poly, ov) for ov in graphvertices(particles(poly, 1), bindingrules(poly))]

# Recompute the root-marked canonical key of a state; the whole polyform is its own ball.
function _rekey!(s::EnvironmentState)
    poly = s.poly
    hm, old2new = _canonmarked(graphrep(poly), (_rootgroup(poly),), _markoffset(bindingrules(poly)))
    s.key = hm
    s.rootvertex = old2new[tocanon(poly, first(graphvertices(particles(poly, 1), bindingrules(poly))))]
    return s
end

mutable struct EnvironmentEnumAux{BS<:BindingSite,G<:AbstractNautyGraph}
    seen::Set{G}                              # offspring keys of the current parent
    pairs::Vector{Tuple{BS,SpeciesSiteLoc}}   # attachments that stay within the ball
    bufs::EnvironmentBuffers
    depth::Int
end
function Base.copy(aux::EnvironmentEnumAux)
    return EnvironmentEnumAux(copy(aux.seen), copy(aux.pairs), copy(aux.bufs), aux.depth)
end

# Attachment options whose owner particle is at distance <= depth-1 from the root, so the new
# particle lands inside the ball.
function _collectballpairs!(aux::EnvironmentEnumAux, poly::Polyform)
    rules = bindingrules(poly)
    empty!(aux.pairs)
    dist = _particledists!(aux.bufs, poly, (1,); maxdepth=aux.depth)
    for (p, part) in enumerate(poly.particles)
        0 <= dist[p] <= aux.depth - 1 || continue
        for k in 1:nsites(part, rules)
            site = bindingsite(part, rules, k)
            _isbound_vertex(poly, part, first(site.vertices); canonidxs=false) && continue
            isinert(rules, color(site)) && continue
            for siteloc in compatible_speciessitelocs(rules, color(site))
                push!(aux.pairs, (site, siteloc))
            end
        end
    end
    return aux.pairs
end

function _adjenv!(u::EnvironmentState, v::EnvironmentState, j::Integer, aux::EnvironmentEnumAux)
    rules = bindingrules(v.poly)
    if nparticles(v.poly) == 0
        j > nspecies(rules) && return nothing
        copy!(u.poly, Polyform(rules, j))
        return _rekey!(u)
    end

    j == 1 && _collectballpairs!(aux, v.poly)
    j > length(aux.pairs) && return nothing

    site, siteloc = aux.pairs[j]
    copy!(u.poly, v.poly)
    ismissing(raise!(u.poly, site, siteloc)) && return missing
    _rekey!(u)
    u.key in aux.seen && return missing
    push!(aux.seen, copy(u.key))
    return u
end

# The parent drops the canonically-last particle among those at maximal distance from the root.
# Removing a maximal-distance particle keeps a ball a ball with every other distance unchanged: a
# shortest root-path passes only through vertices of strictly smaller distance, so it never routes
# through a maximal-distance one.
function _lsenv!(w::EnvironmentState, v::EnvironmentState, bufs::EnvironmentBuffers)
    n = nparticles(v.poly)
    rules = bindingrules(v.poly)
    if n <= 1
        copy!(w.poly, Polyform(rules))
        copy!(w.key, graphrep(w.poly))
        w.rootvertex = 0
        return w
    end

    dist = _particledists!(bufs, v.poly, (1,); maxdepth=typemax(Int))
    D = maximum(dist)
    _, old2new = _canonmarked(graphrep(v.poly), (_rootgroup(v.poly),), _markoffset(rules))
    canonpos(p) = minimum(old2new[tocanon(v.poly, ov)]
                          for ov in graphvertices(particles(v.poly, p), rules))
    drop = argmax(canonpos, (p for p in 1:n if dist[p] == D))

    w.poly = subpolyform(v.poly, [p for p in 1:n if p != drop])
    return _rekey!(w)
end

"""
    particleenvironments(f, rules::BindingRules; depth, maxsize=Inf, maxstrs=Inf)

Enumerate every particle environment of `rules` at radius `depth`, streaming each to
`f(env, nparticles)`. Each class is visited exactly once; `f` may return `REJECT` to prune all
extensions of an environment, or `BREAK` to stop.
"""
function particleenvironments(f, rules::BindingRules; depth::Integer, maxsize=Inf, maxstrs=Inf,
                              kwargs...)
    v₀ = EnvironmentState(rules)
    BS = BindingSite{posetype(rules),numtype(rules)}
    G = typeof(graphrep(v₀.poly))
    aux = EnvironmentEnumAux(Set{G}(), Tuple{BS,SpeciesSiteLoc}[], EnvironmentBuffers(), Int(depth))
    lsbufs = EnvironmentBuffers()
    rsys = RSSystem((w, v) -> _lsenv!(w, v, lsbufs), _adjenv!, v₀;
                    compare=(a, b) -> a.key == b.key, aux)

    frs = (s, _) -> f(ParticleEnvironment(copy(s.key), s.rootvertex, rules), nparticles(s.poly))
    return reversesearch(frs, rsys; maxdepth=maxsize, maxverts=maxstrs + 1, kwargs...)
end

"""
    particleenvironments(rules::BindingRules; depth, kwargs...)

Return all particle environments of `rules` at radius `depth`.
"""
function particleenvironments(rules::BindingRules; depth::Integer, kwargs...)
    G = typeof(graphrep(Polyform(rules)))
    envs = ParticleEnvironment{G,typeof(rules)}[]
    particleenvironments((e, _) -> (push!(envs, e); ACCEPT), rules; depth, kwargs...)
    return envs
end

# ------------------------------------------------------------------------------------------------
# Graph surgery on an environment: particles are the components of the interior edges, bonds are
# the bidirectional edges between components, and crops are component-level breadth-first
# searches. No poses are involved.

# Partition the vertices of an environment graph into particles: bond edges are bidirectional,
# interior edges are not. (Same assumption as `exterior_edges`: no double edges inside one
# particle's graph.)
function _components(g::AbstractNautyGraph)
    nvg = nv(g)
    comp = zeros(Int, nvg)
    ncomp = 0
    stack = Int[]
    for v₀ in 1:nvg
        comp[v₀] == 0 || continue
        ncomp += 1
        comp[v₀] = ncomp
        push!(stack, v₀)
        while !isempty(stack)
            v = pop!(stack)
            for w in outneighbors(g, v)
                (has_edge(g, w, v) || comp[w] != 0) && continue
                comp[w] = ncomp
                push!(stack, w)
            end
            for w in inneighbors(g, v)
                (has_edge(g, v, w) || comp[w] != 0) && continue
                comp[w] = ncomp
                push!(stack, w)
            end
        end
    end
    compverts = [Int[] for _ in 1:ncomp]
    for v in 1:nvg
        push!(compverts[comp[v]], v)
    end
    return comp, compverts
end

# Particle-level adjacency of an environment graph; duplicate entries are harmless to `bfs!`.
function _componentadjacency(g, comp, ncomp)
    adj = [Int[] for _ in 1:ncomp]
    for v in vertices(g), w in outneighbors(g, v)
        (has_edge(g, w, v) && comp[v] != comp[w]) || continue
        push!(adj[comp[v]], comp[w])
    end
    return adj
end

# Vertices of the components reached by a crop search, with their positions in the kept order.
function _keptvertices(compverts, dist, nvg)
    verts = Int[]
    pos = zeros(Int, nvg)
    for c in eachindex(compverts)
        dist[c] >= 0 || continue
        for v in compverts[c]
            push!(verts, v)
            pos[v] = length(verts)
        end
    end
    return verts, pos
end

# The induced subgraph on `verts` with every old root mark stripped and the components `first`
# and `second` marked once and twice, canonized.
function _markedcrop(g, verts, comp, offset, first, second)
    h = g[verts]
    for (k, v) in enumerate(verts)
        l = mod1(labels(h)[k], offset)
        setlabel!(h, k, l + (comp[v] == first ? offset : comp[v] == second ? 2offset : 0))
    end
    canonize!(h)
    return h
end

"""
    rootbonds(env::ParticleEnvironment, depth)

Return one entry per bond of the root of `env`, as a named tuple `(view, reversed, sites)`.

  - `view`: the canonical graph of every particle within `depth - 1` of either end of the bond,
    with the root marked first and its partner second
  - `reversed`: the same crop with the two marks swapped
  - `sites`: the labels of the two bonded sites, root first

Two root bonds, of the same or of different environments, can be the two ends of one bond in a
structure exactly when the `view` of one equals the `reversed` of the other. The crop within the
ball equals the crop within any structure containing it, since every path of length at most
`depth - 1` from a root neighbor stays inside the ball.
"""
function rootbonds(env::ParticleEnvironment, depth::Integer)
    depth >= 1 || throw(ArgumentError("root bonds need `depth >= 1`"))
    g = env.graph
    offset = _markoffset(env.rules)
    comp, compverts = _components(g)
    adj = _componentadjacency(g, comp, length(compverts))
    rootcomp = comp[env.rootvertex]

    dist = fill(-1, length(compverts))
    queue = zeros(Int, length(compverts))
    out = NamedTuple{(:view, :reversed, :sites),Tuple{typeof(g),typeof(g),Tuple{Int,Int}}}[]
    for u in compverts[rootcomp], w in outneighbors(g, u)
        (has_edge(g, w, u) && comp[w] != rootcomp) || continue
        partner = comp[w]
        fill!(dist, -1)
        bfs!(c -> adj[c], dist, queue, (rootcomp, partner); maxdepth=depth - 1)
        verts, _ = _keptvertices(compverts, dist, nv(g))
        push!(out, (view=_markedcrop(g, verts, comp, offset, rootcomp, partner),
                    reversed=_markedcrop(g, verts, comp, offset, partner, rootcomp),
                    sites=(mod1(labels(g)[u], offset), mod1(labels(g)[w], offset))))
    end
    return out
end

"""
    crop(env::ParticleEnvironment, depth)

Return the environment of the same root at radius `depth`: the sub-ball of all particles within
`depth` of the root.

Distances inside a ball equal the ambient distances to the root, so the crop of the ball is the
ball the ambient structure would have produced at the smaller radius.
"""
function crop(env::ParticleEnvironment, depth::Integer)
    depth >= 0 || throw(ArgumentError("`depth` must be nonnegative"))
    g = env.graph
    comp, compverts = _components(g)
    adj = _componentadjacency(g, comp, length(compverts))
    rootcomp = comp[env.rootvertex]

    dist = fill(-1, length(compverts))
    queue = zeros(Int, length(compverts))
    bfs!(c -> adj[c], dist, queue, (rootcomp,); maxdepth=depth)

    verts, pos = _keptvertices(compverts, dist, nv(g))
    # the root's mark survives the crop unchanged, so no relabeling is needed
    h = g[verts]
    old2new = invperm(collect(Int, canonize!(h)))
    return ParticleEnvironment(h, old2new[pos[env.rootvertex]], env.rules)
end
