"""
    Environment{N}

The local environment of `N` root particles: all particles within graph distance `depth` of a
root, as an isomorphism class. See [`ParticleEnvironment`](@ref) and [`BondEnvironment`](@ref)
for the two instantiations in use.

  - `graph`: canonized graph representation with the root particles distinguished, in order
  - `rootvertices`: one canonical site vertex per root particle
  - `depth`: crop radius around each root
  - `rules`: the `BindingRules` the labels refer to

Identity (`hash`/`==`) is `(graph, depth)`; `rootvertices` is convenience metadata. Every valid
environment is itself a valid polyform.
"""
struct Environment{N,G<:AbstractNautyGraph,S<:BindingRules}
    graph::G
    rootvertices::NTuple{N,Int}
    depth::Int
    rules::S
end

"""
    ParticleEnvironment

[`Environment`](@ref) of a single particle: the ball of all particles within graph distance
`depth` of the root. `rootvertices` holds the root's first site vertex.
"""
const ParticleEnvironment{G,S} = Environment{1,G,S}

"""
    BondEnvironment

[`Environment`](@ref) of a bond: all particles within graph distance `depth` of either endpoint,
with the endpoints as roots, in order. `rootvertices` holds the anchoring bond's two site
vertices — metadata only, since parallel bonds between the same particle pair share a class.
`reverse` swaps the two roots; a bond environment equal to its own reversal is symmetric.
"""
const BondEnvironment{G,S} = Environment{2,G,S}

Base.hash(e::Environment, h::UInt) = hash(e.graph, hash(e.depth, h))
function Base.:(==)(a::Environment{N}, b::Environment{N}) where {N}
    return a.depth == b.depth && a.rules === b.rules && a.graph == b.graph
end
Base.:(==)(::Environment, ::Environment) = false

_envname(N) = N == 1 ? "ParticleEnvironment" : N == 2 ? "BondEnvironment" : "Environment{$N}"
Base.show(io::Core.IO, e::Environment{N}) where {N} = print(io, _envname(N), "[k=$(e.depth), nv=$(nv(e.graph))]")

# Root marks are added in multiples of this, so marked labels never collide with site labels.
function _markoffset(sys::BindingRules)
    m = 0
    for i in 1:nspecies(sys)
        m = max(m, maximum(labels(graphrep(species(sys, i)))))
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
    sys = bindingrules(poly)
    g = graphrep(poly)

    resize!(bufs.dist, nparticles(poly))
    resize!(bufs.queue, nparticles(poly))
    resize!(bufs.vertex2particle, nv(g))
    for (i, part) in enumerate(poly.particles), ov in graphvertices(part, sys)
        bufs.vertex2particle[ov] = i
    end

    fill!(bufs.dist, -1)
    bfs!(bufs.dist, bufs.queue, roots; maxdepth) do p
        # Any edge leaving the particle is a bond; interior edges stay within it.
        empty!(bufs.neighbors)
        for ov in graphvertices(particles(poly, p), sys)
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

# Induced marked subgraph on the particles with dist in [0, depth], with the `roots` particles
# marked in order. `rootvertices` gives one original site vertex per root; the returned tuple holds
# their canonical indices in the new graph.
function _envgraph(poly::Polyform, dist::Vector{Int}, depth::Integer, roots, rootvertices)
    sys = bindingrules(poly)
    verts = Int[]
    offsets = zeros(Int, nparticles(poly))   # particle -> position of its first vertex in `verts`
    for p in 1:nparticles(poly)
        0 <= dist[p] <= depth || continue
        offsets[p] = length(verts) + 1
        for ov in graphvertices(particles(poly, p), sys)
            push!(verts, tocanon(poly, ov))
        end
    end

    h = graphrep(poly)[verts]
    groups = map(p -> offsets[p]:(offsets[p] + nsites(particles(poly, p), sys) - 1), roots)
    hm, old2new = _canonmarked(h, groups, _markoffset(sys))

    canonrootvertices = map(roots, rootvertices) do p, ov
        old2new[offsets[p] + ov - leading_vertex(particles(poly, p))]
    end
    return hm, canonrootvertices
end

"""
    ParticleEnvironment(poly::Polyform, root::Integer; depth::Integer)

Extract the environment of particle `root` in `poly`: the ball of all particles within graph
distance `depth`.
"""
function ParticleEnvironment(poly::Polyform, root::Integer; depth::Integer, bufs=EnvironmentBuffers())
    dist = _particledists!(bufs, poly, (root,); maxdepth=depth)
    rootvertex = first(graphvertices(particles(poly, root), bindingrules(poly)))
    hm, rootvertices = _envgraph(poly, dist, depth, (root,), (rootvertex,))
    return Environment(hm, rootvertices, Int(depth), bindingrules(poly))
end

"""
    BondEnvironment(poly::Polyform, bond; depth::Integer)

Extract the environment of `bond` (an element of `bonds(poly)`): all particles within graph
distance `depth` of either endpoint. The endpoints become the environment's roots, in the order
given.
"""
function BondEnvironment(poly::Polyform, bond::Pair; depth::Integer, bufs=EnvironmentBuffers())
    sys = bindingrules(poly)
    (p1, s1), (p2, s2) = bond.first, bond.second
    p1 == p2 && throw(ArgumentError("bond endpoints must be distinct particles"))

    dist = _particledists!(bufs, poly, (p1, p2); maxdepth=depth)
    sitevertex(p, s) = first(bindingsites(particles(poly, p), sys, s).vertices)
    hm, rootvertices = _envgraph(poly, dist, depth, (p1, p2),
                                 (sitevertex(p1, s1), sitevertex(p2, s2)))
    return Environment(hm, rootvertices, Int(depth), sys)
end

"""
    particleenvironments(poly::Polyform; depth::Integer)

Return the environment of every particle of `poly`, in particle order.
"""
function particleenvironments(poly::Polyform; depth::Integer)
    bufs = EnvironmentBuffers()
    return [ParticleEnvironment(poly, p; depth, bufs) for p in 1:nparticles(poly)]
end

"""
    reverse(env::BondEnvironment)

Return the environment with the two roots swapped.
"""
function Base.reverse(env::BondEnvironment)
    offset = _markoffset(env.rules)
    h = copy(env.graph)
    for v in vertices(h)
        l = labels(h)[v]
        if l > 2offset
            setlabel!(h, v, l - offset)
        elseif l > offset
            setlabel!(h, v, l + offset)
        end
    end
    old2new = invperm(collect(Int, canonize!(h)))
    return Environment(h, (old2new[env.rootvertices[2]], old2new[env.rootvertices[1]]), env.depth,
                       env.rules)
end

# ------------------------------------------------------------------------------------------------
# Enumeration of all particle environments of a rule set, by reverse search over rooted balls.

mutable struct EnvironmentState{P<:Polyform,G<:AbstractNautyGraph}
    poly::P               # working polyform; the root is always particle 1
    key::G                # root-marked canonical graph, the reverse-search identity
    rootvertex::Int       # canonical index of the root's first site vertex in `key`
end

function EnvironmentState(sys::BindingRules)
    poly = Polyform(sys)
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
    pairs::Vector{Tuple{BS,BindingSiteLoc}}   # attachments that stay within the ball
    bufs::EnvironmentBuffers
    depth::Int
end
function Base.copy(aux::EnvironmentEnumAux)
    return EnvironmentEnumAux(copy(aux.seen), copy(aux.pairs), copy(aux.bufs), aux.depth)
end

# Attachment options whose owner particle is at distance <= depth-1 from the root, so the new
# particle lands inside the ball.
function _collectballpairs!(aux::EnvironmentEnumAux, poly::Polyform)
    sys = bindingrules(poly)
    empty!(aux.pairs)
    dist = _particledists!(aux.bufs, poly, (1,); maxdepth=aux.depth)
    for (p, part) in enumerate(poly.particles)
        0 <= dist[p] <= aux.depth - 1 || continue
        for k in 1:nsites(part, sys)
            site = bindingsites(part, sys, k)
            isbound_vertex(poly, tocanon(poly, first(site.vertices))) && continue
            isinert(sys, color(site)) && continue
            for siteloc in compatible_sitelocs(sys, color(site))
                push!(aux.pairs, (site, siteloc))
            end
        end
    end
    return aux.pairs
end

function _adjenv!(u::EnvironmentState, v::EnvironmentState, j::Integer, aux::EnvironmentEnumAux)
    sys = bindingrules(v.poly)
    if nparticles(v.poly) == 0
        j > nspecies(sys) && return nothing
        copy!(u.poly, Polyform(sys, j))
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
    sys = bindingrules(v.poly)
    if n <= 1
        copy!(w.poly, Polyform(sys))
        copy!(w.key, graphrep(w.poly))
        w.rootvertex = 0
        return w
    end

    dist = _particledists!(bufs, v.poly, (1,); maxdepth=typemax(Int))
    D = maximum(dist)
    _, old2new = _canonmarked(graphrep(v.poly), (_rootgroup(v.poly),), _markoffset(sys))
    canonpos(p) = minimum(old2new[tocanon(v.poly, ov)]
                          for ov in graphvertices(particles(v.poly, p), sys))
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
    aux = EnvironmentEnumAux(Set{G}(), Tuple{BS,BindingSiteLoc}[], EnvironmentBuffers(), Int(depth))
    lsbufs = EnvironmentBuffers()
    rsys = RSSystem((w, v) -> _lsenv!(w, v, lsbufs), _adjenv!, v₀;
                    compare=(a, b) -> a.key == b.key, aux)

    frs = (s, _) -> f(Environment(copy(s.key), (s.rootvertex,), Int(depth), rules), nparticles(s.poly))
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
# Bond environments of a particle environment, by pure graph surgery: particles are the components
# of the interior edges, root-incident bonds are the bidirectional edges leaving the root
# component, and the crop is a component-level breadth-first search. No poses are involved.

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

"""
    bondenvironments(env::ParticleEnvironment)

Return the environments of every root-incident bond of `env`, cropped at `depth - 1`, with the
root as the first endpoint.

The crop inside the ball equals the crop inside any larger structure containing it: every path of
length at most `depth - 1` from a root neighbor stays within the ball.
"""
function bondenvironments(env::ParticleEnvironment)
    env.depth >= 1 || throw(ArgumentError("bond environments need `depth >= 1`"))
    g = env.graph
    offset = _markoffset(env.rules)
    comp, compverts = _components(g)
    ncomp = length(compverts)
    rootcomp = comp[env.rootvertices[1]]

    # particle-level adjacency; duplicate entries are harmless to `bfs!`
    adj = [Int[] for _ in 1:ncomp]
    for v in vertices(g), w in outneighbors(g, v)
        (has_edge(g, w, v) && comp[v] != comp[w]) || continue
        push!(adj[comp[v]], comp[w])
    end

    dist = zeros(Int, ncomp)
    queue = zeros(Int, ncomp)
    pos = zeros(Int, nv(g))
    out = BondEnvironment{typeof(env.graph),typeof(env.rules)}[]
    for u in compverts[rootcomp], w in outneighbors(g, u)
        (has_edge(g, w, u) && comp[w] != rootcomp) || continue
        partner = comp[w]

        # crop: every particle within depth-1 of either endpoint
        fill!(dist, -1)
        bfs!(c -> adj[c], dist, queue, (rootcomp, partner); maxdepth=env.depth - 1)

        verts = Int[]
        for c in 1:ncomp
            dist[c] >= 0 || continue
            for v in compverts[c]
                push!(verts, v)
                pos[v] = length(verts)
            end
        end

        h = g[verts]
        for (k, v) in enumerate(verts)
            l = mod1(labels(h)[k], offset)   # strip the old root mark
            l += comp[v] == rootcomp ? offset : comp[v] == partner ? 2offset : 0
            setlabel!(h, k, l)
        end
        old2new = invperm(collect(Int, canonize!(h)))
        push!(out, Environment(h, (old2new[pos[u]], old2new[pos[w]]), env.depth - 1, env.rules))
    end
    return out
end
