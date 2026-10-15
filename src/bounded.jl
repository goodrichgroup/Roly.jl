"""
    isbounded(rules::BindingRules; maxlength=_maxchainlength(rules) + 1)

Return `true` if every structure `rules` admits is finite, and `false` if it admits arbitrarily
large ones.

`maxlength` bounds the chains tried. Lower it only to trade certainty for speed: below the
default, `true` means no repeat was found within the budget.
"""
function isbounded(rules::BindingRules; maxlength::Integer=_maxchainlength(rules) + 1)
    return _unboundedgrowth(rules; maxlength) === nothing
end

# A rigid motion `g` and a polyform `C` such that `C, g(C), g²(C), …` is an infinite valid
# structure, as `(; generator, cell, period)`, or `nothing` if none exists within `maxlength`.
function _unboundedgrowth(rules::BindingRules; maxlength::Integer=_maxchainlength(rules) + 1)
    for i in 1:nspecies(rules)
        w = _walkchains(Polyform(rules, i), maxlength)
        w === nothing || return w
    end
    return nothing
end

# How long a chain of `rules` can get before it must repeat a state: one state per site per
# distinct twist its mates can arrive with, summed over every species.
function _maxchainlength(rules::BindingRules)
    total = 0
    for i in 1:nspecies(rules)
        ps = species(rules, i)
        for k in 1:nsites(ps)
            site = bindingsite(ps, k)
            phases = 1
            for loc in compatible_speciessitelocs(rules, color(site))
                mate = bindingsite(rules, loc)
                phases = max(phases, _ndistincttwists(mate, site))
            end
            total += phases
        end
    end
    return total
end

# Depth-first over directed chains of up to `maxlength` particles, extending only at the particle
# last added, so the states a chain has passed through are exactly the chain read from its start.
# Returns the first repeat whose motion escapes. `raise!` forms every geometric contact, so a chain
# may cross-link back onto itself; that is a valid structure, and the newest particle is extended
# regardless.
function _walkchains(poly::Polyform, maxlength)
    rules = bindingrules(poly)
    root = (speciesindex(poly.particles[1]), 0, 0)       # the first particle has no incoming bond
    stack = [(poly, [root])]
    while !isempty(stack)
        poly, states = pop!(stack)
        n = nparticles(poly)

        # the stretch between two occurrences of a state repeats under one motion
        j = findlast(==(states[end]), view(states, 1:(n - 1)))
        if !isnothing(j)
            w = _screwescape(poly, j, n)
            isnothing(w) || return w
        end
        n < maxlength || continue

        part = poly.particles[end]
        children = empty(stack)
        for k in 1:nsites(part, rules)
            site = bindingsite(part, rules, k)
            _isbound_vertex(poly, part, first(site.vertices); canonidxs=false) && continue
            isinert(rules, color(site)) && continue
            for loc in compatible_speciessitelocs(rules, color(site))
                mate = bindingsite(rules, loc)
                for r in 0:(_ndistincttwists(site, mate) - 1)
                    child = copy(poly)
                    ismissing(raise!(child, site, loc, r)) && continue
                    push!(children, (child, push!(copy(states), (loc.species, loc.site, Int(r)))))
                end
            end
        end
        append!(stack, reverse!(children))      # first child first, as the recursion had it
    end
    return nothing
end

# The stretch between two occurrences of a state repeats under `g`. Whether that goes on forever
# is decided by how far the screw climbs, and if it can, by finitely many placements.
function _screwescape(poly::Polyform, i::Integer, j::Integer)
    rules = bindingrules(poly)
    cell = poly.particles[i:(j - 1)]
    g = poly.particles[j].pose * inv(poly.particles[i].pose)

    axis, h = _screwclimb(g)
    isnothing(h) && return nothing               # a pure rotation cannot escape its own annulus

    us = (dot(p.pose.x, axis) for p in cell)
    rmax = maximum(bounding_radius(species(rules, speciesindex(p))) for p in cell)
    extent = maximum(us) - minimum(us) + 2 * rmax
    gm = g
    for _ in 1:ceil(Int, extent / abs(h))
        for p in cell
            moved = gm * p
            first(_overlap_and_contacts(cell, moved, rules)) === true && return nothing
        end
        gm = gm * g
    end
    return (; generator=g, cell=cell, period=length(cell))
end

# `(axis, climb)` of a rigid motion: the direction it carries things along, and how far per
# application. `(nothing, nothing)` when it carries them nowhere: a pure rotation in 3D, any
# rotation in 2D, or the identity.
function _screwclimb(g::Pose{D,F}) where {D,F}
    tol = sqrt(eps(F))
    # a composed motion accumulates its angle, so a full turn comes back as 2π rather than 0;
    # wrap to (-π, π] before asking whether it turns at all
    turning = abs(rem(rotation_angle(g.psi), 2 * F(π), RoundNearest)) > tol
    if !turning
        # with no turn there is no rotation axis, so the translation's own direction serves, and
        # the climb along it is the whole translation. In the plane this is the only motion that
        # escapes: a turn there has a center rather than an axis, and nothing to climb along
        n = norm(g.x)
        n > tol || return nothing, nothing        # the identity: no repeat at all
        return g.x / n, n
    end
    D == 3 || return nothing, nothing
    axis = rotation_axis(g.psi)
    climb = dot(g.x, axis)                        # the translation's part along the axis; the rest
    return axis, abs(climb) > tol ? climb : nothing    # only says where the axis sits
end
