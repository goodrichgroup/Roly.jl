mutable struct PolyformAux{BS<:BindingSite}
    seen::Set{NautyDiGraph}
    pairs::Vector{Tuple{BS,BindingSiteLoc}}
end
Base.copy(polyaux::PolyformAux) = typeof(polyaux)(copy(polyaux.seen), copy(polyaux.pairs))

function collect_compatible_pairs!(aux::PolyformAux, poly::Polyform)
    sys = assemblysystem(poly)
    inv_cv = invperm(canonical_vertices(poly))
    empty!(aux.pairs)

    for orig_v in canonical_vertices(poly)
        part = particle(poly, orig_v)
        isnothing(part) && continue

        for k in 1:nsites(part, sys)
            site = bindingsites(part, sys, k)
            isbound_vertex(poly, inv_cv[first(site.vertices)]) && continue
            isinert(sys, color(site)) && continue

            for siteloc in compatible_sitelocs(sys, color(site))
                push!(aux.pairs, (site, siteloc))
            end
        end
    end
    return
end

function ls!(k::Polyform, s::Polyform)
    copy!(k, s)
    return lower!(k)
end

function adj!(u::Polyform, v::Polyform, j::Integer, aux::PolyformAux)
    ### TODO: exploit that canonical labels are in order of color (Nauty User Guide p.4)
    ### to filter possible offspring by the color of the particle that would be removed
    if nparticles(v) == 0
        j > nspecies(assemblysystem(v)) && return nothing
        return copy!(u, Polyform(assemblysystem(v), j))
    end

    j == 1 && collect_compatible_pairs!(aux, v)
    j > length(aux.pairs) && return nothing

    site, siteloc = aux.pairs[j]
    copy!(u, v)
    out = raise!(u, site, siteloc)
    (ismissing(out) || isnothing(out)) && return out

    if graphrep(out) ∈ aux.seen
        return missing
    else
        push!(aux.seen, copy(graphrep(out)))
        return out
    end
end

"""
    polyenum([f], assembly_system::AssemblySystem; maxsize=Inf, maxstrs=Inf, kwargs...)

Iterate over all structures (polyforms) allowed by `assembly_system` using _reverse-search_.
Structures are generated up to size `maxsize`, and the enumeration terminates after `maxstrs`
structures have been generated. Use the function `f` to process generated structures and impose
constraints.

`f(s)` must take a structure `s` and return one of three signals:

- `ACCEPT` (or `true`): enumeration continues as normal.
- `REJECT` (or `false`): the offspring of the current structure will not be generated.
- `BREAK`: the enumeration terminates immediately.

To ensure well-defined behavior, the function `f` may not accept the offspring of structures
that it rejects.

All other keyword arguments are passed to the underlying `reversesearch` routine.

Returns `(num_structures, largest_size)`.
"""
function polyenum(f, sys::AssemblySystem; maxsize=Inf, maxstrs=Inf, kwargs...)
    v₀ = Polyform(sys)
    BS = BindingSite{posetype(sys), numtype(sys)}
    aux = PolyformAux{BS}(Set{NautyDiGraph}(), Tuple{BS,BindingSiteLoc}[])
    rsys = RSSystem(ls!, adj!, v₀; aux)

    frs = isnothing(f) ? nothing : (v, args...) -> f(v, nparticles(v), args...)
    _, nv, maxdepth = reversesearch(frs, rsys; maxdepth=maxsize, maxverts=maxstrs+1, kwargs...)
    return nv-1, maxdepth
end
polyenum(sys::AssemblySystem; kwargs...) = polyenum(nothing, sys; kwargs...)


# """
#     polygen([f::Function], assembly_system::AssemblySystem; maxsize=Inf, maxstrs=Inf)
#
# Generate all structures allowed by `assembly_system` using brute-force enumeration with
# deduplication by graph hash. Returns a list of all generated structures.
# """
# function polygen(f::Function, assembly_system::AssemblySystem; maxsize=Inf, maxstrs=Inf)
#     ...
# end
# polygen(sys::AssemblySystem; kwargs...) = polygen((_,__)->ReverseSearch.ACCEPT, sys; kwargs...)
