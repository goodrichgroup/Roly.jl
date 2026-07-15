mutable struct PolyformAux{BS<:BindingSite}
    seen::Set{NautyDiGraph}
    pairs::Vector{Tuple{BS,BindingSiteLoc}}
end
Base.copy(polyaux::PolyformAux) = typeof(polyaux)(copy(polyaux.seen), copy(polyaux.pairs))

function ls!(k::Polyform, s::Polyform)
    copy!(k, s)
    return lower!(k)
end

function adj!(u::Polyform, v::Polyform, j::Integer, aux::PolyformAux)
    ### TODO: exploit that canonical labels are in order of color (Nauty User Guide p.4)
    ### to filter possible offspring by the color of the particle that would be removed
    if nparticles(v) == 0
        j > nspecies(bindingrules(v)) && return nothing
        return copy!(u, Polyform(bindingrules(v), j))
    end

    j == 1 && collect_compatible_pairs!(aux.pairs, v)
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
    polyenum([f], assembly_system::BindingRules; maxsize=Inf, maxstrs=Inf, kwargs...)

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
function polyenum(f, sys::BindingRules; maxsize=Inf, maxstrs=Inf, kwargs...)
    v₀ = Polyform(sys)
    BS = BindingSite{posetype(sys), numtype(sys)}
    aux = PolyformAux{BS}(Set{NautyDiGraph}(), Tuple{BS,BindingSiteLoc}[])
    rsys = RSSystem(ls!, adj!, v₀; aux)

    frs = isnothing(f) ? nothing : (v, _) -> f(v, nparticles(v))
    result = reversesearch(frs, rsys; maxdepth=maxsize, maxverts=maxstrs+1, kwargs...)
    return result.nvertices - 1, result.depth_reached
end
polyenum(sys::BindingRules; kwargs...) = polyenum(nothing, sys; kwargs...)

function polygen(sys; kwargs...)
    v₀ = Polyform(sys)
    strs = typeof(v₀)[]
    f(s, args...) = (push!(strs, copy(s)); true)
    polyenum(f, sys; kwargs...)
    return sort!(strs; by=Roly.nparticles)
end
