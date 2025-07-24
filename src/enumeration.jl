function f!(k::Polyform{T,F}, s::Polyform{T,F}, assembly_system) where {T,F}
    copy!(k, s)
    lower!(k, assembly_system)
    return k
end

function adj!(u::Polyform{T,F}, v::Polyform{T,F}, j::Integer, hashes::Vector{HashType},
              assembly_system::AssemblySystem) where {T,F}
    ### TODO: We should exploit that canonical labels are guaranteed (see Nauty User Guide p.4) to be in order of color
    ### So, if we know the color of the particle that would be removed, we can filter the possible offspring
    if size(v) == 0
        bblocks = buildingblocks(assembly_system)
        if j > length(bblocks)
            return nothing
        end

        copy!(u, bblocks[j])
        push!(hashes, rhash(u))
        return u
    end

    copy!(u, v)
    success, _ = raise!(u, j, assembly_system)
    !success && return nothing

    h = rhash(u)
    if h ∉ hashes
        push!(hashes, h)
        return u
    else
        return missing
    end
end

"""
    polyenum([f], assembly_system::AssemblySystem{D,T,F}; maxsize=Inf, maxstrs=Inf, kwargs...) where {D,T,F}

Iterate over all structures (polyforms) allowed by `assembly_system` using _reverse-search_. Structures are generated up to size `maxsize`, 
and the enumeration terminates after `maxstrs` structures have been generated. Use the function `f` to process the 
generated structures and to impose constraints that the allowed structures must satisfy.

`f(s)` must take as inputs a structure `s` and must return one of three signals:

- `ACCEPT` (or `true`): enumeration continues as normal.
- `REJECT` (or `false`): the offspring of the current structure will not be generated.
- `BREAK`: the enumeration terminates immediately.

To ensure well-defined behavior, the function `f` may not accept the offspring of structures that it rejects.

All other keyword arguments are passed to the underlying `reversesearch` routine, see its documentation for more information.

Return value consist of the total number of structures generated and the largest structure size observed.
"""
function polyenum(f, assembly_system::AssemblySystem{D,T,F}; maxsize=Inf, maxstrs=Inf, kwargs...) where {D,T,F}
    ls(u, v) = f!(u, v, assembly_system)
    adj(u, v, j, aux) = adj!(u, v, j, aux, assembly_system)

    v₀ = Polyform{D,T,F}()
    aux = HashType[]
    rsys = RSSystem(ls, adj, v₀; compare=is_isomorphic, aux)

    frs = isnothing(f) ? nothing : (v, depth, args...)->f(v, size(v), args...)
    _, nv, maxdepth = reversesearch(frs, rsys; maxdepth=maxsize, maxverts=maxstrs+1, kwargs...)
    return nv-1, maxdepth
end
polyenum(assembly_system::AssemblySystem; kwargs...) = polyenum(nothing, assembly_system; kwargs...)


"""
    polygen([f::Function], assembly_system::AssemblySystem{D,T,F,G}; maxsize=Inf, maxstrs=Inf) where {D,T,F,G}

Generate all structures (polyforms) allowed by `assembly_system` using a brute force enumeration, and remove duplicates by comparing the graph
hashes of structure anatomies. This function may be faster than `polyenum` for small enumerations, but requires much more memory. Structures are 
generated up to size `maxsize`, and the enumeration terminates after `maxstrs` structures have been generated. Use the function `f` to impose additional 
constraints that the generated structures must satisfy.

`f(s)` must take as inputs a structure `s` and must return one of three signals:

- `ACCEPT` (or `true`): `s` is added to the list of structures and the enumeration continues as normal.
- `REJECT` (or `false`): `s` is not added to the list of structures and the offspring of `s` will not be generated.
- `BREAK`: `s` is not added to the list of structures and the enumeration terminates immediately.

To ensure well-defined behavior, the function `f` may not accept offspring of structures that it rejects.

Return value is a list of all generated structures.
"""
function polygen(f::Function, assembly_system::AssemblySystem{D,T,F,G}; maxsize=Inf, maxstrs=Inf) where {D,T,F,G}
    strs = [bblock for bblock in buildingblocks(assembly_system) if f(bblock) == ReverseSearch.ACCEPT]
    hashes = Set{HashType}()
    queue = Queue{Polyform{D,T,F}}()

    for bblock in strs
        hashval = rhash(bblock)
        enqueue!(queue, bblock)
        push!(hashes, hashval)
    end

    u = Polyform{D,T,F}()
    nstrs = length(strs)

    while !isempty(queue) && nstrs < maxstrs
        v = dequeue!(queue)
        n = size(v)
        if n >= maxsize
            continue
        end

        j = 1
        vert, partner_label = open_bond(v, assembly_system, j)
        while !iszero(vert)
            copy!(u, v)

            success = attach_monomer!(u, vert, partner_label, assembly_system, true)
            hashval = rhash(u)

            if success && (hashval ∉ hashes)
                push!(hashes, hashval)
                next = copy(u)
                depth = size(next)

                signal = f(next)
                if signal == ReverseSearch.ACCEPT
                    push!(strs, next)
                    enqueue!(queue, next)
                    nstrs += 1
                end
                if nstrs == maxstrs || signal == ReverseSearch.BREAK
                    break
                end
            end

            j += 1
            vert, partner_label = open_bond(v, assembly_system, j)
        end
    end

    return strs
end
function polygen(assembly_system::AssemblySystem; kwargs...)
    return polygen(_->ReverseSearch.ACCEPT, assembly_system; kwargs...)
end
