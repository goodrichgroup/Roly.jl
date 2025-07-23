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
            return nothing, 1
        end

        copy!(u, bblocks[j])
        push!(hashes, rhash(u))
        return u, 1
    end

    copy!(u, v)
    success, j_new = raise!(u, j, assembly_system)

    while success && rhash(u) ∈ hashes
        copy!(u, v)
        success, j_new = raise!(u, j_new + 1, assembly_system)
    end

    success && push!(hashes, rhash(u))
    return success ? u : nothing, j_new - j
end

function polyenum(f, assembly_system::AssemblySystem{D,T,F}; maxsize=Inf, maxstrs=Inf, cached=true) where {D,T,F}
    ls(u, v) = f!(u, v, assembly_system)
    adj(u, v, j, aux) = adj!(u, v, j, aux, assembly_system)

    rsys = RSSystem(ls, adj, ≃)
    v₀ = Polyform{D,T,F}()
    aux = HashType[]

    _, nv, maxdepth = reversesearch(f, rsys, v₀; cached, aux, maxdepth=maxsize, maxverts=maxstrs+1)
    return nv-1, maxdepth
end
polyenum(assembly_system::AssemblySystem; kwargs...) = polyenum(nothing, assembly_system; kwargs...)


function polygen(f::Function, assembly_system::AssemblySystem{D,T,F,G}; maxsize=Inf, maxstrs=Inf) where {D,T,F,G}
    bblocks = buildingblocks(assembly_system)

    values = [f(bblock) for bblock in bblocks]
    hashes = Set{HashType}()
    queue = Queue{Polyform{D,T,F}}()

    for bblock in bblocks
        hashval = rhash(bblock)
        enqueue!(queue, bblock)
        push!(hashes, hashval)
    end

    u = Polyform{D,T,F}()
    n_strs = length(bblocks)

    while !isempty(queue) && n_strs < maxstrs
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
                next = copy(u)
            
                push!(values, f(next))
                push!(hashes, hashval)
                enqueue!(queue, next)
                
                n_strs += 1
                n_strs == maxstrs && break
            end

            j += 1
            vert, partner_label = open_bond(v, assembly_system, j)
        end
    end

    return values
end
function polygen(assembly_system::AssemblySystem; kwargs...)
    return polygen(identity, assembly_system; kwargs...)
end
