
function open_bond(p::Polyform, assembly_system::AssemblySystem, j::Integer)
    for v in Iterators.reverse(p.canonical_order)
        !iszero(p.bond_partners[v]) && continue
        label = vertex2label(p, assembly_system, v)
        partner_label, Δj = @views find_nth(!iszero, assembly_system.intmat[:, label], j)
        j -= Δj
        if isnothing(partner_label)
            continue
        else
            return v, partner_label
        end
    end

    return 0, 0
end

function get_sitepos(p::Polyform, assembly_system, v)
    particle, side = vertex2particle(p, assembly_system, v)
    spc = species(p)[particle]
    geom = assembly_system.geometries[spc]
    return p.xs[particle] + rotate(geom.xs[side], p.ψs[particle])
end

function attach_buildingblock!(p::Polyform{D}, v::Integer, attached_sitelabel::Integer; canonize::Bool=false) where {D}
    !iszero(p._bonds[v]) && return false
    sys = assemblysystem(p)
    nvp = nvertices(p)

    bbspcs, bbsite = label2siteloc(sys, attached_sitelabel)
    bblock = buildingblocks(sys, bbspcs)

    bbxs, bbψs = attached_coordinates(p, bblock, v, bbsite)
    overlap, extrabonds = overlap_or_extrabonds(p, bblock, assembly_system; xs2=xs2, ψs2=ψs2)
    overlap && return false

    append!(p.bonds, zeros(T, nvertices(bblock)))
    blockdiag!(p.anatomy, bblock.anatomy)
    for (; vi, vj) in extrabonds
        vj += nvp

        add_edge!(p.anatomy, vi, vj)
        add_edge!(p.anatomy, vj, vi)
        p.bonds[vi] = vj
        p.bonds[vj] = vi
    end

    append!(p.xs, xs2)
    append!(p.ψs, ψs2)
    resize!(p.canonical_order, nvp + nvertices(bblock))

    if canonize
        canonize!(p)
    else
        p.σ = 0
        p.canonical_order .= 1:nvertices(p)
    end
    return true
end

function attached_coordinates(p::Polyform, buildingblock::BuildingBlock, v, w)
    #TODO this can probably just operate on coordinates directly, no need to pass polyforms/buildingblocks
    x0 = positions(p)[v]
    ψ0 = orientations(p)[v]

    xs = copy(sites(buildingblock).xs)
    Δx = xs[w] + x0

    ψs = copy(sites(buildingblock).ψs)
    Δψ = inv(ψs[w]) * standard_offset(ψ0) * ψ0

    xs .-= Ref(Δx)
    xs .= Ref(Δψ) .* xs
    ψs .*= Ref(Δψ)
    return xs, ψs
end

# function get_attached_coordinates(p1::Polyform, p2::Polyform, assembly_system::AssemblySystem, v1, v2)
#     geoms = geometries(assembly_system)
#     spcs1 = species(p1)
#     spcs2 = species(p2)

#     part1, site1 = vertex2particle(p1, assembly_system, v1)
#     part2, site2 = vertex2particle(p2, assembly_system, v2)

#     Δx, Δψ = attachment_offset(site1, site2, geoms[spcs1[part1]], geoms[spcs2[part2]])
#     Δx = p1.xs[part1] + rotate(Δx, p1.ψs[part1])
#     Δψ = Δψ * p1.ψs[part1]

#     xs2, ψs2 = copy(p2.xs), copy(p2.ψs)
#     grab!(xs2, ψs2, part2, Δx, Δψ)
#     return xs2, ψs2
# end

"""
    check_contacts(xs1, ψs1, xs2, ψs2; intmat, allow_noninteracting=false, allow_misaligned=false, kwargs...)

Take two sets of binding site coordinates and find all binding site contacts. Return early if any invalid
contacts (noninteracting or misalgined sites, depending on kwargs) are found. Return final status (all valid contacts?)
and a list of contacting binding site indices. The trailing kwargs are passed to `isapprox`.
"""
# function check_contacts(xs1, ψs1, xs2, ψs2; intmat, allow_noninteracting=false, allow_misaligned=false, kwargs...)
#     contacts = NTuple{2,Int}[]

#     for i in eachindex(xs1, ψs1), j in eachindex(xs2, ψs2)
#         x_overlap = isapprox(xs1[i], xs2[j]; kwargs...)
#         if x_overlap
#             si, sj = ...
#             interacting = intmat[si, siteloc2index(, j)]
#             if !allow_noninteracting && !interacting
#                 return false, contacts
#             end

#             ψ_overlap = isapprox(ψs1[j], ψs2[j]; kwargs...)
#             if !allow_misaligned && !ψ_overlap
#                 return false, contacts
#             end
#             push!(contacts, (i, j))
#         end
#     end
#     return true, contacts
# end

function overlap_or_extrabonds(p1::Polyform, p2; xs1=nothing, ψs1=nothing, xs2=nothing, ψs2=nothing)
    geoms = geometries(assembly_system)
    intmat = interactionmatrix(assembly_system)
    
    xs1 = isnothing(xs1) ? p1.xs : xs1
    ψs1 = isnothing(ψs1) ? p1.ψs : ψs1
    xs2 = isnothing(xs2) ? p2.xs : xs2
    ψs2 = isnothing(ψs2) ? p2.ψs : ψs2

    extrabonds = Edge{Cint}[]
    for (i, (xi, ψi)) in enumerate(zip(xs1, ψs1)), (j, (xj, ψj)) in enumerate(zip(xs2, ψs2))
        cstat, bound_sites = contact_status(xj - xi, ψi, ψj, geoms[spcs1[i]], geoms[spcs2[j]])
        # Return if particles overlap
        !cstat && return false, extrabonds

        for sites in bound_sites
            si, sj = sites

            vs_i = particle2multivertex(p1, assembly_system, i, si)
            vs_j = particle2multivertex(p2, assembly_system, j, sj)

            # label_i = p1.anatomy.labels[first(vs_i)]
            # label_j = p2.anatomy.labels[first(vs_j)]
            # TODO: the above causes errors if the color label is different from the site site idx
            # i.e. when some sites have the same label
            label_i = vertex2label(p1, assembly_system, first(vs_i))
            label_j = vertex2label(p2, assembly_system, first(vs_j))
            
            # Return if a disallowed bond is found
            !intmat[label_i, label_j] && return false, extrabonds

            bond_to_edge!(extrabonds, vs_i, vs_j)
        end
    end
    return true, extrabonds
end

function bond_to_edge!(edges::AbstractVector{<:Edge}, vs_i::AbstractVector{<:Integer}, vs_j::AbstractVector{<:Integer})
    ni, nj = length(vs_i), length(vs_j)
    shared_symmetry = max(ni, nj) // min(ni, nj)

    if (ni == nj == 1) || !isinteger(shared_symmetry)
        edge = eltype(edges)(first(vs_i), first(vs_j))
        push!(edges, edge)
        return
    end

    sym = numerator(shared_symmetry)
    if ni <= nj
        vs1, vs2 = vs_i, reverse(vs_j[1:sym:end])
    else
        vs1, vs2 = vs_i[1:sym:end], reverse(vs_j)
    end

    for (v1, v2) in zip(vs1, vs2)
        edge = eltype(edges)(v1, v2)
        push!(edges, edge)
    end
    return
end

function raise!(p::Polyform{D,T,F}, k::Integer, assembly_system::AssemblySystem) where {D,T,F}
    v, partner_label = open_bond(p, assembly_system, k)
    iszero(v) && return false, k

    success = attach_monomer!(p, v, partner_label, assembly_system)
    !success && return raise!(p, k+1, assembly_system)

    canonize!(p)
    return true, k
end
    

function lower!(p::Polyform, assembly_system::AssemblySystem)
    n = size(p)
    n == 0 && return false

    k = 1
    if n > 1
        for l in Iterators.reverse(p.canonical_order)
            part, _ = vertex2particle(p, assembly_system, l)
            vertices = particle2multivertex(p, assembly_system, part)
            if !are_bridge(p.anatomy, vertices)
                k = l
                break
            end
        end
    end

    del_part, _ = vertex2particle(p, assembly_system, k)
    del_vs = particle2multivertex(p, assembly_system, del_part)
    sort!(del_vs)

    deleteat!(p.xs, del_part)
    deleteat!(p.ψs, del_part)
    deleteat!(p.species, del_part)

    # open up bonds that are now free
    partners = p.bond_partners[del_vs]
    for part in partners
        if part == 0
            continue
        end
        p.bond_partners[part] = 0
    end
    deleteat!(p.bond_partners, del_vs)
    
    vertex_shift(v) = sum(x -> x <= v, del_vs)
    @views p.bond_partners .-= vertex_shift.(p.bond_partners)

    rem_vertices!(p.anatomy, del_vs)
    deleteat!(p.canonical_order, del_vs)
    canonize!(p)
    return true
end