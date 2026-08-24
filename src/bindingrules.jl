"""
    BindingSiteLoc

Indicates the location of a binding site in the format `(speciesindex, site_index)`.
"""
const BindingSiteLoc = NTuple{2,Int}

"""
    BindingRules(bonds, particlespecies)

`BindingRules` consists of a list of particle species and their allowed bonds. 

The binding rules are specified by the matrix `bonds`, where every row of the matrix defines a valid bond between particle species.
A bonds is specified in the format `[species1 site1 species2 site2]`. For example, if binding site 3 of particle species 1 binds to site 4 of
particle species 2, the corresponding bond would read `[1 3 2 4;]`. The ordering of binding sites depends on the specific particle species.
"""
struct BindingRules{D,PS<:ParticleSpecies}
    intmat::Symmetric{Bool,Matrix{Bool}}
    particlespecies::Vector{PS}
    nbonds::Int
    nsites::Int
    ncolors::Int
    _color2siteloc::Dict{Int,Vector{BindingSiteLoc}}
    _siteloc2color::Dict{BindingSiteLoc,Int}
    _bondlist::Vector{NTuple{2,Int}}
    _bonded_sites::Vector{NTuple{2,Vector{BindingSiteLoc}}}
    _bonded_species::Vector{NTuple{2,Int}}
    _compatible_sitelocs::Vector{Vector{BindingSiteLoc}}
    _distinct_attachments::Vector{Vector{BindingSiteLoc}}
    _isinert::BitVector
    _onlattice::Bool
end
function BindingRules(bonds, particlespecies::AbstractVector{PS}) where {PS<:ParticleSpecies}
    D = dimension(first(particlespecies))
    for ps in particlespecies
        dimension(ps) != D && throw(ArgumentError("particlespecies do not all have the same dimension"))
    end

    particlespecies = _adjust_labels_and_colors(particlespecies)
    color2siteloc, siteloc2color = _make_bindingsite_lookuptables(particlespecies)

    nsites = length(siteloc2color)
    ncolors = maximum(keys(color2siteloc))

    intmat = _parse_intmat(bonds, siteloc2color, ncolors)
    bondlist = map(findall(intmat)) do cartidx
        (cartidx[1], cartidx[2])
    end
    filter!(bondlist) do bond
        bond[1] <= bond[2]
    end
    sort!(bondlist)

    bondedsites = map(bondlist) do bond
        a, b, = bond[1], bond[2]
        sort((color2siteloc[a], color2siteloc[b]))
    end
    bondedspecies = map(bondedsites) do (sitelocs1, sitelocs2)
        spc1, spc2 = sitelocs1[1][1], sitelocs2[1][1]
        (spc1, spc2)
    end

    nbonds = (sum(intmat) + sum(diagview(intmat))) ÷ 2

    compatible_sitelocs_cache = [Vector{BindingSiteLoc}() for _ in 1:ncolors]
    for c in 1:ncolors
        for c2 in 1:ncolors
            intmat[c2, c] && append!(compatible_sitelocs_cache[c], color2siteloc[c2])
        end
    end
    distinct_attachments = [_first_per_orbit(particlespecies, ls) for ls in compatible_sitelocs_cache]
    isinert_cache = BitVector(!any(intmat[:, c]) for c in 1:ncolors)

    return BindingRules{D,PS}(
        intmat,
        particlespecies,
        nbonds,
        nsites,
        ncolors,
        color2siteloc,
        siteloc2color,
        bondlist,
        bondedsites,
        bondedspecies,
        compatible_sitelocs_cache,
        distinct_attachments,
        isinert_cache,
        _onlattice(particlespecies),
    )
end
function BindingRules(bonds, particlespecies::ParticleSpecies)
    nspcs = _extract_nspecies(bonds)
    particlespecies = [particlespecies for _ in 1:nspcs]
    return BindingRules(bonds, particlespecies)
end

# A `Bool` matrix is an interaction matrix, never a bond list, so the species count has to come
# from its size rather than from the largest species index a bond names.
function BindingRules(intmat::AbstractMatrix{Bool}, particlespecies::ParticleSpecies)
    n = size(intmat, 1)
    n == size(intmat, 2) || throw(ArgumentError("interaction matrix must be square"))
    k = _colorspan(particlespecies)
    n % k == 0 || throw(
        ArgumentError(
            "an interaction matrix of size $n cannot be split among copies of a species " *
            "spanning $k colors; give one species per block instead",
        ),
    )
    return BindingRules(intmat, [particlespecies for _ in 1:(n ÷ k)])
end

# How many columns of the interaction matrix one copy of `ps` claims, matching the contiguous
# range `_adjust_labels_and_colors` gives it.
function _colorspan(ps::ParticleSpecies)
    cols = [color(bindingsite(ps, si)) for si in 1:nsites(ps)]
    return maximum(cols) - minimum(cols) + 1
end

"""
    _onlattice(pss)

Whether every species in `pss` tiles, all with the same cell. Allows the coincident-centers
shortcut in `overlap(::Particle, ::Particle, ::BindingRules)`. See [`_tilingcell`](@ref).
"""
function _onlattice(pss::AbstractVector{<:ParticleSpecies})
    cells = map(_tilingcell, pss)
    any(isnothing, cells) && return false
    sides, edge = first(cells)
    return all(c -> c[1] == sides && isapprox(c[2], edge; rtol=sqrt(eps(float(typeof(edge))))), cells)
end

"""
    interactionmatrix(rules::BindingRules)

Return the interaction matrix of the assembly system `rules`.

The interaction matrix is indexed by the colors of binding sites.
If two binding sites have colors `c1` and `c2`, then the truth value of 
`interactionmatrix(rules)[c1, c2]` determines whether the two sites are able
to bind to each other.
"""
@inline interactionmatrix(rules::BindingRules) = rules.intmat

"""
    species(rules::BindingRules)

Return the list of particle species of the assembly system `rules`.
"""
@inline species(rules::BindingRules) = rules.particlespecies

"""
    species(rules::BindingRules, i::Integer)

Return the `i`th particle species of the assembly system `rules`.
"""
@inline species(rules::BindingRules, i::Integer) = rules.particlespecies[i]

"""
    nspecies(rules::BindingRules)

Return the number of particle species of the assembly system `rules`.
"""
@inline nspecies(rules::BindingRules) = length(rules.particlespecies)

"""
    nbonds(rules::BindingRules)

Return the total number of possible bonds between binding sites of the assembly system `rules`.
"""
@inline nbonds(rules::BindingRules) = rules.nbonds

"""
    nsites(rules::BindingRules)

Return the total number of binding sites across all particle species of the assembly system `rules`.
"""
@inline nsites(rules::BindingRules) = rules.nsites

"""
    ncolors(rules::BindingRules)

Return the total number of binding sites colors across all particle species of the assembly system `rules`.

`ncolors(rules)` is not necessarily the equal to `nsites(rules)`. 
If particle species have symmetry, some binding sites may carry the same color.
"""
@inline ncolors(rules::BindingRules) = rules.ncolors

"""
    dimension(::BindingRules)

Return the spatial dimension of the particles of the assembly system `rules`.
"""
@inline dimension(::BindingRules{D}) where {D} = D
@inline dimension(::Type{<:BindingRules{D}}) where {D} = D

@inline numtype(::BindingRules{D,PS}) where {D,PS} = numtype(PS)
@inline posetype(::BindingRules{D,PS}) where {D,PS} = posetype(PS)
@inline speciestype(::BindingRules{D,PS}) where {D,PS} = PS

"""
    bonded_colors(rules::BindingRules)

Return all pairs of binding site colors that may bind according to the binding rules
of the assembly system `rules`.
"""
@inline function bonded_colors(rules::BindingRules)
    return rules._bondlist
end

"""
    bonded_sites(rules::BindingRules)

Return all pairs of binding site locations that may bind according to the binding rules
of the assembly system `rules`.
"""
@inline function bonded_sites(rules::BindingRules)
    return rules._bonded_sites
end

"""
    bonded_species(rules::BindingRules)

Return all pairs of particle species that may bind according to the binding rules
of the assembly system `rules`.
"""
@inline function bonded_species(rules::BindingRules)
    return rules._bonded_species
end

"""
    siteloc2color(rules::BindingRules, siteloc::BindingSiteLoc)

Return the binding site color associated with the binding site location `siteloc`.
"""
@inline function siteloc2color(rules::BindingRules, siteloc::BindingSiteLoc)
    return rules._siteloc2color[siteloc]
end

"""
    color2siteloc(rules::BindingRules, color::Integer)

Return the (possible multiple) binding site locations associated with the binding site color `color`.
"""
@inline function color2siteloc(rules::BindingRules, color::Integer)
    return rules._color2siteloc[color]
end

"""
    color2species(rules::BindingRules, color::Integer)

Return the particle species that contains the binding site with color `color`.
"""
@inline function color2species(rules::BindingRules, color::Integer)
    return color2siteloc(rules, color)[1][1]
end

"""
    _first_per_orbit(particlespecies, sitelocs)

Return one site per symmetry orbit of the particle species. Of the entries of `sitelocs` that
sit on the same species and carry the same graph label, only the first survives.
"""
function _first_per_orbit(
    particlespecies::AbstractVector{<:ParticleSpecies}, sitelocs::AbstractVector{BindingSiteLoc}
)
    reps = BindingSiteLoc[]
    seen = Set{Tuple{Int,Int}}()
    for (spc, k) in sitelocs
        orbit = (spc, sitelabel(particlespecies[spc], k))
        orbit in seen && continue
        push!(seen, orbit)
        push!(reps, (spc, k))
    end
    return reps
end

"""
    distinct_attachments(rules::BindingRules, color::Integer)

The sites a particle may be attached through to a site of `color`, one per symmetry orbit.

[`compatible_sitelocs`](@ref) contains every site the interaction matrix permits, this only lists the
ones that lead to distinguishable structures. See [`_first_per_orbit`](@ref).
"""
@inline distinct_attachments(rules::BindingRules, color::Integer) = rules._distinct_attachments[color]

"""
    isinert(rules::BindingRules, siteloc::BindingSiteLoc)

Return `true` if the binding site at `siteloc` does not bind to any binding sites.
"""
@inline function isinert(rules::BindingRules, siteloc::BindingSiteLoc)
    color = siteloc2color(rules, siteloc)
    return isinert(rules, color)
end

"""
    isinert(rules::BindingRules, color::Integer)

Return `true` if the binding site with color `color` does not bind to any binding sites.
"""
@inline function isinert(rules::BindingRules, color::Integer)
    return rules._isinert[color]
end

"""
    graphrep(rules::BindingRules)

Construct the graph representation of the ass embly system `rules`.
"""
function graphrep(rules::BindingRules)
    imat = interactionmatrix(rules)
    ns = nsites(rules)
    nb = nbonds(rules)

    A = zeros(Int, ns + nb, ns + nb)
    edge_counter = 1
    for i in axes(imat, 1), j in i:size(imat, 2)
        if imat[i, j] == 0
            continue
        end

        A[i, ns + edge_counter] = 1
        A[j, ns + edge_counter] = 1
        A[ns + edge_counter, i] = 1
        A[ns + edge_counter, j] = 1
        edge_counter += 1
    end

    i = 1
    for bb in species(rules)
        a = adjacency_matrix(graphrep(bb))
        n = size(a, 1)
        A[i:(i + n - 1), i:(i + n - 1)] .+= a
        i += n
    end

    vertex_labels = [zeros(Cint, ns); ones(Cint, nb)]
    g = NautyDiGraph(A; vertex_labels)
    return g
end

function compatible_sitelocs(rules::BindingRules, color::Integer)
    return rules._compatible_sitelocs[color]
end

function Base.show(io::Core.IO, rules::BindingRules)
    print(io, "$(dimension(rules))d BindingRules[n=$(nspecies(rules)), k=$(nbonds(rules))]")
end

function _extract_nspecies(bonds::AbstractMatrix{<:Integer})
    _checkshape(bonds)
    length(bonds) == 0 &&
        throw(ArgumentError("binding rules do not not contain any bonds and cannot be used to determine species count"))
    return maximum(@view bonds[:, [1, 3]]; init=0)
end

"""
    _adjust_labels_and_colors(particlespecies)

Shift each species' colors and graph labels into disjoint global ranges, so that
species built independently do not collide and a site's color indexes the interaction matrix
directly.
"""
function _adjust_labels_and_colors(particlespecies::AbstractVector{PS}) where {PS<:ParticleSpecies}
    particlespecies = PS[copy(ps) for ps in particlespecies]
    c = 1
    l = 1
    for ps in particlespecies
        cols = [color(bindingsite(ps, si)) for si in 1:nsites(ps)]
        setcolors!(ps, cols .- minimum(cols) .+ c)
        c = maximum(cols) - minimum(cols) + c + 1

        g = graphrep(ps)
        offset = l - minimum(labels(g))
        for v in vertices(g)
            setlabel!(g, v, label(g, v) + offset)
        end
        l = maximum(labels(g)) + 1
    end
    return particlespecies
end

function _make_bindingsite_lookuptables(particlespecies::AbstractVector{<:ParticleSpecies})
    color2siteloc = Dict{Int,Vector{BindingSiteLoc}}()
    siteloc2color = Dict{BindingSiteLoc,Int}()

    for (spcs, ps) in enumerate(particlespecies)
        for si in 1:nsites(ps)
            site = bindingsite(ps, si)
            c = color(site)
            spcssite = (spcs, si)
            push!(get!(color2siteloc, c, BindingSiteLoc[]), spcssite)
            siteloc2color[spcssite] = c
        end
    end
    return color2siteloc, siteloc2color
end

function _parse_intmat(bonds, siteloc2color, ncolors)
    return _intmat_from_bonds(bonds, siteloc2color, ncolors)
end
function _parse_intmat(bonds::AbstractMatrix{<:Integer}, siteloc2color, ncolors)
    _checkshape(bonds)
    return _intmat_from_bonds(eachrow(bonds), siteloc2color, ncolors)
end

function _parse_intmat(intmat::AbstractMatrix{Bool}, siteloc2color, ncolors)
    n = size(intmat, 1)
    n == size(intmat, 2) || throw(ArgumentError("interaction matrix must be square"))
    n == ncolors || throw(ArgumentError("interaction matrix size ($n) does not match the number of colors ($ncolors)"))
    intmat == intmat' || throw(ArgumentError("interaction matrix must be symmetric"))
    return Symmetric(Matrix{Bool}(intmat))
end
function _intmat_from_bonds(bonds, siteloc2color, ncolors)
    intmat = zeros(Bool, ncolors, ncolors)
    for bond in bonds
        spcs1, site1, spcs2, site2 = bond
        c1 = siteloc2color[(spcs1, site1)]
        c2 = siteloc2color[(spcs2, site2)]
        intmat[c1, c2] = intmat[c2, c1] = true
    end
    return Symmetric(intmat)
end

function _checkshape(bonds::AbstractMatrix{<:Integer})
    mustbe4 = size(bonds, 2)
    mustbe4 != 4 && throw(ArgumentError("bonds must be defined in the form [spcs1 site1 spcs2 site2]"))
    return nothing
end
function _checkshape(bonds::AbstractVector{AbstractVector{<:Integer}})
    mustbe4 = length(bonds[1])
    mustbe4 != 4 && throw(ArgumentError("bonds must be defined in the form [spcs1 site1 spcs2 site2]"))
    return nothing
end
