"""
    BindingSiteLoc

Indicates the location of a binding site in the format `(species_index, site_index)`.
"""
const BindingSiteLoc = NTuple{2,Int}

"""
    BindingRules(bindingrules, particlespecies)

An `BindingRules` consists of a list of particle species and their binding rules. The structures (polyforms)
that satify the binding rules without any particle overlaps are the _allowed structures_ of an assembly system, and can be enumerated or 
iterated over using `polygen` or `polyenum`, respectively.

The binding rules are specified by the matrix `bindingrules`, where every row of the matrix defines a valid bond between particle species.
A bonds is specified in the format `[species1 site1 species2 site2]`. For example, if binding site 3 of particle species 1 binds to site 4 of
particle species 2, the corresponding bond would read `[1 3 2 4;]`. The ordering of binding sites depends on the building block, for 
polygonal building blocks, sites are numbered in clockwise order.

Even without enumerating the allowed structures, it is possible to determine whether two assembly systems lead to the same set of allowed structures.
If this is the case, then the two systems are identical up a relabeling and/or rotating the building blocks -- they are isomorphic.
This can be checked by comparing the `rhash`es of different assembly systems; if the hashes are equal, the systems are isomorphic. 
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
    _isinert::BitVector
end
function BindingRules(bindingrules, particlespecies::AbstractVector{PS}) where {PS<:ParticleSpecies}
    D = dimension(first(particlespecies))
    for ps in particlespecies
        dimension(ps) != D && throw(ArgumentError("particlespecies do not all have the same dimension"))
    end

    particlespecies = _adjust_labels_and_colors(particlespecies)
    color2siteloc, siteloc2color = _make_bindingsite_lookuptables(particlespecies)

    nsites = length(siteloc2color)
    ncolors = length(color2siteloc)

    intmat = _parse_intmat(bindingrules, siteloc2color)
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
            intmat[c2, c] || continue
            append!(compatible_sitelocs_cache[c], color2siteloc[c2])
        end
    end
    isinert_cache = BitVector(!any(intmat[:, c]) for c in 1:ncolors)

    return BindingRules{D,PS}(intmat, particlespecies, nbonds, nsites, ncolors,
                                color2siteloc, siteloc2color,
                                bondlist, bondedsites, bondedspecies,
                                compatible_sitelocs_cache, isinert_cache)
end
function BindingRules(bindingrules, particlespecies::ParticleSpecies)
    nspcs = _extract_nspecies(bindingrules)
    particlespecies = [particlespecies for _ in 1:nspcs]
    return BindingRules(bindingrules, particlespecies)
end


"""
    interactionmatrix(sys::BindingRules)

Return the interaction matrix of the assembly system `sys`.

The interaction matrix is indexed by the colors of binding sites.
If two binding sites have colors `c1` and `c2`, then the truth value of 
`interactionmatrix(sys)[c1, c2]` determines whether the two sites are able
to bind to each other.
"""
@inline interactionmatrix(sys::BindingRules) = sys.intmat

"""
    species(sys::BindingRules)

Return the list of particle species of the assembly system `sys`.
"""
@inline species(sys::BindingRules) = sys.particlespecies

"""
    species(sys::BindingRules, i::Integer)

Return the `i`th particle species of the assembly system `sys`.
"""
@inline species(sys::BindingRules, i::Integer) = sys.particlespecies[i]

"""
    nspecies(sys::BindingRules)

Return the number of particle species of the assembly system `sys`.
"""
@inline nspecies(sys::BindingRules) = length(sys.particlespecies)

"""
    nbonds(sys::BindingRules)

Return the total number of possible bonds between binding sites of the assembly system `sys`.
"""
@inline nbonds(sys::BindingRules) = sys.nbonds

"""
    nsites(sys::BindingRules)

Return the total number of binding sites across all particle species of the assembly system `sys`.
"""
@inline nsites(sys::BindingRules) = sys.nsites

"""
    ncolors(sys::BindingRules)

Return the total number of binding sites colors across all particle species of the assembly system `sys`.

`ncolors(sys)` is not necessarily the equal to `nsites(sys)`. 
If particle species have symmetry, some binding sites may carry the same color.
"""
@inline ncolors(sys::BindingRules) = sys.ncolors

"""
    dimension(::BindingRules)

Return the spatial dimension of the particles of the assembly system `sys`.
"""
@inline dimension(::BindingRules{D}) where {D} = D
@inline dimension(::Type{<:BindingRules{D}}) where {D} = D

@inline numtype(::BindingRules{D,PS}) where {D,PS} = numtype(PS)
@inline posetype(::BindingRules{D,PS}) where {D,PS} = posetype(PS)
@inline speciestype(::BindingRules{D,PS}) where {D,PS} = PS

"""
    bonded_colors(sys::BindingRules)

Return all pairs of binding site colors that may bind according to the binding rules
of the assembly system `sys`.
"""
@inline function bonded_colors(sys::BindingRules)
    return sys._bondlist
end

"""
    bonded_sites(sys::BindingRules)

Return all pairs of binding site locations that may bind according to the binding rules
of the assembly system `sys`.
"""
@inline function bonded_sites(sys::BindingRules)
    return sys._bonded_sites
end

"""
    bonded_species(sys::BindingRules)

Return all pairs of particle species that may bind according to the binding rules
of the assembly system `sys`.
"""
@inline function bonded_species(sys::BindingRules)
    return sys._bonded_species
end

"""
    siteloc2color(sys::BindingRules, siteloc::BindingSiteLoc)

Return the binding site color associated with the binding site location `siteloc`.
"""
@inline function siteloc2color(sys::BindingRules, siteloc::BindingSiteLoc)
    return sys._siteloc2color[siteloc]
end

"""
    color2siteloc(sys::BindingRules, color::Integer)

Return the (possible multiple) binding site locations associated with the binding site color `color`.
"""
@inline function color2siteloc(sys::BindingRules, color::Integer)
    return sys._color2siteloc[color]
end

"""
    color2species(sys::BindingRules, color::Integer)

Return the particle species that contains the binding site with color `color`.
"""
@inline function color2species(sys::BindingRules, color::Integer)
    return color2siteloc(sys, color)[1][1]
end

"""
    isinert(sys::BindingRules, siteloc::BindingSiteLoc)

Return `true` if the binding site at `siteloc` does not bind to any binding sites.
"""
@inline function isinert(sys::BindingRules, siteloc::BindingSiteLoc)
    color = siteloc2color(sys, siteloc)
    return isinert(sys, color)
end

"""
    isinert(sys::BindingRules, color::Integer)

Return `true` if the binding site with color `color` does not bind to any binding sites.
"""
@inline function isinert(sys::BindingRules, color::Integer)
    return sys._isinert[color]
end

"""
    graphrep(sys::BindingRules)

Construct the graph representation of the assembly system `sys`.
"""
function graphrep(sys::BindingRules)
    imat = interactionmatrix(sys)
    ns = nsites(sys)
    nb = nbonds(sys)

    A = zeros(Int, ns + nb, ns + nb)
    edge_counter = 1
    for i in axes(imat, 1), j in i:size(imat, 2)
        if imat[i, j] == 0
            continue
        end

        A[i, ns+edge_counter] = 1
        A[j, ns+edge_counter] = 1
        A[ns+edge_counter, i] = 1
        A[ns+edge_counter, j] = 1
        edge_counter += 1
    end

    i = 1
    for bb in species(sys)
        a = adjacency_matrix(graphrep(bb))
        n = size(a, 1)
        A[i:i+n-1, i:i+n-1] .+= a
        i += n
    end

    vertex_labels = [zeros(Cint, ns); ones(Cint, nb)]
    g = NautyDiGraph(A; vertex_labels)
    return g
end

function compatible_sitelocs(sys::BindingRules, color::Integer)
    return sys._compatible_sitelocs[color]
end

Base.show(io::Core.IO, sys::BindingRules) = print(io, "$(dimension(sys))d BindingRules[n=$(nspecies(sys)), k=$(nbonds(sys))]")


function _extract_nspecies(bindingrules::AbstractMatrix{<:Integer})
    _checkshape(bindingrules)
    nspcs = maximum(@view bindingrules[:, [1, 3]]; init=0)
    nspcs < 1 && throw(ArgumentError("binding rules do not not contain any bonds and cannot be used to determine species count"))
    return nspcs
end

function _adjust_labels_and_colors(particlespecies::AbstractVector{PS}) where {PS<:ParticleSpecies}
    particlespecies = PS[copy(ps) for ps in particlespecies]
    c = 1
    l = 1
    for ps in particlespecies
        g = graphrep(ps)
        labs = labels(g)
        setlabels!(g, labs .- minimum(labs) .+ l)
        l = maximum(labels(g)) + 1

        setcolors!(ps, collect(c:c+nsites(ps)-1))
        c += nsites(ps)
    end
    return particlespecies
end

function _make_bindingsite_lookuptables(particlespecies::AbstractVector{<:ParticleSpecies})
    color2siteloc = Dict{Int,Vector{BindingSiteLoc}}()
    siteloc2color = Dict{BindingSiteLoc,Int}()

    for (spcs, ps) in enumerate(particlespecies)
        for si in 1:nsites(ps)
            site = bindingsites(ps, si)
            c = color(site)
            spcssite = (spcs, si)
            push!(get!(color2siteloc, c, BindingSiteLoc[]), spcssite)
            siteloc2color[spcssite] = c
        end
    end
    return color2siteloc, siteloc2color
end

function _parse_intmat(bindingrules, siteloc2color)
    return _intmat_from_bonds(bindingrules, siteloc2color)
end
function _parse_intmat(bindingrules::AbstractMatrix{<:Integer}, siteloc2color)
    _checkshape(bindingrules)
    return _intmat_from_bonds(eachrow(bindingrules), siteloc2color)
end
function _intmat_from_bonds(bonds, siteloc2color)
    ncolors = length(siteloc2color)
    intmat = falses(ncolors, ncolors)
    for bond in bonds
        spcs1, site1, spcs2, site2 = bond
        c1 = siteloc2color[(spcs1, site1)]
        c2 = siteloc2color[(spcs2, site2)]
        intmat[c1, c2] = intmat[c2, c1] = true
    end
    return Symmetric(intmat)
end

function _checkshape(bindingrules::AbstractMatrix{<:Integer})
    mustbe4 = size(bindingrules, 2)
    mustbe4 != 4 && throw(ArgumentError("bonds must be defined in the form [spcs1 site1 spcs2 site2]"))
    return
end
function _checkshape(bindingrules::AbstractVector{AbstractVector{<:Integer}})
    mustbe4 = length(bindingrules[1])
    mustbe4 != 4 && throw(ArgumentError("bonds must be defined in the form [spcs1 site1 spcs2 site2]"))
    return
end
