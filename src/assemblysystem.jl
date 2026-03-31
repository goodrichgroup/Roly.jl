const BindingSiteLoc = NTuple{2,Int}

"""
    AssemblySystem(bindingrules, buildingblocks)

Defines an _assembly system_, consisting of a list of building blocks, their geometries, and their binding rules. The structures (polyforms)
that satify the binding rules without any particle overlaps are the _allowed structures_ of an assembly system, and can be enumerated or 
iterated over using `polygen` or `polyenum`, respectively.

The binding rules are specified by the matrix `bindingrules`, where every row of the matrix defines a valid bond between particle species.
A bonds is specified in the format `[species1 site1 species2 site2]`. For example, if binding site 3 of particle species 1 binds to site 4 of
particle species 2, the corresponding bond would read `[1 3 2 4;]`. The ordering of binding sites depends on the particle geometry, for 
polygonal geometries, sites are numbered in clockwise order.

Even without enumerating the allowed structures, it is possible to determine whether two assembly systems lead to the same set of allowed structures.
If this is the case, then the two systems are identical up a relabeling and/or rotating the building blocks -- they are isomorphic.
This can be checked by comparing the `rhash`es of different assembly systems; if the hashes are equal, the systems are isomorphic. 
"""
mutable struct AssemblySystem{D,BB<:BuildingBlock}
    intmat::Symmetric{Bool, Matrix{Bool}}
    buildingblocks::Vector{BB}
    nspcs::Int
    nbonds::Int
    nsites::Int
    _idx2siteloc::Vector{BindingSiteLoc}
    _siteloc2idx::Dict{BindingSiteLoc,Int}
    _label2species::Dict{Int,Int}
    _bondlist::Vector{NTuple{2,Int}}
    _bonded_sites::Vector{NTuple{2,BindingSiteLoc}}
    _bonded_species::Vector{NTuple{2,Int}}
end
function AssemblySystem(bindingrules, buildingblocks::AbstractVector{BB}) where {BB<:BuildingBlock}
    D = dimension(first(buildingblocks))
    for bb in buildingblocks
        dimension(bb) != D && throw(ArgumentError("building blocks do not all have the same dimension"))
    end

    buildingblocks = _shift_sitelabels(buildingblocks)
    idx2siteloc, siteloc2idx, label2species = _make_lookuptables(buildingblocks)

    nsites = length(siteloc2idx)
    nspcs = length(buildingblocks)

    intmat = _parse_intmat(bindingrules, siteloc2idx)
    bondlist = map(findall(intmat)) do cartidx
        (cartidx[1], cartidx[2])
    end
    filter!(bondlist) do bond
        bond[1] <= bond[2]
    end
    sort!(bondlist)
    bondedsites = map(bondlist) do bond
        a, b, = bond[1], bond[2]
        sort((idx2siteloc[a], idx2siteloc[b]))
    end
    bondedspecies = map(bondedsites) do ((spc1, _), (spc2, ))
        (spc1, spc2)
    end

    nbonds = (sum(intmat) + sum(diagview(intmat))) ÷ 2
    return AssemblySystem{D,BB}(intmat, buildingblocks, nspcs, nbonds, nsites, idx2siteloc, siteloc2idx, label2species, bondlist, bondedsites, bondedspecies)
end
function AssemblySystem(bindingrules, buildingblock::BuildingBlock)
    nspcs = _extract_nspecies(bindingrules)
    buildingblocks = [buildingblock for _ in 1:nspcs]
    return AssemblySystem(bindingrules, buildingblocks)
end

intmat(sys::AssemblySystem) = sys.intmat
buildingblocks(sys::AssemblySystem) = sys.buildingblocks

nspecies(sys::AssemblySystem) = sys.nspcs
nbonds(sys::AssemblySystem) = sys.nbonds
nsites(sys::AssemblySystem) = sys.nsites

dimension(::AssemblySystem{D}) where {D} = D
dimension(::Type{<:AssemblySystem{D}}) where {D} = D
numtype(::AssemblySystem{D,BB}) where {D,BB} = numtype(BB)
numtype(::Type{<:AssemblySystem{D,BB}}) where {D,BB} = numtype(BB)

positiontype(::AssemblySystem{D,BB}) where {D,BB} = positiontype(BB)
positiontype(::Type{<:AssemblySystem{D,BB}}) where {D,BB} = positiontype(BB)
orientationtype(::AssemblySystem{D,BB}) where {D,BB} = orientationtype(BB)
orientationtype(::Type{<:AssemblySystem{D,BB}}) where {D,BB} = orientationtype(BB)

Base.size(sys::AssemblySystem) = sys.nspcs, sys.nbonds
Base.size(sys::AssemblySystem, i) = i == 1 ? sys.nspcs : i == 2 ? sys.nbonds : 1

Base.show(io::Core.IO, sys::AssemblySystem) = print(io, "$(dimension(sys))d AssemblySytem [n=$(nspecies(sys)), k=$(nbonds(sys))]")

function bondlist(sys::AssemblySystem)
    return sys._bondlist
end
function bonded_sites(sys::AssemblySystem)
    return sys._bonded_sites
end
function bonded_species(sys::AssemblySystem)
    return sys._bonded_species
end

function siteloc2index(sys::AssemblySystem, site::BindingSiteLoc)
    return sys._siteloc2idx[site]
end
function index2siteloc(sys::AssemblySystem, index::Integer)
    return sys._idx2site[index]
end
function label2species(sys::AssemblySystem, label::Integer)
    return sys._label2species[label]
end

function isinert(sys::AssemblySystem, site::BindingSiteLoc)
    idx = siteloc2index(sys, site)
    return !any(@view intmat(assembly_system)[:, idx])
end

function anatomy(sys::AssemblySystem)
    imat = intmat(sys)
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
    for bb in buildingblocks(sys)
        a = adjacency_matrix(bb.anatomy)
        n = size(a, 1)
        A[i:i+n-1, i:i+n-1] .+= a
        i += n
    end

    vertex_labels = [zeros(Cint, ns); ones(Cint, nb)]
    anatomy = NautyDiGraph(A; vertex_labels)
    return anatomy
end

function canonid(sys::AssemblySystem)
    a = anatomy(sys)
    a_prime = NautyDiGraph(adjacency_matrix(a)'; vertex_labels=labels(a))
    return sort([canonical_id(a), canonical_id(a_prime)])
end

function _extract_nspecies(bindingrules::AbstractMatrix{<:Integer})
    _checkshape(bindingrules)
    return maximum(bindingrules[:, [1, 3]])
end

function _shift_sitelabels(buildingblocks::AbstractVector{<:BuildingBlock})
    buildingblocks = [copy(bb) for bb in buildingblocks]

    l = 1
    for bb in buildingblocks
        labs = labels(bb.anatomy)
        # TODO: make all nonidentical labels contiguous, ie. (1, 1, 3, 7) -> (1, 1, 2, 3)
        setlabels!(bb.anatomy, labs .- minimum(labs) .+ l)
        l = maximum(labels(bb.anatomy)) + 1
    end
    return buildingblocks
end
function _make_lookuptables(buildingblocks::AbstractVector{<:BuildingBlock})
    index2siteloc = BindingSiteLoc[]
    siteloc2index = Dict{BindingSiteLoc,Int}()
    label2species = Dict{Int,Int}()

    for (spcs, bb) in enumerate(buildingblocks)
        for site in 1:nsites(bb)
            spcssite = (spcs, site)
            push!(index2siteloc, spcssite)
            siteloc2index[spcssite] = length(index2siteloc)
        end
        for l in labels(anatomy(bb))
            label2species[l] = spcs
        end
    end
    return index2siteloc, siteloc2index, label2species
end

function _parse_intmat(bindingrules, site2label)
    return _intmat_from_bonds(bindingrules, site2label)
end
function _parse_intmat(bindingrules::AbstractMatrix{<:Integer}, site2label)
    _checkshape(bindingrules)
    return _intmat_from_bonds(eachrow(bindingrules), site2label)
end
function _intmat_from_bonds(bonds, site2label)
    nsites = length(site2label)
    intmat = falses(nsites, nsites)
    for bond in bonds
        spcs1, site1, spcs2, site2 = bond
        label1, label2 = site2label[spcs1, site1], site2label[spcs2, site2]
        intmat[label1, label2] = intmat[label2, label1] = true
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

# """
#     composition(p::Polyform, assembly_system::AssemblySystem)

# Return the composition of the polyform `p`, a vector containing the counts of every particle species in p
# and every bond of p.
# """
# function composition(p::Polyform, assembly_system::AssemblySystem)
#     n, k = size(assembly_system)
#     m = zeros(Int, n + k)

#     spcs = Roly.species(p)
#     for s in spcs
#         m[s] += 1
#     end

#     es = Graphs.edges(p.anatomy)
#     double_bonds = [e for e in es if reverse(e) in es]
#     bonds = []
#     for b in double_bonds
#         if b ∉ bonds && reverse(b) ∉ bonds
#             push!(bonds, b)
#         end
#     end

#     bondlist = findall(Roly.intmat(assembly_system))
#     filter!(x->x[1] <= x[2], bondlist)
#     sort!(bondlist)
    
#     for b in bonds
#         lsrc, ldst = sort!([vertex2label(p, assembly_system, b.src), vertex2label(p, assembly_system, b.dst)])
#         i = findfirst(x->x==CartesianIndex(lsrc, ldst), bondlist)
#         m[n + i] += 1
#     end

#     return m
# end


# """
#     compositions(ps::AbstractVector{<:Polyform}, sys::AssemblySystem)

# Return the compositions of all polyforms in `ps`, stacked into a matrix M, with `size(M) = (length(ps), sum(size(sys)))`.
# """
# compositions(ps::AbstractVector{<:Polyform}, sys::AssemblySystem) = reduce(vcat, composition.(ps, Ref(sys))')
