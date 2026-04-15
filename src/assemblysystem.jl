const BindingSiteLoc = NTuple{2,Int}

"""
    AssemblySystem(bindingrules, particlespecies)

Defines an _assembly system_, consisting of a list of particle species and their binding rules. The structures (polyforms)
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
mutable struct AssemblySystem{D,VPS<:AbstractVector{<:ParticleSpecies}}
    intmat::Symmetric{Bool,Matrix{Bool}}
    particlespecies::VPS
    nbonds::Int
    nsites::Int
    ncolors::Int
    _color2siteloc::Dict{Int,Vector{BindingSiteLoc}}
    _siteloc2color::Dict{BindingSiteLoc,Int}
    _label2color::Dict{Int,Int}
    _bondlist::Vector{NTuple{2,Int}}
    _bonded_sites::Vector{NTuple{2,Vector{BindingSiteLoc}}}
    _bonded_species::Vector{NTuple{2,Int}}
end
function AssemblySystem(bindingrules, particlespecies::AbstractVector{<:ParticleSpecies})
    D = dimension(first(particlespecies))
    for ps in particlespecies
        dimension(ps) != D && throw(ArgumentError("particlespecies do not all have the same dimension"))
    end

    particlespecies = _shift_sitecolors_and_labels(particlespecies)
    color2siteloc, siteloc2color, label2color = _make_bindingsite_lookuptables(particlespecies)

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
    return AssemblySystem{D,typeof(particlespecies)}(intmat, particlespecies, nbonds, nsites, ncolors,
                                color2siteloc, siteloc2color, label2color, 
                                bondlist, bondedsites, bondedspecies)
end
function AssemblySystem(bindingrules, particlespecies::ParticleSpecies)
    nspcs = _extract_nspecies(bindingrules)
    particlespecies = [particlespecies for _ in 1:nspcs]
    return AssemblySystem(bindingrules, particlespecies)
end

interactionmatrix(sys::AssemblySystem) = sys.intmat
species(sys::AssemblySystem) = sys.particlespecies
species(sys::AssemblySystem, i::Integer) = sys.particlespecies[i]
nspecies(sys::AssemblySystem) = length(sys.particlespecies)
nbonds(sys::AssemblySystem) = sys.nbonds
nsites(sys::AssemblySystem) = sys.nsites
ncolors(sys::AssemblySystem) = sys.ncolors

dimension(::AssemblySystem{D}) where {D} = D
dimension(::Type{<:AssemblySystem{D}}) where {D} = D
numtype(::AssemblySystem{D,PS}) where {D,PS} = Float64 #numtype(PS)

posetype(::AssemblySystem{D,PS}) where {D,PS} = Pose{2,Float64,Angle2d{Float64}} #positiontype(BB)
speciestype(::AssemblySystem{D,PS}) where {D,PS} = Pose{2,Float64,Angle2d{Float64}} #positiontype(BB)

Base.size(sys::AssemblySystem) = sys.nspcs, sys.nbonds
Base.size(sys::AssemblySystem, i) = i == 1 ? sys.nspcs : i == 2 ? sys.nbonds : 1

Base.show(io::Core.IO, sys::AssemblySystem) = print(io, "$(dimension(sys))d AssemblySytem[n=$(nspecies(sys)), k=$(nbonds(sys))]")

@inline function bondlist(sys::AssemblySystem)
    return sys._bondlist
end
@inline function bonded_sites(sys::AssemblySystem)
    return sys._bonded_sites
end
@inline function bonded_species(sys::AssemblySystem)
    return sys._bonded_species
end

@inline function siteloc2color(sys::AssemblySystem, site::BindingSiteLoc)
    return sys._siteloc2color[site]
end
@inline function color2siteloc(sys::AssemblySystem, color::Integer)
    return sys._color2siteloc[color]
end
@inline function label2color(sys::AssemblySystem, label::Integer)
    return sys._label2color[label]
end
@inline function color2species(sys::AssemblySystem, color::Integer)
    return color2siteloc(sys, color)[1][1]
end
@inline function label2species(sys::AssemblySystem, label::Integer)
    return color2species(sys, label2color(sys, label))
end


@inline function isinert(sys::AssemblySystem, site::BindingSiteLoc)
    color = siteloc2color(sys, site)
    return isinert(sys, color)
end
@inline function isinert(sys::AssemblySystem, color::Integer)
    return !any(@view interactionmatrix(sys)[:, color])
end

function graphrep(sys::AssemblySystem)
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
        a = adjacency_matrix(bb.graphrep)
        n = size(a, 1)
        A[i:i+n-1, i:i+n-1] .+= a
        i += n
    end

    vertex_labels = [zeros(Cint, ns); ones(Cint, nb)]
    graphrep = NautyDiGraph(A; vertex_labels)
    return graphrep
end

function compatible_sitelocs(sys::AssemblySystem, color::Integer)
    intmat = interactionmatrix(sys)
    cs = 1:ncolors(sys)
    return (siteloc for c in cs if intmat[c, color] for siteloc in color2siteloc(sys, c))
end

function canonid(sys::AssemblySystem)
    a = graphrep(sys)
    a_prime = NautyDiGraph(adjacency_matrix(a)'; vertex_labels=labels(a))
    return sort([canonical_id(a), canonical_id(a_prime)])
end

function _extract_nspecies(bindingrules::AbstractMatrix{<:Integer})
    _checkshape(bindingrules)
    return maximum(bindingrules[:, [1, 3]])
end

function _shift_sitecolors_and_labels(particlespecies::AbstractVector{<:ParticleSpecies})
    particlespecies = [copy(ps) for ps in particlespecies]

    c = 1
    l = 1
    for ps in particlespecies
        g = graphrep(ps)
        labs = labels(g)

        # TODO: make all nonidentical labels contiguous, ie. (1, 1, 3, 7) -> (1, 1, 2, 3)
        setlabels!(g, labs .- minimum(labs) .+ l)
        l = maximum(labels(g)) + 1

        setcolor!(ps, c)
        c += nsites(ps)
    end
    return particlespecies
end
function _make_bindingsite_lookuptables(particlespecies::AbstractVector{<:ParticleSpecies})
    color2siteloc = Dict{Int,Vector{BindingSiteLoc}}()
    siteloc2color = Dict{BindingSiteLoc,Int}()
    label2color = Dict{Int,Int}()

    color = 1
    for (spcs, ps) in enumerate(particlespecies)
        for site_indices in equivalent_site_indices(ps)
            for site_index in site_indices
                spcssite = (spcs, site_index)
                if haskey(color2siteloc, color)
                    push!(color2siteloc[color], spcssite)
                else
                    color2siteloc[color] = [spcssite]
                end
                siteloc2color[spcssite] = color

                site = bindingsite(ps, site_index)
                for v in site.vertices
                    l = label(graphrep(ps), v)
                    l ∈ keys(label2color) && continue
                    label2color[l] = color
                end
            end
            color += 1
        end
    end
    return color2siteloc, siteloc2color, label2color
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

#     es = Graphs.edges(p.graphrep)
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
