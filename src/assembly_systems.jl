"""
    AssemblySystem(bindingrules::AbstractMatrix, geometries::Vector{<:AbstractGeometry{T,F}}, face_labels=nothing) where {T,F}

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
struct AssemblySystem{D, T<:Integer, F<:AbstractFloat, G<:AbstractGeometry{T,F}}
    intmat::BitMatrix
    buildingblocks::Vector{Polyform{D,T,F}}
    geometries::Vector{G}
    n_species::Integer
    n_edges::Integer
    _sites_sum::Vector{T}
end
interactionmatrix(sys::AssemblySystem) = sys.intmat
buildingblocks(sys::AssemblySystem) = sys.buildingblocks
geometries(sys::AssemblySystem) = sys.geometries
Base.size(sys::AssemblySystem) = sys.n_species, sys.n_edges
Base.show(io::Core.IO, A::AssemblySystem{D,T,F}) where {D,T,F} = print(io, "AssemblySytem{$D,$T,$F}[n=$(A.n_species), k=$(A.n_edges)]")

function AssemblySystem(bindingrules::AbstractMatrix, geometries::Vector{<:AbstractGeometry{T,F}}, face_labels=nothing) where {T,F}
    ds = [dimension(g) for g in geometries]
    D = first(ds)
    @assert all(ds .== D)

    n_species = length(geometries)
    sites = [sum(vertices_per_site(geom)) for geom in geometries]

    monomers = Polyform{D,T,F}[]
    sites_sum = cumsum(sites)
    last_label = 0

    for i in 1:n_species
        ls = isnothing(face_labels) ? collect(1:sites[i]) : face_labels[i]
        fl = convert(Vector{T}, ls .+ last_label)
        m = create_monomer(geometries[i], T(i), fl)
        push!(monomers, m)

        last_label = fl[end]
    end

    interaction_matrix = instantiate_interactionmatrix(bindingrules, sites)
    n_edges = sum(triu(interaction_matrix))
    return AssemblySystem{D,T,F,eltype(geometries)}(interaction_matrix, monomers, geometries, n_species, n_edges, sites_sum)
end
function AssemblySystem(interactions::AbstractMatrix{<:Integer}, geometry::AbstractGeometry{T,F}, face_labels=nothing) where {T,F}
    n_species = maximum(interactions[:, [1, 3]])
    geometries = [geometry for _ in 1:n_species]
    return AssemblySystem(interactions, geometries, face_labels)
end

function spcssite2label(spcs::Integer, site::Integer, assembly_system::AssemblySystem)
    return irg_flatten(spcs, site, assembly_system._sites_sum)
end
function label2spcssite(label::Integer, assembly_system::AssemblySystem)
    return irg_unflatten(label, assembly_system._sites_sum)
end

function instantiate_interactionmatrix(interactions::AbstractMatrix{<:Integer}, sites)
    n_sites = sum(sites; init=0)
    sites_sum = cumsum(sites)

    interaction_matrix = falses(n_sites, n_sites)
    for edge in eachrow(interactions)
        a, b, c, d = edge
        i, j = irg_flatten(a, b, sites_sum), irg_flatten(c, d, sites_sum)
        interaction_matrix[i, j] = interaction_matrix[j, i] = true
    end
    return interaction_matrix
end
function instantiate_interactionmatrix(interactions::BitMatrix, args...)
    @assert interactions == interactions'
    interaction_matrix = BitMatrix(interactions)
    return interaction_matrix
end

function anatomy(asys::AssemblySystem)
    imat = interactionmatrix(asys)
    n_sites = size(imat, 1)
    n_edges = sum(triu(imat))

    A = zeros(Int, n_sites + n_edges, n_sites + n_edges)
    edge_counter = 1
    for i in axes(imat, 1), j in i:size(imat, 2)
        if imat[i, j] == 0
            continue
        end

        A[i, n_sites+edge_counter] = 1
        A[j, n_sites+edge_counter] = 1
        A[n_sites+edge_counter, i] = 1
        A[n_sites+edge_counter, j] = 1
        edge_counter += 1
    end
    
    i = 1
    for geom in asys.geometries
        a = adjacency_matrix(geom.anatomy)
        n = size(a, 1)
        A[i:i+n-1, i:i+n-1] .+= a
        i += n
    end

    vertex_labels = [zeros(Cint,n_sites); ones(Cint, n_edges)]

    anatomy = NautyDiGraph(A; vertex_labels)
    return anatomy
end

function rhash(asys::AssemblySystem)
    a = anatomy(asys)
    a_prime = NautyDiGraph(adjacency_matrix(a)'; vertex_labels=labels(a))
    return hash(sort([ghash(a), ghash(a_prime)]))
end


"""
    composition(p::Polyform, assembly_system::AssemblySystem)

Return the composition of the polyform `p`, a vector containing the counts of every particle species in p
and every bond of p.
"""
function composition(p::Polyform, assembly_system::AssemblySystem)
    n, k = size(assembly_system)
    m = zeros(Int, n + k)

    spcs = Roly.species(p)
    for s in spcs
        m[s] += 1
    end

    es = Graphs.edges(p.anatomy)
    double_bonds = [e for e in es if reverse(e) in es]
    bonds = []
    for b in double_bonds
        if b ∉ bonds && reverse(b) ∉ bonds
            push!(bonds, b)
        end
    end

    bondlist = findall(Roly.interactionmatrix(assembly_system))
    filter!(x->x[1] <= x[2], bondlist)
    sort!(bondlist)
    
    for b in bonds
        lsrc, ldst = sort!([vertex2label(p, assembly_system, b.src), vertex2label(p, assembly_system, b.dst)])
        i = findfirst(x->x==CartesianIndex(lsrc, ldst), bondlist)
        m[n + i] += 1
    end

    return m
end


"""
    compositions(ps::AbstractVector{<:Polyform}, sys::AssemblySystem)

Return the compositions of all polyforms in `ps`, stacked into a matrix M, with `size(M) = (length(ps), sum(size(sys)))`.
"""
compositions(ps::AbstractVector{<:Polyform}, sys::AssemblySystem) = reduce(vcat, composition.(ps, Ref(sys))')